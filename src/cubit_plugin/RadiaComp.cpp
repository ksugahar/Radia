#include "RadiaComp.hpp"
#include "ComponentPlugin.hpp"
#include "Broker.hpp"
#include "Claro.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

// Export logic (integrated, no .ccm needed)
#include "ExportGmshCommand.hpp"
#include "ExportNastranCommand.hpp"
#include "ExportVtkCommand.hpp"


#include <direct.h>  // _getcwd

#include <QAction>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLineEdit>
#include <QMessageBox>
#include <QProcess>
#include <QPushButton>
#include <QSpinBox>
#include <QHeaderView>
#include <QLabel>
#include <QTableWidget>
#include <QApplication>
#include <QClipboard>
#include <QTextStream>
#include <QDesktopServices>
#include <QUrl>
#include <QVBoxLayout>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

// Ensure DLLs in plugins/ can be found at runtime.
// With compact_netgen build, nglib/ngcore are statically linked (no DLL needed).
// This path setup is still needed for cubit_geom.dll DELAYLOAD resolution.
static void ensure_netgen_dll_path()
{
#ifdef _WIN32
  static bool done = false;
  if (done) return;
  done = true;

  const char *pd = std::getenv("CUBIT_PLUGIN_DIR");
  if (pd && pd[0]) {
    SetDllDirectoryA(pd);
    return;
  }
  // Derive plugins/ from Cubit bin/ path
  HMODULE hm = nullptr;
  GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                     GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                     (LPCSTR)&ensure_netgen_dll_path, &hm);
  if (hm) {
    char path[MAX_PATH];
    if (GetModuleFileNameA(hm, path, MAX_PATH)) {
      std::string dir(path);
      auto pos = dir.find_last_of("\\/");
      if (pos != std::string::npos) {
        dir = dir.substr(0, pos);
        std::string plugins = dir + "\\plugins";
        DWORD attr = GetFileAttributesA(plugins.c_str());
        if (attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY))
          dir = plugins;
        SetDllDirectoryA(dir.c_str());
      }
    }
  }
#endif
}

// ============================================================
// COMPONENT_PLUGIN exports (Claro framework .ccl)
// ============================================================
COMPONENT_PLUGIN

void register_components()
{
  Broker::instance()->register_component(new RadiaComp());
}

void print_plugin_version()
{
  Broker::command << "Radia Cubit Plugin 1.0\n";
}

void print_plugin_options()
{
}

// ============================================================
// RadiaComp - Component lifecycle
// ============================================================

RadiaComp::RadiaComp()
  : Component("radiacomp"), mMenuInitialized(false)
{}

RadiaComp::~RadiaComp() {}

void RadiaComp::start_up(int withGUI)
{
  if (withGUI)
    setup_menus();
}

void RadiaComp::clean_up()
{
  cleanup_menus();
  clean_up_complete();
}

void RadiaComp::setup_menus()
{
  if (mMenuInitialized)
    return;

  Claro* gui = Claro::instance();
  if (!gui)
    return;

  RadiaMenuHandler* handler = new RadiaMenuHandler();
  std::vector<QAction*> menu_list;

  QAction* a_gmsh = new QAction("GMSH...", handler);
  a_gmsh->setStatusTip("Export mesh to GMSH format (.msh)");
  // Menu order follows lab priority: .vol > .msh > .bdf > .vtk
  QAction* a_netgen = new QAction("Netgen Vol (.vol)...", handler);
  a_netgen->setStatusTip("Export curved mesh with labels (.vol, order 1-5)");
  QObject::connect(a_netgen, SIGNAL(triggered()), handler, SLOT(export_netgen()));
  menu_list.push_back(a_netgen);

  QObject::connect(a_gmsh, SIGNAL(triggered()), handler, SLOT(export_gmsh()));
  menu_list.push_back(a_gmsh);

  QAction* a_nastran = new QAction("Nastran BDF...", handler);
  a_nastran->setStatusTip("Export mesh to Nastran BDF format (.bdf)");
  QObject::connect(a_nastran, SIGNAL(triggered()), handler, SLOT(export_nastran()));
  menu_list.push_back(a_nastran);

  QAction* a_vtk = new QAction("VTK...", handler);
  a_vtk->setStatusTip("Export mesh to VTK format (.vtk)");
  QObject::connect(a_vtk, SIGNAL(triggered()), handler, SLOT(export_vtk()));
  menu_list.push_back(a_vtk);

  // Separator before evaluation tools
  QAction* sep = new QAction(handler);
  sep->setSeparator(true);
  menu_list.push_back(sep);

  QAction* a_vol = new QAction("Mesh Evaluation...", handler);
  a_vol->setStatusTip("Volume/Area p-convergence check (order 1-5) against CAD");
  QObject::connect(a_vol, SIGNAL(triggered()), handler, SLOT(mesh_volume()));
  menu_list.push_back(a_vol);

  gui->add_to_menu("&Export Mesh", menu_list, "radiacomp");

  mMenuInitialized = true;
}

void RadiaComp::cleanup_menus()
{
  if (!mMenuInitialized)
    return;

  Claro* gui = Claro::instance();
  if (gui)
    gui->remove_menu_items("radiacomp");

  mMenuInitialized = false;
}

// ============================================================
// RadiaMenuHandler - show dialog and execute command
// ============================================================

//! Check if Cubit has geometry. If not, prompt for a .jou file.
//! Returns the journal path if one was played (empty if model was already loaded).
static QString ensure_model()
{
  if (CubitInterface::get_volume_count() > 0)
    return QString();

  // No geometry — ask user to load a journal file
  QString jou = QFileDialog::getOpenFileName(
      nullptr, "No model loaded - Select Journal File",
      QString(), "Cubit Journal (*.jou);;All Files (*)");
  if (jou.isEmpty())
    return QString();

  jou.replace("\\", "/");
  std::string cmd = "play \"" + jou.toStdString() + "\"";
  CubitInterface::cmd(cmd.c_str());

  if (CubitInterface::get_volume_count() > 0)
    return jou;
  return QString();
}

// Persistent journal path across exports in one session
static QString s_lastJouPath;

static void run_export(ExportDialog::Format fmt)
{
  if (CubitInterface::get_volume_count() == 0) {
    QString jou = ensure_model();
    if (CubitInterface::get_volume_count() == 0)
      return;
    if (!jou.isEmpty())
      s_lastJouPath = jou;
  }

  // Try get_current_journal_file() first, fall back to saved path
  QString jouPath;
  std::string jf = CubitInterface::get_current_journal_file();
  if (!jf.empty())
    jouPath = QString::fromStdString(jf);
  else
    jouPath = s_lastJouPath;

  // Modal dialog — Cubit command execution requires the GUI thread,
  // so the dialog must block until the user clicks OK/Cancel.
  ExportDialog dlg(fmt, jouPath);
  if (dlg.exec() != QDialog::Accepted)
    return;

  // Set Cubit working directory to output file's directory
  // (prevents intermediate files landing in repo root)
  QString outDir = QFileInfo(dlg.filePath()).absolutePath();
  outDir.replace("\\", "/");
  CubitInterface::silent_cmd(("cd \"" + outDir.toLocal8Bit() + "\"").constData());

  // Execute export via APREPRO command (same codepath as journal playback)
  std::string cmd = dlg.cubitCommand().toStdString();
  if (cmd.empty()) {
    PRINT_ERROR("No export command.\n");
    return;
  }

  CubitInterface::cmd(cmd.c_str());

  QString outFile = dlg.filePath();
  if (!QFile::exists(outFile)) {
    PRINT_ERROR("Export failed: %s\n", outFile.toStdString().c_str());
    return;
  }
  PRINT_INFO("Export complete: %s\n", outFile.toStdString().c_str());

  // Open exported file with OS-associated application
  QDesktopServices::openUrl(QUrl::fromLocalFile(outFile));
}

void RadiaMenuHandler::export_gmsh()    { run_export(ExportDialog::GMSH); }
void RadiaMenuHandler::export_nastran() { run_export(ExportDialog::Nastran); }
void RadiaMenuHandler::export_vtk()     { run_export(ExportDialog::VTK); }
// ============================================================
// Helpers for subprocess-based operations (Netgen export, Mesh Eval)
// ============================================================
static QString find_external_python();
static QString find_calc_script(const QString &python, const QString &name);

void RadiaMenuHandler::export_netgen()
{
  if (CubitInterface::get_volume_count() == 0) {
    ensure_model();
    if (CubitInterface::get_volume_count() == 0) return;
  }

  // Modal dialog — Cubit commands require GUI thread
  QString jouPath;
  std::string jf = CubitInterface::get_current_journal_file();
  if (!jf.empty()) jouPath = QString::fromStdString(jf);
  if (jouPath.isEmpty()) jouPath = s_lastJouPath;

  ExportDialog dlgOpt(ExportDialog::NETGEN_VOL, jouPath);
  if (dlgOpt.exec() != QDialog::Accepted) return;

  int order = dlgOpt.order();
  QString volPath = dlgOpt.filePath();
  volPath.replace("\\", "/");

  // Ensure Netgen DLLs are on search path before export netgen command
  ensure_netgen_dll_path();

  // Set Cubit working directory to output file's directory
  QString outDir = QFileInfo(volPath).absolutePath();
  outDir.replace("\\", "/");
  CubitInterface::silent_cmd(("cd \"" + outDir.toLocal8Bit() + "\"").constData());

  // --- Phase 1: C++ export .vol + companion JSON (no subprocess, no cub5) ---
  PRINT_INFO("Exporting Netgen .vol (order %d)...\n", order);
  std::string cmd = std::string("export netgen \"")
    + volPath.toLocal8Bit().constData()
    + "\" order " + std::to_string(order) + " overwrite";
  CubitInterface::silent_cmd(cmd.c_str());

  // Check .vol was created
  if (!QFile::exists(volPath)) {
    PRINT_ERROR("export netgen command failed.\n");
    return;
  }
  PRINT_INFO("Exported: %s\n", volPath.toLocal8Bit().constData());

  // Open exported .vol with OS-associated application (vol-viewer)
  QDesktopServices::openUrl(QUrl::fromLocalFile(volPath));

  // --- Phase 2: subprocess to verify .vol with NGSolve ---
  QString python = find_external_python();
  QString script = find_calc_script(python, "calc_verify_vol");
  if (script.isEmpty()) {
    PRINT_WARNING("calc_verify_vol.py not found. Skipping NGSolve verification.\n");
    return;
  }

  PRINT_INFO("Verifying mesh with NGSolve...\n");

  QProcess proc;
  QStringList args;
  if (python == "py") args << "-3";
  args << script << "--vol" << volPath;

  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  for (const QString &key : env.keys()) {
    if (key.toUpper().contains("QT") || key.toUpper().contains("PYSIDE"))
      env.remove(key);
  }
  proc.setProcessEnvironment(env);

  proc.start(python, args);
  if (!proc.waitForFinished(300000)) {
    PRINT_WARNING("NGSolve verification timed out.\n");
    return;
  }

  // Capture output before consuming
  QByteArray stdoutData = proc.readAllStandardOutput();
  QByteArray stderrData = proc.readAllStandardError();

  // Write log file for debugging
  {
    QString logDir = "S:/Radia/01_GitHub/logs";
    if (QDir(logDir).exists()) {
      QString ts = QDateTime::currentDateTime().toString("yyyyMMdd_HHmmss");
      QString logFile = QDir(logDir).filePath(
          "export_vol_" + ts + ".log");
      QFile lf(logFile);
      if (lf.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream ls(&lf);
        ls << "=== export netgen .vol log ===\n";
        ls << "time: " << QDateTime::currentDateTime().toString(Qt::ISODate) << "\n";
        ls << "vol: " << volPath << "\n";
        ls << "order: " << order << "\n";
        ls << "exit_code: " << proc.exitCode() << "\n\n";
        ls << "--- stdout ---\n" << stdoutData << "\n";
        ls << "--- stderr ---\n" << stderrData << "\n";
        lf.close();
      }
    }
  }

  if (proc.exitCode() != 0) {
    PRINT_WARNING("NGSolve verification failed:\n%s\n",
      stderrData.left(2000).constData());
    return;
  }

  // Parse JSON result from subprocess
  QString out = QString::fromUtf8(stdoutData);
  QStringList lines = out.split('\n', Qt::SkipEmptyParts);
  if (lines.isEmpty()) { PRINT_WARNING("No verification output.\n"); return; }

  QJsonDocument doc = QJsonDocument::fromJson(lines.last().toUtf8());
  if (!doc.isObject()) { PRINT_WARNING("Invalid verification JSON.\n"); return; }

  QJsonObject r = doc.object();
  if (r.contains("error")) {
    PRINT_WARNING("Verification error: %s\n",
      r["error"].toString().toStdString().c_str());
    return;
  }

  // --- Show consistency check dialog (modeless: user can interact with Cubit) ---
  QJsonArray mats = r["materials"].toArray();
  QJsonArray bnds = r["boundaries"].toArray();
  QJsonArray warns = r["warnings"].toArray();
  int nMats = mats.size();
  int nBnds = bnds.size();

  QDialog *dlg = new QDialog();
  dlg->setAttribute(Qt::WA_DeleteOnClose);
  dlg->setWindowTitle(QString("Netgen Vol Export - Order %1").arg(order));
  dlg->setMinimumSize(700, 400);
  QVBoxLayout *layout = new QVBoxLayout(dlg);

  // File info
  layout->addWidget(new QLabel(
    QString("Output: %1\nElements: %2, Order: %3")
      .arg(volPath).arg(r["n_elements"].toInt()).arg(order)));

  // Helper lambda for table creation
  auto makeTable = [&](const char *label, const QStringList &headers,
                       const QJsonArray &data,
                       std::function<void(QTableWidget*, int, const QJsonObject&)> fillRow) {
    layout->addWidget(new QLabel(label));
    int n = data.size();
    if (n == 0) {
      layout->addWidget(new QLabel("  (none)"));
      return;
    }
    QTableWidget *tbl = new QTableWidget(n, headers.size(), dlg);
    tbl->setHorizontalHeaderLabels(headers);
    tbl->horizontalHeader()->setStretchLastSection(true);
    tbl->verticalHeader()->setVisible(false);
    tbl->setEditTriggers(QAbstractItemView::NoEditTriggers);
    for (int i = 0; i < n; i++)
      fillRow(tbl, i, data[i].toObject());
    tbl->resizeColumnsToContents();
    layout->addWidget(tbl);
  };

  auto errItem = [](double err) {
    return new QTableWidgetItem(QString::number(err, 'e', 2));
  };

  // Material table (Volume)
  makeTable("Materials (Volume):",
    {"Name", "CAD Volume", "NGSolve Volume", "Error [%]"},
    mats, [&](QTableWidget *t, int i, const QJsonObject &m) {
      t->setItem(i, 0, new QTableWidgetItem(m["name"].toString()));
      t->setItem(i, 1, new QTableWidgetItem(
        QString::number(m["cad_volume"].toDouble(), 'e', 6)));
      t->setItem(i, 2, new QTableWidgetItem(
        m.contains("ng_volume") ?
          QString::number(m["ng_volume"].toDouble(), 'e', 6) : "N/A"));
      if (m.contains("error_pct"))
        t->setItem(i, 3, errItem(m["error_pct"].toDouble()));
    });

  // Boundary table (Area)
  makeTable("Boundaries (Surface Area):",
    {"Name", "CAD Area", "NGSolve Area", "Error [%]"},
    bnds, [&](QTableWidget *t, int i, const QJsonObject &b) {
      t->setItem(i, 0, new QTableWidgetItem(b["name"].toString()));
      t->setItem(i, 1, new QTableWidgetItem(
        QString::number(b["cad_area"].toDouble(), 'e', 6)));
      t->setItem(i, 2, new QTableWidgetItem(
        b.contains("ng_area") ?
          QString::number(b["ng_area"].toDouble(), 'e', 6) : "N/A"));
      if (b.contains("error_pct"))
        t->setItem(i, 3, errItem(b["error_pct"].toDouble()));
    });

  // Edge table (Length / BBND) — always shown, even if empty
  QJsonArray edges = r["edges"].toArray();
  makeTable("Edges (BBND Length):",
    {"Name", "CAD Length", "NGSolve Length", "Error [%]"},
    edges, [&](QTableWidget *t, int i, const QJsonObject &e) {
      t->setItem(i, 0, new QTableWidgetItem(e["name"].toString()));
      t->setItem(i, 1, new QTableWidgetItem(
        QString::number(e["cad_length"].toDouble(), 'e', 6)));
      t->setItem(i, 2, new QTableWidgetItem(
        e.contains("ng_length") ?
          QString::number(e["ng_length"].toDouble(), 'e', 6) : "N/A"));
      if (e.contains("error_pct"))
        t->setItem(i, 3, errItem(e["error_pct"].toDouble()));
    });

  // Warnings
  if (!warns.isEmpty()) {
    QString warnText = "<b style='color:red'>Warnings:</b><ul>";
    for (auto w : warns)
      warnText += "<li>" + w.toString() + "</li>";
    warnText += "</ul>";
    QLabel *warnLabel = new QLabel(warnText);
    warnLabel->setTextFormat(Qt::RichText);
    layout->addWidget(warnLabel);
  }

  // Buttons
  QDialogButtonBox *buttons = new QDialogButtonBox(
    QDialogButtonBox::Ok, dlg);
  layout->addWidget(buttons);
  QObject::connect(buttons, &QDialogButtonBox::accepted, dlg, &QDialog::close);

  PRINT_INFO("Exported: %s (order %d)\n",
    volPath.toStdString().c_str(), order);

  dlg->show();  // modeless result dialog: user can interact with Cubit
}

// ============================================================
// Mesh Evaluation — subprocess to calc_mesh_eval.py (p=1..5)
// ============================================================

// Find external Python 3.12 (not Cubit's 3.10)
static QString find_external_python()
{
  QByteArray env = qgetenv("RADIA_PYTHON");
  if (!env.isEmpty()) return QString::fromLocal8Bit(env);
#ifdef _WIN32
  QProcess probe;
  probe.start("py", {"-3", "-c", "print('ok')"});
  if (probe.waitForFinished(5000) && probe.exitCode() == 0) return "py";
#endif
  return "python";
}

// Find calc_mesh_eval.py without importing radia
static QString find_calc_script(const QString &python, const QString &name)
{
  QProcess p;
  QStringList args;
  if (python == "py") args << "-3";
  args << "-c" << QString("import radia.panels.%1 as m; import os; "
                           "print(os.path.abspath(m.__file__))").arg(name);
  p.start(python, args);
  if (p.waitForFinished(10000) && p.exitCode() == 0)
    return QString::fromUtf8(p.readAllStandardOutput()).trimmed();
  return QString();
}

void RadiaMenuHandler::mesh_volume()
{
  if (CubitInterface::get_volume_count() == 0) {
    ensure_model();
    if (CubitInterface::get_volume_count() == 0) return;
  }

  // Ensure Netgen DLLs are on search path before export netgen commands
  ensure_netgen_dll_path();

  // --- Phase 1: Save temp .cub5 for subprocess ---
  QString tmpDir = QDir::tempPath();
  tmpDir.replace("\\", "/");
  if (!tmpDir.endsWith('/')) tmpDir += '/';
  QString cub5Path = tmpDir + "radia_eval_temp.cub5";
  {
    std::string cmd = "save cub5 \"" + cub5Path.toStdString() + "\" overwrite";
    CubitInterface::silent_cmd(cmd.c_str());
    if (!QFile::exists(cub5Path)) {
      PRINT_ERROR("Failed to save temp .cub5\n");
      return;
    }
  }

  // --- Phase 2: subprocess — calc_mesh_eval.py does extract_curved_mesh p=1..5 ---
  QString python = find_external_python();
  QString script = find_calc_script(python, "calc_mesh_eval");
  if (script.isEmpty()) {
    PRINT_ERROR("Cannot find calc_mesh_eval.py. Is radia installed?\n");
    QFile::remove(cub5Path);
    return;
  }

  PRINT_INFO("Evaluating mesh p=1..5 with NGSolve...\n");

  QProcess proc;
  QStringList args;
  if (python == "py") args << "-3";
  args << script << "--cub5" << cub5Path;

  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  for (const QString &key : env.keys()) {
    if (key.toUpper().contains("QT") || key.toUpper().contains("PYSIDE"))
      env.remove(key);
  }
  proc.setProcessEnvironment(env);

  proc.start(python, args);
  if (!proc.waitForFinished(600000)) {
    PRINT_ERROR("Mesh evaluation timed out.\n");
    return;
  }
  if (proc.exitCode() != 0) {
    PRINT_ERROR("Mesh evaluation failed:\n%s\n",
      proc.readAllStandardError().left(2000).constData());
    return;
  }

  // Cleanup temp .cub5
  QFile::remove(cub5Path);

  // Parse JSON from last line
  QString out = QString::fromUtf8(proc.readAllStandardOutput());
  QStringList lines = out.split('\n', Qt::SkipEmptyParts);
  if (lines.isEmpty()) { PRINT_ERROR("No output.\n"); return; }

  QJsonDocument doc = QJsonDocument::fromJson(lines.last().toUtf8());
  if (!doc.isObject()) { PRINT_ERROR("Invalid JSON.\n"); return; }

  QJsonObject r = doc.object();
  if (r.contains("error")) {
    PRINT_ERROR("Error: %s\n", r["error"].toString().toStdString().c_str());
    return;
  }

  // Build result table
  double cad_v = r["cad_vol_total"].toDouble();
  double cad_a = r["cad_area_total"].toDouble();
  double cad_l = r["cad_length_total"].toDouble();

  QJsonArray orders = r["orders"].toArray();
  int nrows = orders.size();

  // --- Dialog with QTableWidget ---
  QDialog dlg;
  dlg.setWindowTitle("Mesh Evaluation - Volume / Area / Length");
  QVBoxLayout *layout = new QVBoxLayout(&dlg);

  // CAD reference labels
  QLabel *cadLabel = new QLabel(
    QString("CAD Volume: %1 m^3     CAD Area: %2 m^2     CAD Length: %3 m")
      .arg(cad_v, 0, 'e', 6).arg(cad_a, 0, 'e', 6).arg(cad_l, 0, 'e', 6), &dlg);
  layout->addWidget(cadLabel);

  // Table: empty column after L err for visual separation from future columns
  // Helper: size table to fit all contents (no scrollbar)
  auto fitTable = [](QTableWidget *t) {
    t->resizeColumnsToContents();
    t->resizeRowsToContents();
    int w = t->verticalHeader()->isVisible() ? t->verticalHeader()->width() : 0;
    for (int c = 0; c < t->columnCount(); c++)
      w += t->columnWidth(c);
    w += t->frameWidth() * 2 + 4;
    int h = t->horizontalHeader()->height();
    for (int r = 0; r < t->rowCount(); r++)
      h += t->rowHeight(r);
    h += t->frameWidth() * 2 + 2;
    t->setFixedSize(w, h);
  };

  QTableWidget *table = new QTableWidget(nrows, 7, &dlg);
  table->setHorizontalHeaderLabels(
    {"p", "Volume", "V err [%]", "Area", "A err [%]", "Length", "L err [%]"});
  table->verticalHeader()->setVisible(false);
  table->setSelectionMode(QAbstractItemView::ContiguousSelection);
  table->setEditTriggers(QAbstractItemView::NoEditTriggers);

  for (int i = 0; i < nrows; i++) {
    QJsonObject o = orders[i].toObject();
    int p = o["order"].toInt();
    table->setItem(i, 0, new QTableWidgetItem(QString::number(p)));

    if (o.contains("error")) {
      table->setItem(i, 1, new QTableWidgetItem(o["error"].toString()));
    } else {
      table->setItem(i, 1, new QTableWidgetItem(
        QString::number(o["ng_volume"].toDouble(), 'e', 6)));
      table->setItem(i, 2, new QTableWidgetItem(
        QString::number(o["vol_error_pct"].toDouble(), 'e', 2)));
      table->setItem(i, 3, new QTableWidgetItem(
        QString::number(o["ng_area"].toDouble(), 'e', 6)));
      table->setItem(i, 4, new QTableWidgetItem(
        QString::number(o["area_error_pct"].toDouble(), 'e', 2)));
      table->setItem(i, 5, new QTableWidgetItem(
        QString::number(o["ng_length"].toDouble(), 'e', 6)));
      table->setItem(i, 6, new QTableWidgetItem(
        QString::number(o["len_error_pct"].toDouble(), 'e', 2)));
    }
  }
  fitTable(table);
  layout->addWidget(table);

  // --- Format round-trip table (GMSH API verification) ---
  QJsonArray fmtArr = r["format_roundtrip"].toArray();
  int nfmt = fmtArr.size();
  if (nfmt > 0) {
    layout->addWidget(new QLabel(
        "\nFormat Round-Trip vs .vol (GMSH API Jacobian verification):"));
    QTableWidget *fmtTable = new QTableWidget(nfmt, 7, &dlg);
    fmtTable->setHorizontalHeaderLabels(
        {"Format", "Order", "Volume", "V err [%]", "Area", "A err [%]", "neg_det"});
    fmtTable->verticalHeader()->setVisible(false);
    fmtTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    for (int i = 0; i < nfmt; i++) {
      QJsonObject fo = fmtArr[i].toObject();
      fmtTable->setItem(i, 0, new QTableWidgetItem(fo["format"].toString()));
      fmtTable->setItem(i, 1, new QTableWidgetItem(
          QString::number(fo["order"].toInt())));
      if (fo.contains("error")) {
        auto *errItem = new QTableWidgetItem(fo["error"].toString());
        errItem->setBackground(QColor(255, 200, 200));
        fmtTable->setItem(i, 2, errItem);
      } else {
        fmtTable->setItem(i, 2, new QTableWidgetItem(
            QString::number(fo["volume"].toDouble(), 'e', 6)));
        fmtTable->setItem(i, 3, new QTableWidgetItem(
            QString::number(fo["vol_error_pct"].toDouble(), 'e', 2)));
        fmtTable->setItem(i, 4, new QTableWidgetItem(
            QString::number(fo["area"].toDouble(), 'e', 6)));
        fmtTable->setItem(i, 5, new QTableWidgetItem(
            QString::number(fo["area_error_pct"].toDouble(), 'e', 2)));
        fmtTable->setItem(i, 6, new QTableWidgetItem(
            QString::number(fo["neg_det"].toInt())));
      }
    }
    fitTable(fmtTable);
    layout->addWidget(fmtTable);
  }

  dlg.setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  layout->setSizeConstraint(QLayout::SetFixedSize);

  // Buttons: OK + Copy
  QDialogButtonBox *buttons = new QDialogButtonBox(
    QDialogButtonBox::Ok | QDialogButtonBox::Save, &dlg);
  buttons->button(QDialogButtonBox::Save)->setText("Copy to Clipboard");
  layout->addWidget(buttons);

  QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
  QObject::connect(buttons->button(QDialogButtonBox::Save), &QPushButton::clicked,
    [&, table, nrows, fmtArr, nfmt]() {
      QString tsv = "p\tVolume\tV err [%]\tArea\tA err [%]\tLength\tL err [%]\n";
      for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < 7; j++) {
          if (j > 0) tsv += "\t";
          QTableWidgetItem *item = table->item(i, j);
          tsv += item ? item->text() : "";
        }
        tsv += "\n";
      }
      if (nfmt > 0) {
        tsv += "\nFormat\tOrder\tVolume\tV err [%]\tArea\tA err [%]\tneg_det\n";
        for (int i = 0; i < nfmt; i++) {
          QJsonObject fo = fmtArr[i].toObject();
          tsv += fo["format"].toString() + "\t";
          tsv += QString::number(fo["order"].toInt()) + "\t";
          if (fo.contains("error")) {
            tsv += fo["error"].toString() + "\t\t\t\t\n";
          } else {
            tsv += QString::number(fo["volume"].toDouble(), 'e', 6) + "\t";
            tsv += QString::number(fo["vol_error_pct"].toDouble(), 'e', 2) + "\t";
            tsv += QString::number(fo["area"].toDouble(), 'e', 6) + "\t";
            tsv += QString::number(fo["area_error_pct"].toDouble(), 'e', 2) + "\t";
            tsv += QString::number(fo["neg_det"].toInt()) + "\n";
          }
        }
      }
      QApplication::clipboard()->setText(tsv);
    });

  dlg.exec();
}

// ============================================================
// ExportDialog - format-specific options
// ============================================================

// --- Settings persistence ---
static QString settingsDir()
{
  // Use AppData/Roaming/Radia/ (not ~/.cubit which is a file on this system)
  QString appdata = QDir::homePath() + "/AppData/Roaming/Radia";
  if (!QDir(appdata).exists())
    QDir().mkpath(appdata);
  return appdata;
}

static QString settingsPath()
{
  return settingsDir() + "/export_settings.json";
}

static QJsonObject loadSettings()
{
  QFile f(settingsPath());
  if (!f.open(QIODevice::ReadOnly)) return {};
  return QJsonDocument::fromJson(f.readAll()).object();
}

static void saveSettings(const QJsonObject &obj)
{
  QFile f(settingsPath());
  if (f.open(QIODevice::WriteOnly))
    f.write(QJsonDocument(obj).toJson(QJsonDocument::Indented));
}

ExportDialog::ExportDialog(Format format, const QString &jouPath, QWidget* parent)
  : QDialog(parent), mFormat(format),
    mVersion(nullptr), mNoPyramid(nullptr)
{
  // Window title
  const char* titles[] = {"Export Netgen Vol", "Export GMSH", "Export Nastran BDF", "Export VTK"};
  setWindowTitle(titles[format]);
  setMinimumWidth(500);

  QVBoxLayout* layout = new QVBoxLayout(this);
  QFormLayout* form = new QFormLayout();

  // Determine default directory and base name
  const char* exts[] = {".vol", ".msh", ".bdf", ".vtk"};
  QString defaultDir;
  QString baseName = "ExportedMesh";

  if (!jouPath.isEmpty()) {
    // Use journal file's directory and basename
    QString jp = jouPath;
    jp.replace("\\", "/");
    int lastSlash = jp.lastIndexOf('/');
    if (lastSlash >= 0)
      defaultDir = jp.left(lastSlash);
    // Extract basename without extension
    QString fname = (lastSlash >= 0) ? jp.mid(lastSlash + 1) : jp;
    int dot = fname.lastIndexOf('.');
    if (dot > 0)
      baseName = fname.left(dot);
  }

  if (defaultDir.isEmpty()) {
    // Try default_dir from export_settings.json (set by Python startup)
    QJsonObject all = loadSettings();
    if (all.contains("default_dir")) {
      QString dd = all["default_dir"].toString();
      if (QDir(dd).exists())
        defaultDir = dd;
    }
    // Final fallback: CWD
    if (defaultDir.isEmpty()) {
      char cwd[1024];
      if (_getcwd(cwd, sizeof(cwd))) {
        defaultDir = QString::fromLocal8Bit(cwd);
        defaultDir.replace("\\", "/");
      }
    }
  }

  // Directory row: text + browse button
  QHBoxLayout* dirRow = new QHBoxLayout();
  mDir = new QLineEdit(defaultDir);
  QPushButton* browseBtn = new QPushButton("...");
  browseBtn->setFixedWidth(30);
  connect(browseBtn, SIGNAL(clicked()), this, SLOT(browseDir()));
  connect(mDir, SIGNAL(textChanged(QString)), this, SLOT(updatePreview()));
  dirRow->addWidget(mDir);
  dirRow->addWidget(browseBtn);
  form->addRow("Directory:", dirRow);

  // Filename row: text input (no browse, user types)
  mFileName = new QLineEdit(baseName + exts[format]);
  connect(mFileName, SIGNAL(textChanged(QString)), this, SLOT(updatePreview()));
  form->addRow("Filename:", mFileName);

  // Order
  mOrderCombo = new QComboBox();
  if (format == NETGEN_VOL) {
    mOrderCombo->addItems({"1", "2", "3", "4", "5"});
    mOrderCombo->setCurrentIndex(2);  // default order 3
  } else {
    mOrderCombo->addItems({"1 (linear)", "2 (quadratic)"});
  }
  connect(mOrderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
  form->addRow("Order:", mOrderCombo);

  // Version: GMSH is always v2.2 (v4.1 is for post-processing only)

  // Dimension (not for NETGEN_VOL)
  mDimension = new QComboBox();
  mDimension->addItems({"3D", "2D"});
  connect(mDimension, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
  if (format != NETGEN_VOL) {
    form->addRow("Dimension:", mDimension);
  }

  // NoPyramid (Nastran only)
  if (format == Nastran) {
    mNoPyramid = new QComboBox();
    mNoPyramid->addItems({"Keep pyramids", "Convert to degenerate hex (JMAG)"});
    connect(mNoPyramid, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
    form->addRow("Pyramids:", mNoPyramid);
  }

  layout->addLayout(form);

  // Command preview: left-aligned, click to select all for easy copy
  QFormLayout* previewForm = new QFormLayout();
  mPreview = new QLineEdit();
  mPreview->setReadOnly(true);
  mPreview->setAlignment(Qt::AlignLeft);
  connect(mPreview, &QLineEdit::cursorPositionChanged, [this](int, int) {
    if (!mPreview->hasSelectedText())
      mPreview->setCursorPosition(0);
  });
  // Focus = select all (click or tab)
  auto *prevPtr = mPreview;  // capture for lambda
  connect(qApp, &QApplication::focusChanged, [prevPtr](QWidget*, QWidget *now) {
    if (now == prevPtr) prevPtr->selectAll();
  });
  previewForm->addRow("Command:", mPreview);
  layout->addLayout(previewForm);

  // OK / Cancel
  QDialogButtonBox* buttons = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
  connect(buttons, &QDialogButtonBox::accepted, [this]() {
    // Save settings before accepting
    const char* keys[] = {"netgen_vol", "gmsh", "nastran", "vtk"};
    QJsonObject all = loadSettings();
    QJsonObject s;
    s["dir"] = mDir->text();
    s["order"] = mOrderCombo->currentIndex();
    if (mDimension) s["dimension"] = mDimension->currentIndex();
    if (mNoPyramid) s["nopyramid"] = mNoPyramid->currentIndex();
    all[keys[mFormat]] = s;
    saveSettings(all);
    accept();
  });
  layout->addWidget(buttons);

  // Load saved settings (override defaults with previous values)
  {
    const char* keys[] = {"netgen_vol", "gmsh", "nastran", "vtk"};
    QJsonObject all = loadSettings();
    QJsonObject s = all[keys[format]].toObject();
    if (s.contains("dir")) mDir->setText(s["dir"].toString());
    if (s.contains("order")) mOrderCombo->setCurrentIndex(s["order"].toInt());
    if (s.contains("dimension") && mDimension)
      mDimension->setCurrentIndex(s["dimension"].toInt());
    if (s.contains("nopyramid") && mNoPyramid)
      mNoPyramid->setCurrentIndex(s["nopyramid"].toInt());
  }

  updatePreview();
}

void ExportDialog::browseDir()
{
  QString dir = QFileDialog::getExistingDirectory(this, "Select Output Directory",
                                                   mDir->text());
  if (!dir.isEmpty()) {
    dir.replace("\\", "/");
    mDir->setText(dir);
  }
}

QString ExportDialog::filePath() const
{
  QString dir = mDir->text();
  dir.replace("\\", "/");
  if (!dir.endsWith('/')) dir += '/';
  return dir + mFileName->text();
}

QString ExportDialog::cubitCommand() const
{
  QString file = filePath();
  int order = mOrderCombo->currentIndex() + 1;

  QString cmd;
  switch (mFormat) {
    case NETGEN_VOL: {
      // Not a Cubit APREPRO command — shown as preview only
      cmd = QString("export netgen \"%1\" order %2").arg(file).arg(order);
      break;
    }
    case GMSH: {
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("export gmsh \"%1\" order %2 version 2 dimension %3")
                .arg(file).arg(order).arg(dim);
      break;
    }
    case Nastran: {
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("export radia_nastran \"%1\" order %2 dimension %3")
                .arg(file).arg(order).arg(dim);
      if (mNoPyramid && mNoPyramid->currentIndex() == 1)
        cmd += " nopyramid";
      break;
    }
    case VTK: {
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("export vtk \"%1\" order %2 dimension %3")
                .arg(file).arg(order).arg(dim);
      break;
    }
  }
  cmd += " overwrite";
  return cmd;
}

void ExportDialog::updatePreview()
{
  mPreview->setText(cubitCommand());
  mPreview->setCursorPosition(0);  // show beginning of command
}
