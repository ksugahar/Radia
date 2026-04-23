#include "RadiaComp.hpp"
#include "ComponentPlugin.hpp"
#include "Broker.hpp"
#include "Claro.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include "LauncherLogic.hpp"

// Export logic (integrated, no .ccm needed)
#include "ExportGmshCommand.hpp"
#include "ExportNastranCommand.hpp"
#include "ExportVtkCommand.hpp"


#include <direct.h>  // _getcwd

#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLineEdit>
#include <QMessageBox>
#include <QProcess>
#include <QRegularExpression>
#include <QPushButton>
#include <QSet>
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

// ============================================================
// CCL debug log — writes to C:\radia_ccl_log.txt
// Append-mode, timestamped. Safe to leave on permanently.
// ============================================================
static void ccl_log(const char *fmt, ...)
{
  FILE *f = fopen("C:\\radia_ccl_log.txt", "a");
  if (!f) return;
  // Timestamp + username
  QDateTime now = QDateTime::currentDateTime();
  QString user = qEnvironmentVariable("USERNAME", "?");
  fprintf(f, "[%s %s] ",
          now.toString("yyyy-MM-dd HH:mm:ss").toLocal8Bit().constData(),
          user.toLocal8Bit().constData());
  va_list ap;
  va_start(ap, fmt);
  vfprintf(f, fmt, ap);
  va_end(ap);
  fprintf(f, "\n");
  fclose(f);
}

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
  ccl_log("=== Radia Cubit Plugin start_up (GUI=%d) ===", withGUI);
  ccl_log("  version: 1.0 (built " __DATE__ " " __TIME__ ")");
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

  QAction* a_femeem = new QAction("FEMEEM...", handler);
  a_femeem->setStatusTip("Export mesh to FEMEEM format (1st-order tet, directory output)");
  QObject::connect(a_femeem, SIGNAL(triggered()), handler, SLOT(export_femeem()));
  menu_list.push_back(a_femeem);

  QAction* a_meg = new QAction("MEG (ELF/MAGIC)...", handler);
  a_meg->setStatusTip("Export mesh to ELF/MAGIC MEG format (1st-order)");
  QObject::connect(a_meg, SIGNAL(triggered()), handler, SLOT(export_meg()));
  menu_list.push_back(a_meg);

  // Separator before evaluation tools
  QAction* sep = new QAction(handler);
  sep->setSeparator(true);
  menu_list.push_back(sep);

  QAction* a_vol = new QAction("Mesh Evaluation...", handler);
  a_vol->setStatusTip("Volume/Area p-convergence check (order 1-5) against CAD");
  QObject::connect(a_vol, SIGNAL(triggered()), handler, SLOT(mesh_volume()));
  menu_list.push_back(a_vol);

  // Separator before analysis launcher
  QAction* sep2 = new QAction(handler);
  sep2->setSeparator(true);
  menu_list.push_back(sep2);

  // Radia-NGSolve analysis launcher (was in Python Solve menu, now C++)
  QAction* a_ngs = new QAction("Radia-NGSolve...", handler);
  a_ngs->setStatusTip("Export .vol and launch analysis window (IH / EM / PCB)");
  QObject::connect(a_ngs, SIGNAL(triggered()), handler, SLOT(launch_radia_ngsolve()));
  menu_list.push_back(a_ngs);

  gui->add_to_menu("&Export Mesh", menu_list, "radiacomp");

  // Python Solve menu is removed (2026-04-14). Radia-NGSolve is here in C++.

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
// Forward declarations for settings (defined later, after ExportDialog)
static QString settingsDir();
static QString settingsPath();
static QJsonObject loadSettings();

static QString default_start_dir()
{
  // Read default_dir from export_settings.json (set by Python startup)
  QJsonObject all = loadSettings();
  if (all.contains("default_dir")) {
    QString dd = all["default_dir"].toString();
    if (QDir(dd).exists())
      return dd;
  }
  return QString();
}

static QString ensure_model()
{
  if (CubitInterface::get_volume_count() > 0)
    return QString();

  ccl_log("ensure_model: no volumes, prompting for .jou");
  // No geometry — ask user to load a journal file
  QString jou = QFileDialog::getOpenFileName(
      nullptr, "No model loaded - Select Journal File",
      default_start_dir(), "Cubit Journal (*.jou);;All Files (*)");
  if (jou.isEmpty()) {
    ccl_log("ensure_model: cancelled");
    return QString();
  }

  jou.replace("\\", "/");
  ccl_log("ensure_model: playing %s", jou.toLocal8Bit().constData());
  std::string cmd = "play \"" + jou.toStdString() + "\"";
  CubitInterface::cmd(cmd.c_str());

  if (CubitInterface::get_volume_count() > 0) {
    ccl_log("ensure_model: OK (%d volumes)", CubitInterface::get_volume_count());
    return jou;
  }
  ccl_log("ensure_model: FAIL (0 volumes after play)");
  return QString();
}

// Persistent journal path across exports in one session
static QString s_lastJouPath;

//! Ensure a .jou path is known.  If no journal is associated with
//! the current model, prompt the user to save one.  Returns the
//! .jou path (empty if cancelled).
static bool _is_valid_jou_path(const QString& p)
{
  // Reject empty, just-a-separator, or paths missing a basename.
  if (p.isEmpty()) return false;
  if (p == "/" || p == "\\" || p == "." || p == "..") return false;
  // Must have a non-empty basename after the last slash and (optional)
  // .jou suffix, otherwise volPath derivation gives "/.vol".
  QFileInfo fi(p);
  if (fi.completeBaseName().isEmpty()) return false;
  return true;
}

static QString ensure_jou_path()
{
  // 1. Cubit's own journal path
  std::string jf = CubitInterface::get_current_journal_file();
  if (!jf.empty()) {
    QString p = QString::fromStdString(jf).replace("\\", "/");
    if (_is_valid_jou_path(p)) {
      ccl_log("ensure_jou_path: from Cubit journal: %s", jf.c_str());
      return p;
    }
    ccl_log("ensure_jou_path: Cubit returned invalid journal '%s' "
            "— forcing save prompt.", jf.c_str());
  }

  // 2. Previously loaded/saved path in this session
  if (!s_lastJouPath.isEmpty() && _is_valid_jou_path(s_lastJouPath)) {
    ccl_log("ensure_jou_path: from session: %s",
            s_lastJouPath.toLocal8Bit().constData());
    return s_lastJouPath;
  }

  // 3. No journal — prompt user to save one
  ccl_log("ensure_jou_path: no journal known, prompting save");
  QString defaultDir = default_start_dir();
  if (defaultDir.isEmpty()) {
    char cwd[1024];
    if (_getcwd(cwd, sizeof(cwd)))
      defaultDir = QString::fromLocal8Bit(cwd).replace("\\", "/");
  }
  QString jou = QFileDialog::getSaveFileName(
      nullptr, "Save Journal (determines output filenames)",
      defaultDir + "/model.jou",
      "Cubit Journal (*.jou);;All Files (*)");
  if (jou.isEmpty()) {
    ccl_log("ensure_jou_path: cancelled");
    return QString();
  }

  jou.replace("\\", "/");

  // Ensure parent directory exists
  QDir().mkpath(QFileInfo(jou).absolutePath());

  // Save journal only (.cub5 is NOT saved — user manages session manually)
  std::string cmd = "save journal \"" + jou.toStdString() + "\" overwrite";
  CubitInterface::silent_cmd(cmd.c_str());
  s_lastJouPath = jou;
  ccl_log("ensure_jou_path: saved %s", jou.toLocal8Bit().constData());
  PRINT_INFO("Saved journal: %s\n", jou.toLocal8Bit().constData());
  return jou;
}

static void run_export(ExportDialog::Format fmt)
{
  if (CubitInterface::get_volume_count() == 0) {
    QString jou = ensure_model();
    if (CubitInterface::get_volume_count() == 0)
      return;
    if (!jou.isEmpty())
      s_lastJouPath = jou;
  }

  // Ensure .jou is saved — determines output basename
  QString jouPath = ensure_jou_path();
  if (jouPath.isEmpty())
    return;

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

  ccl_log("run_export: cmd=[%s]", cmd.c_str());
  CubitInterface::cmd(cmd.c_str());

  QString outFile = dlg.filePath();
  bool exists = (fmt == ExportDialog::FEMEEM)
                  ? QDir(outFile).exists()
                  : QFile::exists(outFile);
  if (!exists) {
    ccl_log("run_export: FAIL %s", outFile.toStdString().c_str());
    PRINT_ERROR("Export failed: %s\n", outFile.toStdString().c_str());
    return;
  }
  ccl_log("run_export: OK %s", outFile.toStdString().c_str());
  PRINT_INFO("Export complete: %s\n", outFile.toStdString().c_str());

  // Open exported file/directory with OS-associated application
  QDesktopServices::openUrl(QUrl::fromLocalFile(outFile));
}

void RadiaMenuHandler::export_gmsh()    { run_export(ExportDialog::GMSH); }
void RadiaMenuHandler::export_nastran() { run_export(ExportDialog::Nastran); }
void RadiaMenuHandler::export_vtk()     { run_export(ExportDialog::VTK); }
void RadiaMenuHandler::export_femeem()  { run_export(ExportDialog::FEMEEM); }
void RadiaMenuHandler::export_meg()     { run_export(ExportDialog::MEG); }
// ============================================================
// Helpers for subprocess-based operations (Netgen export, Mesh Eval)
// ============================================================
static QString find_external_python();
static QString find_calc_script(const QString &python, const QString &name);

void RadiaMenuHandler::export_netgen()
{
  if (CubitInterface::get_volume_count() == 0) {
    QString jou = ensure_model();
    if (CubitInterface::get_volume_count() == 0) return;
    if (!jou.isEmpty())
      s_lastJouPath = jou;
  }

  // Ensure .jou is saved — determines output basename
  QString jouPath = ensure_jou_path();
  if (jouPath.isEmpty()) return;

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
  ccl_log("export_netgen: jou=%s order=%d",
          jouPath.toLocal8Bit().constData(), order);
  PRINT_INFO("Exporting Netgen .vol (order %d)...\n", order);
  std::string cmd = std::string("radia_export netgen \"")
    + volPath.toLocal8Bit().constData()
    + "\" order " + std::to_string(order) + " overwrite";
  CubitInterface::silent_cmd(cmd.c_str());

  // Check .vol was created
  if (!QFile::exists(volPath)) {
    ccl_log("export_netgen: FAIL %s", volPath.toLocal8Bit().constData());
    PRINT_ERROR("radia_export netgen command failed.\n");
    return;
  }
  ccl_log("export_netgen: OK %s", volPath.toLocal8Bit().constData());
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
    QString jou = ensure_model();
    if (CubitInterface::get_volume_count() == 0) return;
    if (!jou.isEmpty())
      s_lastJouPath = jou;
  }

  // Ensure .jou is saved — determines output basename
  QString jouPath = ensure_jou_path();
  if (jouPath.isEmpty()) return;

  ensure_netgen_dll_path();

  // Derive basename: /path/to/model.jou -> /path/to/model
  QString base = jouPath;
  if (base.endsWith(".jou", Qt::CaseInsensitive))
    base.chop(4);
  QString outDir = QFileInfo(base).absolutePath();
  QDir().mkpath(outDir);
  CubitInterface::silent_cmd(
      ("cd \"" + outDir.toLocal8Bit() + "\"").constData());

  // --- Phase 1: C++ exports (all Cubit operations here) ---
  // p-convergence: .vol at order 1-5 (each produces companion .vol.json)
  int maxOrder = 5;
  ccl_log("mesh_eval: jou=%s base=%s maxOrder=%d",
          jouPath.toLocal8Bit().constData(),
          base.toLocal8Bit().constData(), maxOrder);
  PRINT_INFO("Exporting .vol order 1-%d...\n", maxOrder);
  for (int p = 1; p <= maxOrder; p++) {
    QString volPath = QString("%1_p%2.vol").arg(base).arg(p);
    std::string cmd = "radia_export netgen \"" +
        volPath.toLocal8Bit().toStdString() +
        "\" order " + std::to_string(p) + " overwrite";
    CubitInterface::silent_cmd(cmd.c_str());
    if (!QFile::exists(volPath))
      PRINT_WARNING("order %d export failed\n", p);
    else
      PRINT_INFO("  p=%d: %s\n", p, volPath.toLocal8Bit().constData());
  }

  // Format QA: .msh (order 1-3), .bdf/.vtk (order 1-2)
  // GMSH order 4-5 is an error (NetgenCurver face/volume interior node
  // extraction is unreliable at p>=4).  Use radia_export netgen for p>=4.
  struct FmtSpec { const char *name; const char *ext; const char *cmd; int maxOrd; };
  int gmshMaxOrd = maxOrder < 3 ? maxOrder : 3;
  int bdfMaxOrd  = maxOrder < 2 ? maxOrder : 2;
  FmtSpec fmts[] = {
    {"gmsh",    "msh", "radia_export gmsh",    gmshMaxOrd},
    {"nastran", "bdf", "radia_export nastran",  bdfMaxOrd},
    {"vtk",     "vtk", "radia_export vtk",      bdfMaxOrd},
  };
  for (auto &fmt : fmts) {
    for (int p = 1; p <= fmt.maxOrd; p++) {
      QString fpath = QString("%1_qa_%2_o%3.%4")
          .arg(base).arg(fmt.name).arg(p).arg(fmt.ext);
      std::string cmd = std::string(fmt.cmd) + " \"" +
          fpath.toLocal8Bit().constData() + "\"";
      if (p >= 2)
        cmd += " order " + std::to_string(p);
      cmd += " overwrite";
      CubitInterface::silent_cmd(cmd.c_str());
    }
  }

  // --- Phase 2: subprocess — calc_mesh_eval.py reads .vol + .vol.json ---
  // NO import cubit in Python — all Cubit operations done above
  QString python = find_external_python();
  QString script = find_calc_script(python, "calc_mesh_eval");
  if (script.isEmpty()) {
    PRINT_ERROR("Cannot find calc_mesh_eval.py. Is radia installed?\n");
    return;
  }

  PRINT_INFO("Evaluating mesh with NGSolve + GMSH...\n");

  QProcess proc;
  QStringList args;
  if (python == "py") args << "-3";
  args << script << "--vol-base" << base
       << "--max-order" << QString::number(maxOrder);

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

  // --- Format QA table (GMSH API Jacobian vs CAD) ---
  QJsonArray fmtArr = r["format_qa"].toArray();
  int nfmt = fmtArr.size();
  if (nfmt > 0) {
    layout->addWidget(new QLabel(
        "\nFormat QA vs CAD (GMSH API Jacobian verification):"));
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
// Radia-NGSolve launcher — replaces Python register_toolbar.py
// ============================================================

// Metadata parsed from radia_*.py scripts
struct ModeScript {
  QString title;       // TITLE = "Induction Heating"
  QString fileName;    // radia_ih.py
  QStringList required; // REQUIRED_LABELS
  QStringList optional; // OPTIONAL_LABELS
  bool needsVol = true; // NEEDS_VOL — default true.  When false, the
                        // launcher hides the `.vol` row and skips the
                        // `radia_export netgen` call (e.g. PEEC-inductance
                        // which takes a STEP file, not a mesh).
};

// Find radia package directory via external Python
static QString find_radia_package_dir(const QString &python)
{
  QProcess p;
  QStringList args;
  if (python == "py") args << "-3";
  args << "-c"
       << "import importlib.util, os; "
          "spec = importlib.util.find_spec('radia'); "
          "print(os.path.dirname(spec.origin) if spec else '')";
  p.start(python, args);
  if (p.waitForFinished(10000) && p.exitCode() == 0) {
    QString dir = QString::fromUtf8(p.readAllStandardOutput()).trimmed();
    if (!dir.isEmpty() && QDir(dir).exists())
      return dir;
  }
  return QString();
}

// Discover radia_*.py analysis windows and their metadata
static QList<ModeScript> discover_mode_scripts(const QString &pkgDir)
{
  QList<ModeScript> modes;
  if (pkgDir.isEmpty()) return modes;

  QDir dir(pkgDir);
  QStringList files = dir.entryList({"radia_*.py"}, QDir::Files, QDir::Name);

  // Match Python string assignment: TITLE = "..." or TITLE = '...'
  // Use backreference to ensure matching quote type (handles apostrophes).
  QRegularExpression reTitle(R"(^TITLE\s*=\s*(["'])(.+?)\1)");
  QRegularExpression reReq(R"(^REQUIRED_LABELS\s*=\s*\[([^\]]*)\])");
  QRegularExpression reOpt(R"(^OPTIONAL_LABELS\s*=\s*\[([^\]]*)\])");
  QRegularExpression reNeedsVol(R"(^NEEDS_VOL\s*=\s*(True|False))");

  auto parseList = [](const QString &raw) -> QStringList {
    QStringList result;
    // Parse Python list literal: ["a", "b", 'c']
    // Backreference \1 ensures matching quote type.
    QRegularExpression reItem(R"((["'])([^"']+?)\1)");
    auto it = reItem.globalMatch(raw);
    while (it.hasNext())
      result << it.next().captured(2);
    return result;
  };

  for (const QString &fname : files) {
    if (fname == "radia_gui_base.py") continue;

    QFile f(dir.filePath(fname));
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) continue;

    ModeScript ms;
    ms.fileName = fname;

    while (!f.atEnd()) {
      QString line = QString::fromUtf8(f.readLine());
      auto m = reTitle.match(line);
      if (m.hasMatch()) ms.title = m.captured(2);
      m = reReq.match(line);
      if (m.hasMatch()) ms.required = parseList(m.captured(1));
      m = reOpt.match(line);
      if (m.hasMatch()) ms.optional = parseList(m.captured(1));
      m = reNeedsVol.match(line);
      if (m.hasMatch()) ms.needsVol = (m.captured(1) == "True");
    }
    if (!ms.title.isEmpty())
      modes.append(ms);
  }
  return modes;
}

// Thin Qt wrapper around RadiaLauncherLogic::get_model_labels.
// All actual API calls live in LauncherLogic.hpp so that the headless
// `radia_export verify_launcher` command tests the SAME code that the
// dialog runs (2026-04-14 hardening — see feedback_gui_test_with_deploy).
static QStringList get_model_labels()
{
  QStringList out;
  for (auto &n : RadiaLauncherLogic::get_model_labels())
    out.append(QString::fromStdString(n));
  return out;
}

// Launcher settings persistence (separate from export settings)
static QString launcherSettingsPath()
{
  QString home = QDir::homePath() + "/.radia";
  QDir().mkpath(home);
  return home + "/radia_launcher.json";
}

static QJsonObject loadLauncherSettings()
{
  QFile f(launcherSettingsPath());
  if (!f.open(QIODevice::ReadOnly)) return {};
  return QJsonDocument::fromJson(f.readAll()).object();
}

static void saveLauncherSettings(const QJsonObject &obj)
{
  QJsonObject existing = loadLauncherSettings();
  for (auto it = obj.begin(); it != obj.end(); ++it)
    existing[it.key()] = it.value();
  QFile f(launcherSettingsPath());
  if (f.open(QIODevice::WriteOnly))
    f.write(QJsonDocument(existing).toJson(QJsonDocument::Indented));
}

void RadiaMenuHandler::launch_radia_ngsolve()
{
  // --- Ensure model loaded ---
  if (CubitInterface::get_volume_count() == 0) {
    QString jou = ensure_model();
    if (CubitInterface::get_volume_count() == 0) return;
    if (!jou.isEmpty())
      s_lastJouPath = jou;
  }

  ensure_netgen_dll_path();

  // --- Find external Python and radia package ---
  QString python = find_external_python();
  QString pkgDir = find_radia_package_dir(python);
  if (pkgDir.isEmpty()) {
    PRINT_ERROR("Cannot find radia package. Is radia installed? "
                "(pip install radia)\n");
    return;
  }

  // --- Discover analysis windows ---
  QList<ModeScript> modes = discover_mode_scripts(pkgDir);
  if (modes.isEmpty()) {
    PRINT_ERROR("No radia_*.py analysis windows found in %s\n",
                pkgDir.toStdString().c_str());
    return;
  }

  // Ensure .jou is saved — determines output basename
  QString jouPath = ensure_jou_path();
  if (jouPath.isEmpty()) return;

  // --- Load saved settings ---
  QJsonObject prev = loadLauncherSettings();
  QString lastMode = prev["last_mode"].toString();
  int lastOrder = prev["last_order"].toInt(2);

  // --- Build dialog ---
  QDialog dlg;
  dlg.setWindowTitle("Radia-NGSolve");
  dlg.setMinimumWidth(500);
  QVBoxLayout *layout = new QVBoxLayout(&dlg);

  // Mode combo
  QHBoxLayout *hMode = new QHBoxLayout();
  hMode->addWidget(new QLabel("Analysis:"));
  QComboBox *modeCombo = new QComboBox();
  for (const auto &ms : modes)
    modeCombo->addItem(ms.title);
  if (!lastMode.isEmpty()) {
    int idx = modeCombo->findText(lastMode);
    if (idx >= 0) modeCombo->setCurrentIndex(idx);
  }
  hMode->addWidget(modeCombo);
  hMode->addStretch();
  layout->addLayout(hMode);

  // Order combo (wrapped in a container so the whole row can be
  // hidden for STEP-only modes that do not export a mesh).
  QWidget *orderRow = new QWidget();
  QHBoxLayout *hOrder = new QHBoxLayout(orderRow);
  hOrder->setContentsMargins(0, 0, 0, 0);
  hOrder->addWidget(new QLabel("Mesh order:"));
  QComboBox *orderCombo = new QComboBox();
  for (int p = 1; p <= 5; p++)
    orderCombo->addItem(QString::number(p));
  orderCombo->setCurrentIndex(qBound(0, lastOrder - 1, 4));
  hOrder->addWidget(orderCombo);
  hOrder->addStretch();
  layout->addWidget(orderRow);

  // Auto-Kelvin checkbox.  When checked (default) and the current
  // Cubit model has an "air" block but no "kelvin" block, the launcher
  // injects `add_kelvin_cubit()` via the auto_add_kelvin_from_current_
  // model() helper just before `radia_export netgen`.  Uncheck for
  // no-Kelvin baselines (Dirichlet truncation on the outer air
  // boundary) or for BEM modes where Kelvin is unnecessary.
  // Settings persist via launcher.json alongside last_mode / last_order.
  QWidget *kelvinRow = new QWidget();
  QHBoxLayout *hKelvin = new QHBoxLayout(kelvinRow);
  hKelvin->setContentsMargins(0, 0, 0, 0);
  QCheckBox *kelvinCheck = new QCheckBox(
      "Add Kelvin open boundary (auto)");
  bool lastKelvinAuto = prev.contains("last_kelvin_auto")
      ? prev["last_kelvin_auto"].toBool() : true;
  kelvinCheck->setChecked(lastKelvinAuto);
  kelvinCheck->setToolTip(
      "Before exporting the .vol, auto-add a Kelvin exterior sphere "
      "to the Cubit model (needs an 'air' block).  Idempotent — skips "
      "if a 'kelvin' block already exists.  Uncheck for Dirichlet-"
      "truncation baselines.");
  hKelvin->addWidget(kelvinCheck);

  // Kelvin exterior mesh size (optional).  Blank means "auto" — let
  // add_kelvin_cubit inherit the size from the air outer surface via
  // copy-mesh (the current default).  A typical override is a value
  // 2-3x the air surface size, since the Kelvin region is just an
  // open-boundary trick and can be much coarser than the physical
  // domain.  Value in meters.  Persisted across sessions via
  // launcher.json / last_kelvin_mesh.
  hKelvin->addSpacing(20);
  hKelvin->addWidget(new QLabel("Mesh size [m]:"));
  QLineEdit *kelvinMeshEdit = new QLineEdit();
  kelvinMeshEdit->setPlaceholderText("auto (blank) or e.g. 0.03");
  kelvinMeshEdit->setMaximumWidth(130);
  kelvinMeshEdit->setToolTip(
      "Tet edge length for the Kelvin exterior sphere.  Blank = auto "
      "(inherit from air surface).  For faster meshing / smaller "
      ".vol, specify a value coarser than your air mesh.");
  QString lastKelvinMesh = prev.contains("last_kelvin_mesh")
      ? prev["last_kelvin_mesh"].toString() : QString();
  kelvinMeshEdit->setText(lastKelvinMesh);
  hKelvin->addWidget(kelvinMeshEdit);
  hKelvin->addStretch();
  layout->addWidget(kelvinRow);

  // Label validation area
  QLabel *labelWidget = new QLabel();
  labelWidget->setWordWrap(true);
  labelWidget->setTextFormat(Qt::RichText);
  QGroupBox *labelGroup = new QGroupBox("Labels (blocks / sidesets)");
  QVBoxLayout *labelLayout = new QVBoxLayout(labelGroup);
  labelLayout->addWidget(labelWidget);
  layout->addWidget(labelGroup);

  // Output .vol path: derived from the journal path (.jou -> .vol in the
  // same directory). Editable + Browse so the user can override.
  // 2026-04-14: ensure_jou_path now refuses to return "/" or other
  // unusable paths, so this derivation cannot produce "/.vol" anymore.
  QString volBase = jouPath;
  if (volBase.endsWith(".jou", Qt::CaseInsensitive))
    volBase.chop(4);
  QString defaultVol = volBase + ".vol";
  QWidget *volRow = new QWidget();
  QHBoxLayout *hVol = new QHBoxLayout(volRow);
  hVol->setContentsMargins(0, 0, 0, 0);
  hVol->addWidget(new QLabel("Output .vol:"));
  QLineEdit *volEdit = new QLineEdit(defaultVol);
  hVol->addWidget(volEdit, 1);
  QPushButton *volBrowseBtn = new QPushButton("Browse...");
  hVol->addWidget(volBrowseBtn);
  QObject::connect(volBrowseBtn, &QPushButton::clicked, [&dlg, volEdit]() {
    QString picked = QFileDialog::getSaveFileName(
        &dlg, "Output .vol path", volEdit->text(),
        "Netgen Vol (*.vol);;All Files (*)");
    if (!picked.isEmpty()) {
      picked.replace("\\", "/");
      volEdit->setText(picked);
    }
  });
  layout->addWidget(volRow);

  // OK / Cancel
  QDialogButtonBox *buttons = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dlg);
  QPushButton *okBtn = buttons->button(QDialogButtonBox::Ok);
  QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &dlg, &QDialog::reject);
  layout->addWidget(buttons);

  // Label update lambda
  auto updateLabels = [&]() {
    int idx = modeCombo->currentIndex();
    if (idx < 0 || idx >= modes.size()) return;
    const ModeScript &ms = modes[idx];

    // Hide the `.vol` row, Mesh order, and label-requirements entirely
    // when this mode does not need a mesh (e.g. PEEC-inductance takes a
    // STEP, not a mesh).  Showing widgets the user cannot act on ("Mesh
    // order: 2" while no mesh is exported) is confusing UX.
    volRow->setVisible(ms.needsVol);
    orderRow->setVisible(ms.needsVol);
    labelGroup->setVisible(ms.needsVol);

    QStringList modelLabels = get_model_labels();
    QSet<QString> labelSet;
    for (const QString &l : modelLabels)
      labelSet.insert(l.toLower());

    QStringList lines;
    bool allOk = true;
    for (const QString &lbl : ms.required) {
      if (labelSet.contains(lbl.toLower())) {
        lines << QString("<span style=\"color:green\">"
                         "[OK] %1 (required)</span>").arg(lbl);
      } else {
        lines << QString("<span style=\"color:red; font-weight:bold\">"
                         "[MISSING] %1 (required)</span>").arg(lbl);
        allOk = false;
      }
    }
    for (const QString &lbl : ms.optional) {
      if (labelSet.contains(lbl.toLower())) {
        lines << QString("<span style=\"color:green\">"
                         "[OK] %1 (optional)</span>").arg(lbl);
      } else {
        lines << QString("<span style=\"color:gray\">"
                         "[ - ] %1 (optional)</span>").arg(lbl);
      }
    }
    if (ms.required.isEmpty() && ms.optional.isEmpty())
      lines << "<span style=\"color:gray\">No label requirements</span>";

    labelWidget->setText(lines.join("<br>"));
    okBtn->setEnabled(allOk);
  };

  QObject::connect(modeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                   [&](int) { updateLabels(); });
  updateLabels();

  // --- Execute dialog ---
  if (dlg.exec() != QDialog::Accepted)
    return;

  int modeIdx = modeCombo->currentIndex();
  const ModeScript &ms = modes[modeIdx];
  int order = orderCombo->currentIndex() + 1;

  QString volPath;
  QString outDir;
  if (ms.needsVol) {
    // Use the (possibly user-edited) Output .vol path.
    volPath = volEdit->text().trimmed().replace("\\", "/");
    if (volPath.isEmpty() || volPath == "/" || volPath == "/.vol") {
      PRINT_ERROR("Output .vol path is empty or invalid: %s\n",
                  volPath.toLocal8Bit().constData());
      return;
    }
    outDir = QFileInfo(volPath).absolutePath();
  } else {
    // STEP-only modes (e.g. PEEC-inductance): no .vol, the analysis
    // window runs in the journal's directory.
    outDir = QFileInfo(jouPath).absolutePath();
  }
  if (outDir.isEmpty()) outDir = ".";

  bool kelvinAuto = kelvinCheck->isChecked();
  QString kelvinMeshStr = kelvinMeshEdit->text().trimmed();
  // Parse kelvin mesh size: empty -> auto (null in JSON), else float.
  bool hasMeshSize = false;
  double kelvinMeshSize = 0.0;
  if (!kelvinMeshStr.isEmpty()) {
    bool ok = false;
    double v = kelvinMeshStr.toDouble(&ok);
    if (ok && v > 0.0) { hasMeshSize = true; kelvinMeshSize = v; }
    else {
      PRINT_WARNING("Kelvin mesh size '%s' is not a positive number; "
                    "using auto (blank).\n",
                    kelvinMeshStr.toLocal8Bit().constData());
    }
  }

  // Save settings (persist last values for next dialog open)
  saveLauncherSettings({
      {"last_mode", ms.title},
      {"last_order", order},
      {"last_jou_dir", outDir},
      {"last_kelvin_auto", kelvinAuto},
      {"last_kelvin_mesh", kelvinMeshStr},  // store raw text so blank = blank
  });

  // --- Export .vol (only when the mode needs it) ---
  if (ms.needsVol) {
    if (!QDir().mkpath(outDir)) {
      PRINT_ERROR("Cannot create output directory: %s\n",
                  outDir.toLocal8Bit().constData());
      return;
    }
    CubitInterface::silent_cmd(
        ("cd \"" + outDir.toLocal8Bit() + "\"").constData());

    // --- Auto-Kelvin: add open-boundary sphere if requested ---
    // Plays radia/panels/auto_kelvin_entry.py, a tiny wrapper that
    // imports add_kelvin.auto_add_kelvin_from_current_model() and
    // invokes it.  The helper is idempotent (skips if a "kelvin"
    // block already exists) and degrades gracefully (warns + returns
    // None) when no "air" block is present or detection fails.
    // Runs BEFORE `radia_export netgen` so the Kelvin block ends up
    // in the .vol output.
    //
    // We use `play "<path>.py"` rather than `python "stmt"` because
    // the `python` keyword is not recognised at the top level of
    // Cubit's journal parser (verified 2026-04-23), while `play`
    // works for both .jou and .py files.
    if (kelvinAuto) {
      QString panelsDirFs = pkgDir + "/panels";
      panelsDirFs.replace("\\", "/");
      QString entryScript = panelsDirFs + "/auto_kelvin_entry.py";
      if (!QFile::exists(entryScript)) {
        PRINT_WARNING("Auto-Kelvin skipped: %s not found\n",
                      entryScript.toLocal8Bit().constData());
      } else {
        // Write the entry-script config JSON so the Python helper
        // picks up dialog state (add_kelvin + kelvin_mesh_size).
        // Schema matches tests/panels/test_auto_kelvin_cli.py exactly.
        QString cfgPath = outDir + "/radia_launcher_config.json";
        cfgPath.replace("\\", "/");
        {
          QJsonObject cfg;
          cfg["add_kelvin"]        = true;   // kelvinAuto == true here
          cfg["kelvin_air_block"]  = QString("air");
          cfg["kelvin_block_name"] = QString("kelvin");
          if (hasMeshSize)
            cfg["kelvin_mesh_size"] = kelvinMeshSize;
          else
            cfg["kelvin_mesh_size"] = QJsonValue();  // null = auto
          QFile f(cfgPath);
          if (f.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
            f.write(QJsonDocument(cfg).toJson(QJsonDocument::Compact));
            f.close();
            PRINT_INFO("Auto-Kelvin: config -> %s (mesh_size=%s)\n",
                       cfgPath.toLocal8Bit().constData(),
                       hasMeshSize
                         ? QString::number(kelvinMeshSize).toLocal8Bit()
                                   .constData()
                         : "auto");
          } else {
            PRINT_WARNING("Auto-Kelvin: cannot write config %s, "
                          "entry will run with defaults\n",
                          cfgPath.toLocal8Bit().constData());
          }
        }

        // Expose the config path to Cubit's embedded Python via the
        // process environment block.  auto_kelvin_entry.py reads
        // $RADIA_LAUNCHER_CONFIG.  SetEnvironmentVariable mutates the
        // current process, which is what os.environ reads from.
        SetEnvironmentVariableW(L"RADIA_LAUNCHER_CONFIG",
                                (LPCWSTR)cfgPath.utf16());
        // Cubit's `play` exec's scripts without binding __file__,
        // so auto_kelvin_entry.py cannot locate its sibling
        // add_kelvin.py via os.path.abspath(__file__).  Pass the
        // panels directory explicitly via RADIA_PANELS_DIR.
        SetEnvironmentVariableW(L"RADIA_PANELS_DIR",
                                (LPCWSTR)panelsDirFs.utf16());

        PRINT_INFO("Auto-Kelvin: playing %s\n",
                   entryScript.toLocal8Bit().constData());
        std::string playCmd = "play \"" +
            entryScript.toLocal8Bit().toStdString() + "\"";
        CubitInterface::cmd(playCmd.c_str());

        // Clear both env vars so they don't leak into subsequent
        // invocations (user may re-open the dialog with different
        // settings in the same Cubit session).
        SetEnvironmentVariableW(L"RADIA_LAUNCHER_CONFIG", nullptr);
        SetEnvironmentVariableW(L"RADIA_PANELS_DIR", nullptr);
      }
    }

    PRINT_INFO("Exporting Netgen .vol (order %d) to %s...\n",
               order, volPath.toLocal8Bit().constData());

    std::string exportCmd = "radia_export netgen \"" +
        volPath.toLocal8Bit().toStdString() +
        "\" order " + std::to_string(order) + " overwrite";
    CubitInterface::cmd(exportCmd.c_str());

    if (!QFile::exists(volPath)) {
      PRINT_ERROR("radia_export netgen failed. Check blocks/sidesets.\n");
      return;
    }
    PRINT_INFO("Exported: %s (order %d)\n",
               volPath.toLocal8Bit().constData(), order);
  } else {
    PRINT_INFO("Launching %s (STEP-only mode, no .vol export)...\n",
               ms.title.toLocal8Bit().constData());
  }

  // --- Launch analysis window subprocess ---
  QString scriptPath = pkgDir + "/" + ms.fileName;
  if (!QFile::exists(scriptPath)) {
    PRINT_ERROR("Analysis script not found: %s\n",
                scriptPath.toStdString().c_str());
    return;
  }

  QStringList args;
  if (python == "py") args << "-3";
  args << scriptPath;
  // Only pass the .vol path when the mode uses it; STEP-only modes
  // launch with no positional argument (empty vol_path in Python).
  if (ms.needsVol) args << volPath;

  PRINT_INFO("Launching: %s %s\n",
             python.toStdString().c_str(),
             args.join(" ").toStdString().c_str());

  // Static startDetached: no heap QProcess needed, no memory leak.
  // Environment cleaning (Qt/PySide) is handled by the PySide6 app
  // itself (it imports its own Qt, ignoring Cubit's env vars).
  qint64 pid = 0;
  QProcess::startDetached(python, args, outDir, &pid);
  if (pid > 0)
    PRINT_INFO("PID: %lld\n", pid);
  else
    PRINT_WARNING("Failed to launch analysis window.\n");

  // Open the exported .vol with OS viewer (only for mesh modes).
  if (ms.needsVol)
    QDesktopServices::openUrl(QUrl::fromLocalFile(volPath));
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
    mVersion(nullptr), mNoPyramid(nullptr), mScale(nullptr)
{
  // Window title
  const char* titles[] = {"Export Netgen Vol", "Export GMSH", "Export Nastran BDF", "Export VTK", "Export FEMEEM", "Export MEG (ELF/MAGIC)"};
  setWindowTitle(titles[format]);
  setMinimumWidth(500);

  QVBoxLayout* layout = new QVBoxLayout(this);
  QFormLayout* form = new QFormLayout();

  // Determine default directory and base name
  // FEMEEM: no extension (directory name), others: file extension
  const char* exts[] = {".vol", ".msh", ".bdf", ".vtk", "", ".meg"};
  QString defaultDir;
  QString baseName = (format == FEMEEM) ? "femeem_output" : "ExportedMesh";

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

  if (format == FEMEEM) {
    // FEMEEM: single directory path (files in.dat etc. are fixed names)
    QHBoxLayout* dirRow = new QHBoxLayout();
    mDir = new QLineEdit(defaultDir + "/" + baseName);
    QPushButton* browseBtn = new QPushButton("...");
    browseBtn->setFixedWidth(30);
    connect(browseBtn, SIGNAL(clicked()), this, SLOT(browseDir()));
    connect(mDir, SIGNAL(textChanged(QString)), this, SLOT(updatePreview()));
    dirRow->addWidget(mDir);
    dirRow->addWidget(browseBtn);
    form->addRow("Directory:", dirRow);

    // Hidden mFileName (filePath() uses dir + filename, keep empty)
    mFileName = new QLineEdit("");
    mFileName->setVisible(false);
  } else {
    // Other formats: directory + filename
    QHBoxLayout* dirRow = new QHBoxLayout();
    mDir = new QLineEdit(defaultDir);
    QPushButton* browseBtn = new QPushButton("...");
    browseBtn->setFixedWidth(30);
    connect(browseBtn, SIGNAL(clicked()), this, SLOT(browseDir()));
    connect(mDir, SIGNAL(textChanged(QString)), this, SLOT(updatePreview()));
    dirRow->addWidget(mDir);
    dirRow->addWidget(browseBtn);
    form->addRow("Directory:", dirRow);

    mFileName = new QLineEdit(baseName + exts[format]);
    connect(mFileName, SIGNAL(textChanged(QString)), this, SLOT(updatePreview()));
    form->addRow("Filename:", mFileName);
  }

  // Order (not for FEMEEM/MEG — always 1st order)
  mOrderCombo = new QComboBox();
  if (format == FEMEEM || format == MEG) {
    mOrderCombo->addItems({"1 (linear)"});
    // No form row — hidden, always order 1
  } else if (format == NETGEN_VOL) {
    mOrderCombo->addItems({"1", "2", "3", "4", "5"});
    mOrderCombo->setCurrentIndex(2);  // default 3
    connect(mOrderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
    form->addRow("Order:", mOrderCombo);
  } else if (format == GMSH) {
    // GMSH: tet/hex up to order 3, wedge up to order 2.
    // Order 4-5 unsupported (NetgenCurver interior node extraction
    // unreliable at p>=4). Use Netgen Vol for order 4-5.
    mOrderCombo->addItems({"1", "2", "3 (tet/hex only, wedge stays 2)"});
    mOrderCombo->setCurrentIndex(0);  // default 1
    connect(mOrderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
    form->addRow("Order:", mOrderCombo);
  } else {
    // Nastran, VTK: order 1-2 (format specification limit)
    mOrderCombo->addItems({"1 (linear)", "2 (quadratic)"});
    connect(mOrderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
    form->addRow("Order:", mOrderCombo);
  }

  // GMSH output: always v4.1 (lab-wide standard, 2026-04).

  // Scale (FEMEEM only)
  if (format == FEMEEM) {
    mScale = new QLineEdit("1.0");
    connect(mScale, SIGNAL(textChanged(QString)), this, SLOT(updatePreview()));
    form->addRow("Scale:", mScale);
  }

  // Dimension
  mDimension = new QComboBox();
  if (format == MEG) {
    mDimension->addItems({"3D (threed)", "2D (twod)", "Axisymmetric"});
  } else {
    mDimension->addItems({"3D", "2D"});
  }
  connect(mDimension, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
  if (format != NETGEN_VOL && format != FEMEEM) {
    form->addRow("Dimension:", mDimension);
  }

  // NoPyramid (Nastran only)
  if (format == Nastran) {
    mNoPyramid = new QComboBox();
    mNoPyramid->addItems({"Keep pyramids", "Convert to degenerate hex (JMAG)"});
    connect(mNoPyramid, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
    form->addRow("Pyramids:", mNoPyramid);
  }

  // Block label table (MEG only)
  mBlockTable = nullptr;
  if (format == MEG) {
    layout->addLayout(form);
    QLabel *lbl = new QLabel("Block -> ELF Label:");
    layout->addWidget(lbl);
    mBlockTable = new QTableWidget(0, 3);
    mBlockTable->setHorizontalHeaderLabels({"Block", "Name", "ELF Type"});
    mBlockTable->horizontalHeader()->setStretchLastSection(true);
    mBlockTable->setColumnWidth(0, 50);
    mBlockTable->setColumnWidth(1, 140);
    mBlockTable->setMinimumHeight(120);
    mBlockTable->setMaximumHeight(250);
    layout->addWidget(mBlockTable);
    populateBlockTable();
    // Update labels when DIM changes
    connect(mDimension, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
  } else {
    layout->addLayout(form);
  }

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
    const char* keys[] = {"netgen_vol", "gmsh", "nastran", "vtk", "femeem", "meg"};
    QJsonObject all = loadSettings();
    QJsonObject s;
    s["dir"] = mDir->text();
    s["order"] = mOrderCombo->currentIndex();
    if (mDimension) s["dimension"] = mDimension->currentIndex();
    if (mNoPyramid) s["nopyramid"] = mNoPyramid->currentIndex();
    if (mScale) s["scale"] = mScale->text();
    all[keys[mFormat]] = s;
    saveSettings(all);
    accept();
  });
  layout->addWidget(buttons);

  // Load saved settings (override defaults with previous values)
  {
    const char* keys[] = {"netgen_vol", "gmsh", "nastran", "vtk", "femeem", "meg"};
    QJsonObject all = loadSettings();
    QJsonObject s = all[keys[format]].toObject();
    if (s.contains("dir")) mDir->setText(s["dir"].toString());
    if (s.contains("order")) mOrderCombo->setCurrentIndex(s["order"].toInt());
    if (s.contains("dimension") && mDimension)
      mDimension->setCurrentIndex(s["dimension"].toInt());
    if (s.contains("nopyramid") && mNoPyramid)
      mNoPyramid->setCurrentIndex(s["nopyramid"].toInt());
    if (s.contains("scale") && mScale)
      mScale->setText(s["scale"].toString());
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
  if (mFormat == FEMEEM) {
    // FEMEEM: mDir IS the full output directory path
    while (dir.endsWith('/')) dir.chop(1);
    return dir;
  }
  if (!dir.endsWith('/')) dir += '/';
  return dir + mFileName->text();
}

// --- ELF label helpers (must precede cubitCommand) ---

// ELF label options per DIM
static QStringList elfLabelsForDim(int dimIdx)
{
  // dimIdx: 0=3D(T), 1=2D(K), 2=Axisym(R)
  // "(none)" = skip this block in MEG export (e.g. air regions)
  if (dimIdx == 0) {
    return {"(none) - Skip export",
            "MMB - Magnetic body",
            "MMS - Magnetic shell",
            "MMT - Magnetic thin",
            "MMP - Permanent magnet",
            "MWL - Winding/coil",
            "MWV - Vector winding",
            "MCL - Coil line",
            "MCO - Current conductor",
            "MAB - Armature block",
            "MAT - Armature thin",
            "MBB - Boundary"};
  } else {
    return {"(none) - Skip export",
            "MMB - Magnetic body",
            "MMT - Magnetic thin",
            "MMP - Permanent magnet",
            "MWL - Winding/coil",
            "MWV - Vector winding",
            "MCL - Coil line",
            "MCO - Current conductor",
            "MAB - Armature block",
            "MAT - Armature thin",
            "MBB - Boundary"};
  }
}

// Extract 3-char prefix from combo text "MMB - Magnetic body" -> "MMB"
static QString elfPrefixFromCombo(const QString &text)
{
  return text.left(3);
}

// --- end ELF label helpers ---

QString ExportDialog::cubitCommand() const
{
  QString file = filePath();
  int order = mOrderCombo->currentIndex() + 1;

  QString cmd;
  switch (mFormat) {
    case NETGEN_VOL: {
      // Not a Cubit APREPRO command — shown as preview only
      cmd = QString("radia_export netgen \"%1\" order %2").arg(file).arg(order);
      break;
    }
    case GMSH: {
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("radia_export gmsh \"%1\" order %2 dimension %3")
                .arg(file).arg(order).arg(dim);
      break;
    }
    case Nastran: {
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("radia_export nastran \"%1\" order %2 dimension %3")
                .arg(file).arg(order).arg(dim);
      if (mNoPyramid && mNoPyramid->currentIndex() == 1)
        cmd += " nopyramid";
      break;
    }
    case VTK: {
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("radia_export vtk \"%1\" order %2 dimension %3")
                .arg(file).arg(order).arg(dim);
      break;
    }
    case FEMEEM: {
      QString sc = mScale ? mScale->text() : "1.0";
      cmd = QString("radia_export femeem \"%1\" scale %2").arg(file).arg(sc);
      break;
    }
    case MEG: {
      // MEG dimension: 3D -> threed, 2D -> twod, Axisymmetric -> axisymmetric
      QString dimOpt = "threed";
      if (mDimension) {
        int idx = mDimension->currentIndex();
        if (idx == 1) dimOpt = "twod";
        else if (idx == 2) dimOpt = "axisymmetric";
      }
      cmd = QString("radia_export meg \"%1\" %2").arg(file).arg(dimOpt);
      // Encode per-block labels: "1:MMB,2:MWL,3:MCO"
      // Blocks with "(none)" are skipped (not exported)
      if (mBlockTable && mBlockTable->rowCount() > 0) {
        QStringList pairs;
        for (int r = 0; r < mBlockTable->rowCount(); r++) {
          auto *idItem = mBlockTable->item(r, 0);
          auto *combo = qobject_cast<QComboBox*>(mBlockTable->cellWidget(r, 2));
          if (idItem && combo) {
            QString prefix = elfPrefixFromCombo(combo->currentText());
            if (prefix == "(no")  // "(none) - Skip export"
              continue;
            QString bid = idItem->text();
            pairs << bid + ":" + prefix;
          }
        }
        if (!pairs.isEmpty())
          cmd += QString(" labels \"%1\"").arg(pairs.join(","));
      }
      break;
    }
  }
  cmd += " overwrite";
  return cmd;
}

double ExportDialog::scale() const
{
  if (!mScale) return 1.0;
  bool ok;
  double v = mScale->text().toDouble(&ok);
  return (ok && v > 0.0) ? v : 1.0;
}

void ExportDialog::populateBlockTable()
{
  if (!mBlockTable) return;

  std::vector<int> block_ids = CubitInterface::parse_cubit_list("block", "all");
  int nrows = (int)block_ids.size();
  mBlockTable->setRowCount(nrows);

  int dimIdx = mDimension ? mDimension->currentIndex() : 0;
  QStringList labels = elfLabelsForDim(dimIdx);

  for (int r = 0; r < nrows; r++) {
    int bid = block_ids[r];
    // Use get_block_name (NOT get_entity_name("block", ...)) — the
    // generic accessor calls get_mref_entity which raises
    // "Invalid list type for the get_mref_entity function: block"
    // because blocks are groups, not mref entities.
    std::string name = CubitInterface::get_block_name(bid);
    if (name.empty()) name = "(unnamed)";

    // Block ID (read-only)
    auto *idItem = new QTableWidgetItem(QString::number(bid));
    idItem->setFlags(idItem->flags() & ~Qt::ItemIsEditable);
    mBlockTable->setItem(r, 0, idItem);

    // Block name (read-only)
    auto *nameItem = new QTableWidgetItem(QString::fromStdString(name));
    nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
    mBlockTable->setItem(r, 1, nameItem);

    // ELF type dropdown
    auto *combo = new QComboBox();
    combo->addItems(labels);

    // Auto-detect ELF label from block name.
    // "air", "kelvin" etc. -> (none) (skip export).
    // "MMB_iron", "MWL_coil" etc. -> match by 3-char prefix.
    // Unmatched -> MMB (first real label after "(none)").
    QString qname = QString::fromStdString(name).toUpper();
    int defaultIdx = 1;  // MMB (index 1, after "(none)" at index 0)
    // Check for names that should be skipped
    if (qname == "AIR" || qname == "KELVIN" || qname.startsWith("AIR_")
        || qname.startsWith("KELVIN_")) {
      defaultIdx = 0;  // (none)
    } else {
      for (int i = 1; i < labels.size(); i++) {
        if (qname.startsWith(elfPrefixFromCombo(labels[i]))) {
          defaultIdx = i;
          break;
        }
      }
    }
    combo->setCurrentIndex(defaultIdx);
    connect(combo, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
    mBlockTable->setCellWidget(r, 2, combo);
  }
}

void ExportDialog::updateBlockLabels()
{
  if (!mBlockTable || !mDimension) return;

  int dimIdx = mDimension->currentIndex();
  QStringList labels = elfLabelsForDim(dimIdx);

  for (int r = 0; r < mBlockTable->rowCount(); r++) {
    auto *combo = qobject_cast<QComboBox*>(mBlockTable->cellWidget(r, 2));
    if (!combo) continue;
    QString current = elfPrefixFromCombo(combo->currentText());
    combo->blockSignals(true);
    combo->clear();
    combo->addItems(labels);
    // Restore previous selection if still valid
    for (int i = 0; i < labels.size(); i++) {
      if (elfPrefixFromCombo(labels[i]) == current) {
        combo->setCurrentIndex(i);
        break;
      }
    }
    combo->blockSignals(false);
  }
}

void ExportDialog::updatePreview()
{
  // Update block labels when DIM changes (MEG only)
  if (mFormat == MEG && mBlockTable)
    updateBlockLabels();

  mPreview->setText(cubitCommand());
  mPreview->setCursorPosition(0);  // show beginning of command
}
