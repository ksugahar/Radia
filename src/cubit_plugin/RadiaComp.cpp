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
#include "ExportMegCommand.hpp"

#include <direct.h>  // _getcwd

#include <QAction>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QDir>
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
#include <QVBoxLayout>
#include <vector>

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
  // Menu order follows lab priority: .vol > .msh > .bdf > .vtk > .meg
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

  QAction* a_meg = new QAction("MEG...", handler);
  a_meg->setStatusTip("Export mesh to MEG/ELF format (.meg)");
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

  ExportDialog dlg(fmt, jouPath);
  if (dlg.exec() != QDialog::Accepted)
    return;

  std::string file = dlg.filePath().toStdString();
  int order = dlg.order();

  bool ok = false;
  switch (fmt) {
    case ExportDialog::GMSH: {
      ExportGmshCommand cmd;
      std::string ver = dlg.gmshVersion() >= 4 ? "4.1" : "2.2";
      ok = cmd.write_gmsh(file, ver, order);
      break;
    }
    case ExportDialog::Nastran: {
      ExportNastranCommand cmd;
      std::string dim = dlg.dimension() == 2 ? "2" : "3";
      ok = cmd.write_nastran(file, dim, !dlg.noPyramid(), order);
      break;
    }
    case ExportDialog::VTK: {
      ExportVtkCommand cmd;
      std::string dim = dlg.dimension() == 2 ? "2" : "3";
      ok = cmd.write_vtk(file, dim, order);
      break;
    }
    case ExportDialog::MEG: {
      ExportMegCommand cmd;
      // Extract T/K/R from combo text (first char)
      char megDim = dlg.megDimension();
      ok = cmd.write_meg(file, megDim);
      break;
    }
  }

  if (ok) {
    PRINT_INFO("Export complete: %s\n", file.c_str());

    // Auto-save .cub5 alongside the export (same base name)
    QString qfile = QString::fromStdString(file);
    int dot = qfile.lastIndexOf('.');
    if (dot > qfile.lastIndexOf('/')) {
      QString cub5 = qfile.left(dot) + ".cub5";
      std::string save_cmd = "save cub5 \"" + cub5.toStdString() + "\" overwrite";
      CubitInterface::cmd(save_cmd.c_str());
    }
  } else {
    PRINT_ERROR("Export failed: %s\n", file.c_str());
  }
}

void RadiaMenuHandler::export_gmsh()    { run_export(ExportDialog::GMSH); }
void RadiaMenuHandler::export_nastran() { run_export(ExportDialog::Nastran); }
void RadiaMenuHandler::export_vtk()     { run_export(ExportDialog::VTK); }
void RadiaMenuHandler::export_meg()     { run_export(ExportDialog::MEG); }

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

  // Show export dialog (shared with GMSH/Nastran/VTK/MEG)
  QString jouPath;
  std::string jf = CubitInterface::get_current_journal_file();
  if (!jf.empty()) jouPath = QString::fromStdString(jf);
  if (jouPath.isEmpty()) jouPath = s_lastJouPath;

  ExportDialog dlgOpt(ExportDialog::NETGEN_VOL, jouPath);
  if (dlgOpt.exec() != QDialog::Accepted) return;

  int order = dlgOpt.order();
  QString volPath = dlgOpt.filePath();
  volPath.replace("\\", "/");

  // --- Phase 1: C++ export .vol + companion JSON (no subprocess, no cub5) ---
  PRINT_INFO("Exporting Netgen .vol (order %d)...\n", order);
  std::string cmd = "export netgen \"" + volPath.toStdString()
    + "\" order " + std::to_string(order) + " overwrite";
  CubitInterface::silent_cmd(cmd.c_str());

  // Check .vol was created
  if (!QFile::exists(volPath)) {
    PRINT_ERROR("export netgen command failed.\n");
    return;
  }
  PRINT_INFO("Exported: %s\n", volPath.toStdString().c_str());

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
  if (proc.exitCode() != 0) {
    PRINT_WARNING("NGSolve verification failed:\n%s\n",
      proc.readAllStandardError().left(2000).constData());
    return;
  }

  // Parse JSON result from subprocess
  QString out = QString::fromUtf8(proc.readAllStandardOutput());
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
    auto *item = new QTableWidgetItem(QString::number(err, 'e', 2));
    if (std::abs(err) > 1.0)
      item->setBackground(QColor(255, 200, 200));
    return item;
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

  dlg->show();  // modeless: user can interact with Cubit
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

  // --- Phase 1: C++ export .vol for p=1..5 ---
  PRINT_INFO("Exporting .vol for p=1..5...\n");
  QString tmpDir = QDir::tempPath();
  tmpDir.replace("\\", "/");
  if (!tmpDir.endsWith('/')) tmpDir += '/';

  QStringList volPaths;
  for (int p = 1; p <= 5; p++) {
    QString vp = tmpDir + QString("radia_eval_p%1.vol").arg(p);
    std::string cmd = "export netgen \"" + vp.toStdString()
      + "\" order " + std::to_string(p) + " overwrite";
    PRINT_INFO("  p=%d: exporting...\n", p);
    CubitInterface::silent_cmd(cmd.c_str());
    if (QFile::exists(vp)) {
      volPaths << vp;
    } else {
      PRINT_WARNING("  p=%d: export failed, skipping.\n", p);
      volPaths << "";  // placeholder
    }
  }

  // --- Phase 2: subprocess to verify all .vol with NGSolve ---
  QString python = find_external_python();
  QString script = find_calc_script(python, "calc_mesh_eval");
  if (script.isEmpty()) {
    PRINT_ERROR("Cannot find calc_mesh_eval.py. Is radia installed?\n");
    return;
  }

  PRINT_INFO("Verifying with NGSolve...\n");

  QProcess proc;
  QStringList args;
  if (python == "py") args << "-3";
  args << script << "--vol";
  for (auto &vp : volPaths) {
    if (!vp.isEmpty()) args << vp;
  }

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

  // Cleanup temp .vol files
  for (auto &vp : volPaths) {
    if (!vp.isEmpty()) {
      QFile::remove(vp);
      QFile::remove(vp + ".json");
    }
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
  dlg.setMinimumSize(820, 320);
  QVBoxLayout *layout = new QVBoxLayout(&dlg);

  // CAD reference labels
  QLabel *cadLabel = new QLabel(
    QString("CAD Volume: %1 m^3     CAD Area: %2 m^2     CAD Length: %3 m")
      .arg(cad_v, 0, 'e', 6).arg(cad_a, 0, 'e', 6).arg(cad_l, 0, 'e', 6), &dlg);
  layout->addWidget(cadLabel);

  // Table: empty column after L err for visual separation from future columns
  QTableWidget *table = new QTableWidget(nrows, 8, &dlg);
  table->setHorizontalHeaderLabels(
    {"p", "Volume", "V err [%]", "Area", "A err [%]", "Length", "L err [%]", ""});
  table->horizontalHeader()->setStretchLastSection(true);
  table->setColumnWidth(7, 10);  // narrow separator column
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
      // col 7 = empty separator
    }
  }
  table->resizeColumnsToContents();
  layout->addWidget(table);

  // Buttons: OK + Copy
  QDialogButtonBox *buttons = new QDialogButtonBox(
    QDialogButtonBox::Ok | QDialogButtonBox::Save, &dlg);
  buttons->button(QDialogButtonBox::Save)->setText("Copy to Clipboard");
  layout->addWidget(buttons);

  QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
  QObject::connect(buttons->button(QDialogButtonBox::Save), &QPushButton::clicked,
    [&]() {
      QString tsv = "p\tVolume\tV err [%]\tArea\tA err [%]\tLength\tL err [%]\n";
      for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < 7; j++) {  // skip col 7 (empty separator)
          if (j > 0) tsv += "\t";
          QTableWidgetItem *item = table->item(i, j);
          tsv += item ? item->text() : "";
        }
        tsv += "\n";
      }
      QApplication::clipboard()->setText(tsv);
    });

  dlg.exec();
}

// ============================================================
// ExportDialog - format-specific options
// ============================================================

// --- Settings persistence ---
static QString settingsPath()
{
  return QDir::homePath() + "/.cubit/radia_export_settings.json";
}

static QJsonObject loadSettings()
{
  QFile f(settingsPath());
  if (!f.open(QIODevice::ReadOnly)) return {};
  return QJsonDocument::fromJson(f.readAll()).object();
}

static void saveSettings(const QJsonObject &obj)
{
  QDir().mkpath(QDir::homePath() + "/.cubit");
  QFile f(settingsPath());
  if (f.open(QIODevice::WriteOnly))
    f.write(QJsonDocument(obj).toJson(QJsonDocument::Indented));
}

ExportDialog::ExportDialog(Format format, const QString &jouPath, QWidget* parent)
  : QDialog(parent), mFormat(format),
    mVersion(nullptr), mNoPyramid(nullptr)
{
  // Window title
  const char* titles[] = {"Export Netgen Vol", "Export GMSH", "Export Nastran BDF", "Export VTK", "Export MEG"};
  setWindowTitle(titles[format]);
  setMinimumWidth(500);

  QVBoxLayout* layout = new QVBoxLayout(this);
  QFormLayout* form = new QFormLayout();

  // Determine default directory and base name
  const char* exts[] = {".vol", ".msh", ".bdf", ".vtk", ".meg"};
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
    char cwd[1024];
    if (_getcwd(cwd, sizeof(cwd))) {
      defaultDir = QString::fromLocal8Bit(cwd);
      defaultDir.replace("\\", "/");
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
  if (format != MEG) {
    form->addRow("Order:", mOrderCombo);
  }

  // Version: GMSH is always v2.2 (v4.1 is for post-processing only)

  // Dimension (not for NETGEN_VOL)
  mDimension = new QComboBox();
  if (format == MEG) {
    mDimension->addItems({"T (3D)", "K (2D)", "R (Axisymmetric)"});
  } else {
    mDimension->addItems({"3D", "2D"});
  }
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
    const char* keys[] = {"netgen_vol", "gmsh", "nastran", "vtk", "meg"};
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
    const char* keys[] = {"netgen_vol", "gmsh", "nastran", "vtk", "meg"};
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
      cmd = QString("export nastran \"%1\" order %2 dimension %3")
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
    case MEG: {
      // Map combo selection to keyword
      QChar c = mDimension->currentText().at(0);
      QString kw = (c == 'K') ? "twod" : (c == 'R') ? "axisymmetric" : "threed";
      cmd = QString("export meg \"%1\" %2").arg(file).arg(kw);
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
