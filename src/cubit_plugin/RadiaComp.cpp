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
#include <QTextEdit>
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

  QAction* a_netgen = new QAction("Netgen Vol + Pkl...", handler);
  a_netgen->setStatusTip("Export mesh as Netgen .vol (linear) + .pkl (curved)");
  QObject::connect(a_netgen, SIGNAL(triggered()), handler, SLOT(export_netgen()));
  menu_list.push_back(a_netgen);

  gui->add_to_menu("&Export Mesh", menu_list, "radiacomp");

  // --- Mesh Evaluation menu ---
  std::vector<QAction*> eval_list;

  QAction* a_vol = new QAction("Volume / Surface Area...", handler);
  a_vol->setStatusTip("Compute mesh volume and surface area using NGSolve integration");
  QObject::connect(a_vol, SIGNAL(triggered()), handler, SLOT(mesh_volume()));
  eval_list.push_back(a_vol);

  gui->add_to_menu("Mesh &Evaluation", eval_list, "radiacomp");

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

  // Ask for order
  bool ok;
  int order = QInputDialog::getInt(nullptr, "Netgen Export",
    "Curve order (2-5):", 3, 2, 5, 1, &ok);
  if (!ok) return;

  // Get output path
  QString jouPath;
  std::string jf = CubitInterface::get_current_journal_file();
  if (!jf.empty()) jouPath = QString::fromStdString(jf);
  if (jouPath.isEmpty()) jouPath = s_lastJouPath;

  QString defaultDir = jouPath.isEmpty()
    ? QDir::currentPath() : QFileInfo(jouPath).absolutePath();
  QString baseName = jouPath.isEmpty()
    ? "ExportedMesh" : QFileInfo(jouPath).completeBaseName();

  QString volPath = defaultDir + "/" + baseName + ".vol";
  QString pklPath = defaultDir + "/" + baseName + "_curved.pkl";
  volPath.replace("\\", "/");
  pklPath.replace("\\", "/");

  // Save temp .cub5
  QString cub5 = QDir::tempPath() + "/radia_netgen_export.cub5";
  cub5.replace("\\", "/");
  CubitInterface::cmd(("save cub5 \"" + cub5.toStdString() + "\" overwrite").c_str());

  QString python = find_external_python();
  QString script = find_calc_script(python, "calc_export_netgen");
  if (script.isEmpty()) {
    PRINT_ERROR("Cannot find calc_export_netgen.py.\n");
    return;
  }

  PRINT_INFO("Exporting Netgen .vol + .pkl (order %d)...\n", order);

  QProcess proc;
  QStringList args;
  if (python == "py") args << "-3";
  args << script << "--cub5" << cub5 << "--order" << QString::number(order)
       << "--vol" << volPath << "--pkl" << pklPath;

  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  for (const QString &key : env.keys()) {
    if (key.toUpper().contains("QT") || key.toUpper().contains("PYSIDE"))
      env.remove(key);
  }
  proc.setProcessEnvironment(env);

  proc.start(python, args);
  if (!proc.waitForFinished(300000)) {
    PRINT_ERROR("Netgen export timed out.\n");
    return;
  }
  if (proc.exitCode() != 0) {
    PRINT_ERROR("Netgen export failed:\n%s\n",
      proc.readAllStandardError().left(2000).constData());
    return;
  }

  PRINT_INFO("Exported:\n  %s (linear)\n  %s (curved, order %d)\n",
    volPath.toStdString().c_str(), pklPath.toStdString().c_str(), order);
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

  // Save to temp .cub5
  QString cub5 = QDir::tempPath() + "/radia_mesh_eval.cub5";
  cub5.replace("\\", "/");
  CubitInterface::cmd(("save cub5 \"" + cub5.toStdString() + "\" overwrite").c_str());

  QString python = find_external_python();
  QString script = find_calc_script(python, "calc_mesh_eval");
  if (script.isEmpty()) {
    PRINT_ERROR("Cannot find calc_mesh_eval.py. Is radia installed?\n");
    return;
  }

  // Run calc_mesh_eval.py --cub5 ... --max-order 5
  PRINT_INFO("Running mesh evaluation (p=1..5)...\n");

  QProcess proc;
  QStringList args;
  if (python == "py") args << "-3";
  args << script << "--cub5" << cub5 << "--max-order" << "5";

  // Clean Qt env to avoid conflicts
  QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
  for (const QString &key : env.keys()) {
    if (key.toUpper().contains("QT") || key.toUpper().contains("PYSIDE"))
      env.remove(key);
  }
  proc.setProcessEnvironment(env);

  proc.start(python, args);
  if (!proc.waitForFinished(600000)) {  // 10 min
    PRINT_ERROR("calc_mesh_eval.py timed out.\n");
    return;
  }
  if (proc.exitCode() != 0) {
    PRINT_ERROR("calc_mesh_eval.py failed:\n%s\n",
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

  // Build TSV table (tab-separated, copy-pasteable to Excel)
  QString tsv;
  tsv += QString("CAD Volume\t%1\tm^3\n").arg(cad_v, 0, 'e', 6);
  tsv += QString("CAD Area\t%1\tm^2\n\n").arg(cad_a, 0, 'e', 6);
  tsv += "p\tVolume\tV err [%]\tArea\tA err [%]\n";

  QJsonArray orders = r["orders"].toArray();
  for (int i = 0; i < orders.size(); i++) {
    QJsonObject o = orders[i].toObject();
    if (o.contains("error")) {
      tsv += QString("%1\t%2\n").arg(o["order"].toInt()).arg(o["error"].toString());
      continue;
    }
    tsv += QString("%1\t%2\t%3\t%4\t%5\n")
      .arg(o["order"].toInt())
      .arg(o["ng_volume"].toDouble(), 0, 'e', 6)
      .arg(o["vol_error_pct"].toDouble(), 0, 'f', 5)
      .arg(o["ng_area"].toDouble(), 0, 'e', 6)
      .arg(o["area_error_pct"].toDouble(), 0, 'f', 5);
  }

  // Show in a dialog with selectable/copyable QTextEdit
  QDialog dlg;
  dlg.setWindowTitle("Mesh Evaluation - Volume / Surface Area");
  dlg.setMinimumSize(550, 300);
  QVBoxLayout *layout = new QVBoxLayout(&dlg);

  QTextEdit *text = new QTextEdit(&dlg);
  text->setReadOnly(true);
  text->setFont(QFont("Consolas", 10));
  text->setPlainText(tsv);
  layout->addWidget(text);

  QDialogButtonBox *buttons = new QDialogButtonBox(
    QDialogButtonBox::Ok | QDialogButtonBox::Save, &dlg);
  buttons->button(QDialogButtonBox::Save)->setText("Copy to Clipboard");
  layout->addWidget(buttons);

  QObject::connect(buttons, &QDialogButtonBox::accepted, &dlg, &QDialog::accept);
  QObject::connect(buttons->button(QDialogButtonBox::Save), &QPushButton::clicked,
    [&tsv]() {
      QApplication::clipboard()->setText(tsv);
    });

  dlg.exec();
}

// ============================================================
// ExportDialog - format-specific options
// ============================================================

ExportDialog::ExportDialog(Format format, const QString &jouPath, QWidget* parent)
  : QDialog(parent), mFormat(format),
    mVersion(nullptr), mNoPyramid(nullptr)
{
  // Window title
  const char* titles[] = {"Export GMSH", "Export Nastran BDF", "Export VTK", "Export MEG"};
  setWindowTitle(titles[format]);
  setMinimumWidth(500);

  QVBoxLayout* layout = new QVBoxLayout(this);
  QFormLayout* form = new QFormLayout();

  // Determine default directory and base name
  const char* exts[] = {".msh", ".bdf", ".vtk", ".meg"};
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

  // Order (not for MEG)
  mOrderCombo = new QComboBox();
  mOrderCombo->addItems({"1 (linear)", "2 (quadratic)"});
  connect(mOrderCombo, SIGNAL(currentIndexChanged(int)), this, SLOT(updatePreview()));
  if (format != MEG) {
    form->addRow("Order:", mOrderCombo);
  }

  // Version (GMSH only)
  if (format == GMSH) {
    mVersion = new QComboBox();
    mVersion->addItems({"2.2", "4.1"});
    connect(mVersion, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
    form->addRow("Version:", mVersion);
  }

  // Dimension
  mDimension = new QComboBox();
  if (format == MEG) {
    mDimension->addItems({"T (3D)", "K (2D)", "R (Axisymmetric)"});
  } else {
    mDimension->addItems({"3D", "2D"});
  }
  connect(mDimension, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
  form->addRow("Dimension:", mDimension);

  // NoPyramid (Nastran only)
  if (format == Nastran) {
    mNoPyramid = new QComboBox();
    mNoPyramid->addItems({"Keep pyramids", "Convert to degenerate hex (JMAG)"});
    connect(mNoPyramid, SIGNAL(currentTextChanged(QString)), this, SLOT(updatePreview()));
    form->addRow("Pyramids:", mNoPyramid);
  }

  layout->addLayout(form);

  // Command preview
  QFormLayout* previewForm = new QFormLayout();
  mPreview = new QLineEdit();
  mPreview->setReadOnly(true);
  previewForm->addRow("Command:", mPreview);
  layout->addLayout(previewForm);

  // OK / Cancel
  QDialogButtonBox* buttons = new QDialogButtonBox(
      QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
  connect(buttons, SIGNAL(accepted()), this, SLOT(accept()));
  connect(buttons, SIGNAL(rejected()), this, SLOT(reject()));
  layout->addWidget(buttons);

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
    case GMSH: {
      QString ver = (mVersion && mVersion->currentText() == "4.1") ? "4" : "2";
      QString dim = (mDimension->currentText() == "2D") ? "2" : "3";
      cmd = QString("export gmsh \"%1\" order %2 version %3 dimension %4")
                .arg(file).arg(order).arg(ver).arg(dim);
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
}
