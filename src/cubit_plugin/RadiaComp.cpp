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
#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
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
// Mesh Volume / Surface Area — subprocess to calc_volume.py
// ============================================================
#include <QMessageBox>
#include <QInputDialog>
#include <QProcess>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QTemporaryFile>

void RadiaMenuHandler::mesh_volume()
{
  if (CubitInterface::get_volume_count() == 0) {
    ensure_model();
    if (CubitInterface::get_volume_count() == 0)
      return;
  }

  // Ask for curve order
  bool ok;
  int order = QInputDialog::getInt(nullptr, "Mesh Evaluation",
    "Curve order (1=linear, 2-5=high-order):", 3, 1, 5, 1, &ok);
  if (!ok) return;

  // Save to temporary .cub5
  QString tmpDir = QDir::tempPath();
  QString cub5 = tmpDir + "/radia_mesh_eval.cub5";
  cub5.replace("\\", "/");
  std::string save_cmd = "save cub5 \"" + cub5.toStdString() + "\" overwrite";
  CubitInterface::cmd(save_cmd.c_str());

  // Find external Python (system Python 3.12, not Cubit's 3.10)
  QString python;
#ifdef _WIN32
  // Try RADIA_PYTHON env, then py -3, then python
  QByteArray env = qgetenv("RADIA_PYTHON");
  if (!env.isEmpty()) {
    python = QString::fromLocal8Bit(env);
  } else {
    // Try "py -3" first (Windows launcher)
    QProcess probe;
    probe.start("py", {"-3", "-c", "print('ok')"});
    if (probe.waitForFinished(5000) && probe.exitCode() == 0) {
      python = "py";
    } else {
      python = "python";
    }
  }
#else
  python = "python3";
#endif

  // Find calc_volume.py (installed via pip in radia/panels/)
  QString script;
  QProcess findScript;
  QStringList findArgs;
  if (python == "py")
    findArgs << "-3";
  findArgs << "-c" << "import radia.panels.calc_volume as m; import os; print(os.path.abspath(m.__file__))";
  findScript.start(python, findArgs);
  if (findScript.waitForFinished(10000) && findScript.exitCode() == 0) {
    script = QString::fromUtf8(findScript.readAllStandardOutput()).trimmed();
  }
  if (script.isEmpty() || !QFile::exists(script)) {
    QMessageBox::critical(nullptr, "Error",
      "Cannot find calc_volume.py. Is radia installed?\n"
      "  pip install radia");
    return;
  }

  // Run calc_volume.py as subprocess
  PRINT_INFO("Running: %s calc_volume.py --cub5 %s --order %d\n",
             python.toStdString().c_str(), cub5.toStdString().c_str(), order);

  QProcess proc;
  QStringList args;
  if (python == "py")
    args << "-3";
  args << script << "--cub5" << cub5 << "--order" << QString::number(order);
  proc.start(python, args);
  if (!proc.waitForFinished(300000)) {  // 5 min timeout
    QMessageBox::critical(nullptr, "Error", "calc_volume.py timed out (5 min).");
    return;
  }

  if (proc.exitCode() != 0) {
    QString err = QString::fromUtf8(proc.readAllStandardError());
    QMessageBox::critical(nullptr, "Error",
      "calc_volume.py failed:\n" + err.left(2000));
    return;
  }

  // Parse JSON from last line of stdout
  QString out = QString::fromUtf8(proc.readAllStandardOutput());
  QStringList lines = out.split('\n', Qt::SkipEmptyParts);
  if (lines.isEmpty()) {
    QMessageBox::critical(nullptr, "Error", "No output from calc_volume.py");
    return;
  }

  QJsonDocument doc = QJsonDocument::fromJson(lines.last().toUtf8());
  if (doc.isNull() || !doc.isObject()) {
    QMessageBox::critical(nullptr, "Error",
      "Invalid JSON from calc_volume.py:\n" + lines.last().left(500));
    return;
  }

  QJsonObject result = doc.object();
  if (result.contains("error")) {
    QMessageBox::warning(nullptr, "Mesh Evaluation",
      "Error: " + result["error"].toString());
    return;
  }

  // Build result message
  QString msg;
  double cad_total = result["cad_total"].toDouble();
  double ng_total = result["ngsolve_total"].toDouble();
  double ng_area = result["ngsolve_area"].toDouble();
  double err_pct = (cad_total > 0) ? (ng_total - cad_total) / cad_total * 100.0 : 0;

  msg += QString("Curve order: %1\n\n").arg(order);
  msg += QString("CAD Volume:     %1\n").arg(cad_total, 0, 'e', 6);
  msg += QString("NGSolve Volume: %1\n").arg(ng_total, 0, 'e', 6);
  msg += QString("Volume Error:   %1%\n\n").arg(err_pct, 0, 'f', 4);
  msg += QString("Surface Area:   %1\n").arg(ng_area, 0, 'e', 6);

  QJsonArray vols = result["volumes"].toArray();
  if (vols.size() > 1) {
    msg += "\nPer volume:\n";
    for (int i = 0; i < vols.size(); i++) {
      QJsonObject v = vols[i].toObject();
      msg += QString("  %1: CAD=%2  NGSolve=%3\n")
        .arg(v["name"].toString())
        .arg(v["cad_volume"].toDouble(), 0, 'e', 4)
        .arg(v.contains("ngsolve_volume") ? v["ngsolve_volume"].toDouble() : 0, 0, 'e', 4);
    }
  }

  QMessageBox::information(nullptr, "Mesh Evaluation - Volume", msg);
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
