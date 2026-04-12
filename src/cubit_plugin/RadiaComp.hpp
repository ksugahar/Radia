#ifndef RADIACOMP_HPP
#define RADIACOMP_HPP

#include "Component.hpp"
#include <QObject>
#include <QDialog>
#include <QComboBox>
#include <QLineEdit>
#include <QTableWidget>

class QAction;

//! Radia Component for Coreform Cubit GUI.
//!
//! Adds "Export Mesh" menu with GMSH/Nastran/VTK export.
//! .ccl auto-loaded from bin/ — no -commandplugindir needed.
class RadiaComp : public Component
{
public:
  RadiaComp();
  ~RadiaComp() override;

  void start_up(int withGUI) override;
  void clean_up() override;

private:
  void setup_menus();
  void cleanup_menus();
  bool mMenuInitialized;
};

//! QObject-based helper for menu signal/slot connections.
class RadiaMenuHandler : public QObject
{
  Q_OBJECT

public:
  RadiaMenuHandler(QObject* parent = nullptr) : QObject(parent) {}

public slots:
  void export_gmsh();
  void export_nastran();
  void export_vtk();
  void export_netgen();
  void export_femeem();
  void export_meg();
  void mesh_volume();
  void launch_radia_ngsolve();
};

//! Export options dialog.
//!
//! Directory: browse button (default: CWD)
//! Filename: text input (default: journal basename or "ExportedMesh")
class ExportDialog : public QDialog
{
  Q_OBJECT

public:
  enum Format { NETGEN_VOL, GMSH, Nastran, VTK, FEMEEM, MEG };

  ExportDialog(Format format, const QString &jouPath = QString(),
               QWidget* parent = nullptr);

  //! Full output path (dir + filename)
  QString filePath() const;

  //! Preview command string
  QString cubitCommand() const;

  //! Accessors
  int order() const { return mOrderCombo->currentIndex() + 1; }
  int gmshVersion() const { return 2; }  // always v2.2 (v4.1 is post-processing only)
  int dimension() const { return (mDimension->currentText() == "2D") ? 2 : 3; }
  bool noPyramid() const { return mNoPyramid && mNoPyramid->currentIndex() == 1; }
  double scale() const;

private slots:
  void browseDir();
  void updatePreview();

private:
  Format     mFormat;
  QLineEdit* mDir;       // output directory
  QLineEdit* mFileName;  // filename (no path)
  QComboBox* mOrderCombo;
  QComboBox* mVersion;    // GMSH only
  QComboBox* mDimension;
  QComboBox* mNoPyramid;  // Nastran only
  QLineEdit* mScale;      // FEMEEM only
  QTableWidget* mBlockTable;  // MEG only: per-block ELF label
  QLineEdit* mPreview;

  void populateBlockTable();
  void updateBlockLabels();  // called when DIM changes
};

#endif // RADIACOMP_HPP
