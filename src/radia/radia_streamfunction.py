"""radia_streamfunction panel -- Stream-Function coil design (generator-driven).

Layer-3 PySide6 panel for calc_streamfunction.py.  A composite top-level panel
holds the COIL surface .vol Browse + a Method combo + a QStackedWidget over the
three sub-panels (Design / Pareto / Manufacture), each a ModePanel that
bind_argparser()s the SAME calc_streamfunction.build_argparser() with a
mode-specific skip set.  Per the new-panel contract the panel, the subprocess
wiring, the .log capture, the Open-GMSH button and the JSON result are all
derived from that one argparser (no hand-assembled argv, no orphan widgets).

  primary input  : coil surface .vol      -> --coil-vol  (the wp_vol Browse)
  per-mode inputs: eval region .vol (Browse) + target CoefficientFunction expr
                   + the mode's solver / pareto / manufacture knobs
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PySide6.QtWidgets import QStackedWidget

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, msh_output, json_output, run_app,
)

TITLE = "Radia - Stream-Function Coil Design"

_EVAL_BROWSE = {"eval_vol": ("Eval region .vol",
                             "Netgen Vol (*.vol);;All (*)")}
_REG_LABELS = {"regularize": "Regularisation:", "confine": "Current confine:"}
_REG_CHOICES = {"regularize": [("L2 (min |psi|)", "l2"),
                               ("H1 (min current energy)", "h1")],
                "confine": [("Abe edge-equipotential (recommended)", "abe"),
                            ("Unconfined", "off"),
                            ("psi=0 on former edge", "on")]}

# argparser dests that the top panel / build_command own, never a sub widget.
_BASE_SKIP = ("coil_vol", "output", "msh_output", "method")
# mode-specific dest groups (so each sub-panel hides the other modes' knobs)
_PARETO_DESTS = ("pareto_lever", "alpha_min", "alpha_max", "n_alpha",
                 "linf_iter", "geom_scale_min", "geom_scale_max")
_MANUFACTURE_DESTS = ("nlevels", "contour_sub", "chain", "chain_ncut",
                      "chain_passes", "distort", "distort_grid", "distort_iter",
                      "step_output", "peec", "wire_diam", "peec_freq",
                      "flux_plot", "flux_plane", "steps_plot",
                      "target_inductance", "resonance_cap", "nlevels_max")


def _argparser():
    from radia.panels.calc_streamfunction import build_argparser
    return build_argparser()


def _volume_argparser():
    from radia.panels.calc_streamfunction_volume import build_argparser
    return build_argparser()


# Volume-3D mode (foliated-Clebsch volume stream function) uses its OWN calc
# script + argparser, so its skip/label sets are independent of the surface
# modes above.  Routed inputs: --vol via wp_vol; --output / --msh-output via
# build_command; --foliation / --n-threads kept at calc defaults.
_VOLUME_SKIP = ("vol", "output", "msh_output", "foliation", "n_threads")
_VOLUME_LABELS = {"target_bz": "Target Bz [T]:", "n_leaves": "mu-leaves:",
                  "fes_order": "lambda H1 order:", "aca_eps": "ACA+ eps:",
                  "rings_per_span": "Rings / lambda span:",
                  "n_targets": "Axial targets:",
                  "nphi": "Contour nphi:", "nz": "Contour nz:"}


class _DesignPanel(ModePanel):
    """target -> psi -> homogeneity + peak current density."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bind_argparser(
            _argparser(),
            skip=_BASE_SKIP + _PARETO_DESTS + _MANUFACTURE_DESTS,
            file_browse=_EVAL_BROWSE,
            choice_map=_REG_CHOICES, labels=_REG_LABELS,
        )

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path, vol_flag="--coil-vol",
            script_path=calc_script("calc_streamfunction.py"),
            output_path=json_output(vol_path, "_sf_design"),
            extra=["--method", "design",
                   "--msh-output", msh_output(vol_path, "_sf_design")],
        )


class _ParetoPanel(ModePanel):
    """alpha-sweep -> (homogeneity, peak current density) Pareto front."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bind_argparser(
            _argparser(),
            skip=_BASE_SKIP + _MANUFACTURE_DESTS,
            file_browse=_EVAL_BROWSE,
            choice_map=_REG_CHOICES, labels=_REG_LABELS,
        )

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path, vol_flag="--coil-vol",
            script_path=calc_script("calc_streamfunction.py"),
            output_path=json_output(vol_path, "_sf_pareto"),
            extra=["--method", "pareto"],
        )


class _ManufacturePanel(ModePanel):
    """single-stroke chain + sheet-metal distortion + CAD(STEP) + PEEC."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bind_argparser(
            _argparser(),
            skip=_BASE_SKIP + _PARETO_DESTS,
            file_browse=_EVAL_BROWSE,
            choice_map=_REG_CHOICES, labels=_REG_LABELS,
        )

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path, vol_flag="--coil-vol",
            script_path=calc_script("calc_streamfunction.py"),
            output_path=json_output(vol_path, "_sf_manufacture"),
            extra=["--method", "manufacture",
                   "--msh-output", msh_output(vol_path, "_sf_manufacture")],
        )


class _Volume3DPanel(ModePanel):
    """3D VOLUME stream function (foliated Clebsch J=grad(lambda)xgrad(mu)):
    conductor tube .vol + target axial Bz -> equal-current windable wires."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.bind_argparser(
            _volume_argparser(),
            skip=_VOLUME_SKIP,
            labels=_VOLUME_LABELS,
        )

    def build_command(self, vol_path):
        return self.build_command_from_parser(
            vol_path=vol_path, vol_flag="--vol",
            script_path=calc_script("calc_streamfunction_volume.py"),
            output_path=json_output(vol_path, "_sf_volume"),
            extra=["--msh-output", msh_output(vol_path, "_sf_volume")],
        )


_METHOD_FACTORIES = [
    ("Design", _DesignPanel),
    ("Pareto", _ParetoPanel),
    ("Manufacture", _ManufacturePanel),
    ("Volume 3D", _Volume3DPanel),
]


class StreamFunctionPanel(ModePanel):
    """Composite: coil .vol Browse + Method combo + QStackedWidget."""
    def __init__(self, parent=None):
        super().__init__(parent)
        # primary input = the COIL surface .vol (AnalysisWindow.wp_vol contract)
        # surface modes -> coil SURFACE .vol; Volume 3D -> conductor VOLUME .vol
        self.add_browse("wp_vol", "Coil / conductor .vol:",
                        filter_str="Netgen Vol (*.vol);;All (*)")
        self._method_combo = self.add_combo(
            "method", "Method:", [name for name, _ in _METHOD_FACTORIES])
        self._sub_panels = {name: factory()
                            for name, factory in _METHOD_FACTORIES}
        self._stack = QStackedWidget()
        for name, _ in _METHOD_FACTORIES:
            self._stack.addWidget(self._sub_panels[name])
        self._form.addRow(self._stack)
        self._method_combo.currentTextChanged.connect(self._on_method_changed)
        self._on_method_changed(_METHOD_FACTORIES[0][0])

    def _current_sub(self):
        return self._sub_panels[self._method_combo.currentText()]

    def _on_method_changed(self, text):
        if text in self._sub_panels:
            self._stack.setCurrentWidget(self._sub_panels[text])
        cb = getattr(self, "validationChanged", None)
        if callable(cb):
            cb()

    def is_runnable(self):
        return self._current_sub().is_runnable()

    def build_command(self, vol_path):
        return self._current_sub().build_command(vol_path)

    def wp_vol_path(self):
        return self.val("wp_vol") if "wp_vol" in self._widgets else ""

    def save_state(self):
        state = super().save_state()
        for name, sub in self._sub_panels.items():
            for k, v in sub.save_state().items():
                state[f"{name}/{k}"] = v
        return state

    def restore_state(self, state):
        if not state:
            return
        super().restore_state({k: v for k, v in state.items() if "/" not in k})
        for name, sub in self._sub_panels.items():
            prefix = f"{name}/"
            sub.restore_state({k[len(prefix):]: v for k, v in state.items()
                               if k.startswith(prefix)})
        self._on_method_changed(self._method_combo.currentText())


class StreamFunctionWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__(TITLE, vol_path, settings_key="streamfunction")
        panel = StreamFunctionPanel()
        self._set_panel(panel)
        # Restore last session's widget state FIRST, then let an explicit
        # launcher vol_path override the restored coil wp_vol.  Restoring
        # AFTER setText (or not restoring at all) would silently drop the
        # Cubit Solve-menu's current .vol -- the radia_em / radia_ih
        # restore-order bug class.
        self._restore_settings()
        if vol_path and "wp_vol" in panel._widgets:
            panel._widgets["wp_vol"].setText(self.display_path(vol_path))


if __name__ == "__main__":
    run_app(StreamFunctionWindow, sys.argv)
