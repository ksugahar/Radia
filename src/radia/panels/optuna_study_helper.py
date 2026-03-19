"""
Optuna-compatible SQLite study storage for Cubit panel computations.

Each "Solve" click in the panel = 1 Optuna trial.
Results are stored in SQLite next to the cub5 file,
viewable via optuna-dashboard.

Usage:
    from optuna_study_helper import InductanceStudy

    study = InductanceStudy("model.cub5")
    study.record_trial(
        params={"curve_order": 2, "conductivity": 5.8e7, ...},
        result={"inductance_H": 3.7e-7, "n_free_dofs": 269, ...},
    )
    study.launch_dashboard()  # opens browser

Does NOT depend on Qt. Safe to import from both GUI and CLI.
"""

import os
import sys
from datetime import datetime

try:
	import optuna
	from optuna.distributions import (
		IntDistribution,
		FloatDistribution,
		CategoricalDistribution,
	)
	from optuna.trial import FrozenTrial, TrialState
	HAS_OPTUNA = True
except ImportError:
	HAS_OPTUNA = False


# Parameter definitions (distributions for Optuna compatibility)
PARAM_DISTRIBUTIONS = {
	"curve_order": lambda v: IntDistribution(1, 5),
	"conductivity": lambda v: FloatDistribution(1e-2, 1e9, log=True),
	"frequency": lambda v: FloatDistribution(0.0, 1e12),
	"source_block": lambda v: CategoricalDistribution([v]),
	"sink_block": lambda v: CategoricalDistribution([v]),
}


def _db_path_from_cub5(cub5_path):
	"""Derive SQLite path from cub5 file path.

	Example: /path/to/inductor.cub5 -> /path/to/inductor_optuna.db
	"""
	base = os.path.splitext(os.path.abspath(cub5_path))[0]
	return base + "_optuna.db"


def _study_name_from_cub5(cub5_path):
	"""Derive study name from cub5 filename.

	Example: inductor.cub5 -> "inductor"
	"""
	return os.path.splitext(os.path.basename(cub5_path))[0]


class InductanceStudy:
	"""Optuna study wrapper for inductance extraction trials."""

	def __init__(self, cub5_path):
		"""Create or load study for the given cub5 model.

		Args:
			cub5_path: Path to .cub5 file. SQLite DB is placed alongside.
		"""
		if not HAS_OPTUNA:
			self._study = None
			self._db_path = None
			return

		self._db_path = _db_path_from_cub5(cub5_path)
		self._cub5_path = os.path.abspath(cub5_path)
		study_name = _study_name_from_cub5(cub5_path)
		storage = f"sqlite:///{self._db_path}"

		# Suppress Optuna info logs
		optuna.logging.set_verbosity(optuna.logging.WARNING)

		self._study = optuna.create_study(
			study_name=study_name,
			storage=storage,
			direction="maximize",  # higher inductance = stronger coupling
			load_if_exists=True,
		)

	@property
	def available(self):
		"""True if Optuna is installed and study is loaded."""
		return self._study is not None

	@property
	def db_path(self):
		return self._db_path

	@property
	def n_trials(self):
		if not self.available:
			return 0
		return len(self._study.trials)

	def record_trial(self, params, result, start_time=None):
		"""Record a completed computation as an Optuna trial.

		Args:
			params: dict with keys matching PARAM_DISTRIBUTIONS
				(curve_order, conductivity, frequency, source_block, sink_block)
			result: dict from calc_inductance.py JSON output
				(inductance_H, n_free_dofs, neg_diag, surface_area, gmsh_file, ...)
			start_time: datetime when computation started (optional)

		Returns:
			Trial number, or -1 if Optuna not available.
		"""
		if not self.available:
			return -1

		# Build distributions to match params
		distributions = {}
		for key, value in params.items():
			if key in PARAM_DISTRIBUTIONS:
				distributions[key] = PARAM_DISTRIBUTIONS[key](value)

		# Build user_attrs (metadata not in param space)
		inductance_H = result.get("inductance_H", 0.0)
		user_attrs = {
			"n_free_dofs": result.get("n_free_dofs", 0),
			"neg_diag": result.get("neg_diag", 0),
			"surface_area": result.get("surface_area", 0.0),
			"gmsh_file": result.get("gmsh_file", ""),
			"cub5_file": self._cub5_path,
			"inductance_display": _format_inductance(inductance_H),
		}

		now = datetime.now()
		trial = FrozenTrial(
			number=0,  # auto-assigned by add_trial
			state=TrialState.COMPLETE,
			values=[inductance_H],
			datetime_start=start_time or now,
			datetime_complete=now,
			params=params,
			distributions=distributions,
			user_attrs=user_attrs,
			system_attrs={},
			intermediate_values={},
			trial_id=0,  # auto-assigned
			value=None,
		)

		self._study.add_trial(trial)
		return len(self._study.trials) - 1

	def get_trials_dataframe(self):
		"""Return trials as pandas DataFrame (for display)."""
		if not self.available:
			return None
		try:
			return self._study.trials_dataframe()
		except Exception:
			return None

	def launch_dashboard(self, port=8080):
		"""Launch optuna-dashboard as background process.

		Returns:
			URL string, or None if failed.
		"""
		if not self._db_path or not os.path.exists(self._db_path):
			return None

		import subprocess as sp
		import shutil

		storage_url = f"sqlite:///{self._db_path}"

		# Try optuna-dashboard executable, then module
		exe = shutil.which("optuna-dashboard")
		if exe:
			cmd = [exe, storage_url, "--port", str(port)]
		else:
			cmd = [sys.executable, "-m", "optuna_dashboard",
			       storage_url, "--port", str(port)]

		try:
			# DETACHED_PROCESS on Windows
			creationflags = 0x00000008 if sys.platform == "win32" else 0
			sp.Popen(cmd, creationflags=creationflags)
			return f"http://localhost:{port}"
		except Exception:
			return None


def _format_inductance(L):
	"""Format inductance value as human-readable string."""
	L = abs(L)
	if L >= 1e-3:
		return f"{L*1e3:.4f} mH"
	elif L >= 1e-6:
		return f"{L*1e6:.4f} uH"
	elif L >= 1e-9:
		return f"{L*1e9:.4f} nH"
	elif L >= 1e-12:
		return f"{L*1e12:.4f} pH"
	else:
		return f"{L:.4e} H"
