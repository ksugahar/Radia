"""
Test: Optuna study storage for inductance extraction.

Verifies that trials are recorded correctly and optuna-dashboard
can read the SQLite database.
"""
import os
import sys
import tempfile

# Setup path
radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'src', 'radia', 'panels')
sys.path.insert(0, os.path.abspath(radia_src))

from optuna_study_helper import InductanceStudy, HAS_OPTUNA

if not HAS_OPTUNA:
	print("SKIP: optuna not installed")
	sys.exit(0)

# Create temp cub5 path (no actual file needed for study test)
tmpdir = tempfile.mkdtemp(prefix="radia_optuna_test_")
cub5_path = os.path.join(tmpdir, "test_inductor.cub5")
with open(cub5_path, "w") as f:
	f.write("")  # dummy

# Create study
study = InductanceStudy(cub5_path)
assert study.available, "Study should be available"
assert study.n_trials == 0, "Fresh study should have 0 trials"

db_path = study.db_path
assert db_path.endswith("_optuna.db"), f"Unexpected db path: {db_path}"
print(f"DB: {db_path}")

# Record trial 1: order=1
params1 = {
	"curve_order": 1,
	"conductivity": 5.8e7,
	"frequency": 0.0,
	"source_block": "conductor",
	"sink_block": "(closed loop)",
}
result1 = {
	"inductance_H": 3.691e-7,
	"n_free_dofs": 269,
	"neg_diag": 0,
	"surface_area": 1.97e-3,
	"gmsh_file": "/tmp/test.msh",
}
t1 = study.record_trial(params1, result1)
assert t1 == 0, f"First trial should be #0, got #{t1}"
assert study.n_trials == 1

# Record trial 2: order=2
params2 = {**params1, "curve_order": 2}
result2 = {**result1, "inductance_H": 3.45e-7, "n_free_dofs": 1073}
t2 = study.record_trial(params2, result2)
assert t2 == 1, f"Second trial should be #1, got #{t2}"
assert study.n_trials == 2

# Record trial 3: different frequency
params3 = {**params1, "frequency": 1e6}
result3 = {**result1, "inductance_H": 3.50e-7}
t3 = study.record_trial(params3, result3)
assert study.n_trials == 3

# Reload study from same DB (persistence check)
study2 = InductanceStudy(cub5_path)
assert study2.n_trials == 3, f"Reloaded study should have 3 trials, got {study2.n_trials}"

# Check DataFrame
df = study2.get_trials_dataframe()
if df is not None:
	print(f"\nTrials DataFrame ({len(df)} rows):")
	print(df.to_string())
else:
	print("DataFrame not available (pandas missing?)")

# Verify via raw Optuna API
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
raw_study = optuna.load_study(
	study_name="test_inductor",
	storage=f"sqlite:///{db_path}",
)
assert len(raw_study.trials) == 3

for i, trial in enumerate(raw_study.trials):
	print(f"\nTrial #{trial.number}:")
	print(f"  params: {trial.params}")
	print(f"  value:  {trial.values[0]:.4e} H")
	print(f"  attrs:  n_free_dofs={trial.user_attrs.get('n_free_dofs')}")
	print(f"          display={trial.user_attrs.get('inductance_display')}")

# Cleanup (SQLite may hold lock on Windows, ignore errors)
import shutil
try:
	shutil.rmtree(tmpdir, ignore_errors=True)
except Exception:
	pass

print("\n=== ALL TESTS PASSED ===")
