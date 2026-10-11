"""Small installation utilities; numerical execution stays in mexcuda."""
import argparse
import json
from . import matlab_path, verify


def path_main():
    print(matlab_path())


def doctor_main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    options = parser.parse_args()
    report = verify()
    print(json.dumps(report, indent=2) if options.json else f"radia-gpu {report['version']}: payload verified")
