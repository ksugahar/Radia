"""Guards: the compact writer's paired q1/q2/q3 faces share a translated map.

Assumes: the native synthetic fixture uses dyadic coefficients and reference
points, making these exact arithmetic checks independent of CAD tolerances.
Full Netgen reads the file; no Curve, Save, Identify or coefficient repair.
"""
from pathlib import Path
from fractions import Fraction
import ctypes
from ctypes import wintypes
import faulthandler
import hashlib
import itertools
import json
import os
import sys
import time

PHASE_FILE = None


def phase(name, **fields):
    row = {"phase": name, "monotonic_seconds": time.monotonic(), **fields}
    text = json.dumps(row, ensure_ascii=True)
    print(text, flush=True)
    if PHASE_FILE is not None:
        with PHASE_FILE.open("a", encoding="utf-8") as stream:
            stream.write(text + "\n")
            stream.flush()


def module_identity():
    paths = {Path(m.__file__).resolve() for name, m in list(sys.modules.items())
             if name.split(".")[0] in {"netgen", "ngsolve"} and getattr(m, "__file__", None)}
    if os.name == "nt":
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel.GetCurrentProcess.restype = wintypes.HANDLE
        enum = psapi.EnumProcessModulesEx
        enum.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD,
                         ctypes.POINTER(wintypes.DWORD), wintypes.DWORD]
        enum.restype = wintypes.BOOL
        filename = psapi.GetModuleFileNameExW
        filename.argtypes = [wintypes.HANDLE, wintypes.HMODULE, wintypes.LPWSTR, wintypes.DWORD]
        filename.restype = wintypes.DWORD
        modules = (wintypes.HMODULE * 4096)()
        needed = wintypes.DWORD()
        process = kernel.GetCurrentProcess()
        if not enum(process, modules, ctypes.sizeof(modules), ctypes.byref(needed), 3):
            raise ctypes.WinError(ctypes.get_last_error())
        if needed.value > ctypes.sizeof(modules):
            raise RuntimeError("Loaded-module inventory exceeded its declared capacity")
        for module in modules[:needed.value // ctypes.sizeof(wintypes.HMODULE)]:
            buffer = ctypes.create_unicode_buffer(32768)
            size = filename(process, module, buffer, len(buffer))
            if not size or size >= len(buffer):
                raise RuntimeError("Cannot resolve loaded native module path")
            path = Path(buffer.value).resolve()
            if path.name.lower() in {"ngcore.dll", "nglib.dll", "ngsolve.dll"} or any(
                    part.lower() in {"netgen", "ngsolve", "netgen_mesher.libs"} for part in path.parts):
                paths.add(path)
        if not {"ngcore.dll", "nglib.dll"}.issubset({p.name.lower() for p in paths}):
            raise RuntimeError("Required loaded Netgen DLL identities are absent")
    return [{"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in sorted(paths)]


def verify_raw_exterior(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    def rows(name, convert):
        i = lines.index(name)
        return [list(map(convert, row.split())) for row in lines[i+2:i+2+int(lines[i+1])]]
    tets = rows("volumeelements", int)
    surfaces = rows("surfaceelements", int)
    descriptors = rows("facedescriptors", int)
    points = rows("points", Fraction)
    incidence = {}
    opposite = {}
    for tet in tets:
        if tet[:2] != [1, 4]:
            raise AssertionError(f"Expected domain1 TET4, actual={tet}")
        for face in itertools.combinations(tet[2:], 3):
            key = tuple(sorted(face)); incidence[key] = incidence.get(key, 0) + 1
            opposite[key] = next(p for p in tet[2:] if p not in face)
    expected = {key for key, count in incidence.items() if count == 1}
    actual = {tuple(sorted(face[-3:])) for face in surfaces}
    if len(surfaces) != 8 or len(actual) != 8 or actual != expected:
        raise AssertionError(f"Exterior faces actual={sorted(actual)}, expected={sorted(expected)}")
    if any(fd[1:3] != [1, 0] for fd in descriptors) or any(s[2:5] != [1, 0, 3] for s in surfaces):
        raise AssertionError("Expected every exterior descriptor/triangle to own domain1->0")
    for face in surfaces:
        key = tuple(sorted(face[-3:]))
        a, b, c = [points[p-1] for p in face[-3:]]
        inside = points[opposite[key]-1]
        u = [b[k]-a[k] for k in range(3)]; v = [c[k]-a[k] for k in range(3)]
        normal = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0])
        if sum(normal[k]*(inside[k]-a[k]) for k in range(3)) >= 0:
            raise AssertionError(f"Exterior face is not outward: {face[-3:]}")


# Independent evaluation of the declared fixture's polynomial seed (canonical
# master edges 12,13,14,23,24,34; only edge14 is inactive; face123 and face234).
# These dyadic coordinates include nonzero even/odd modes and face bubbles;
# translation equality alone would also accept incorrectly flattened output.
MASTER_COORDINATES = {
    1: (("0", "1/4", "1/2"), ("1/2", "3/4", "3/4")),
    2: (("-13/8192", "1011/4096", "4057/8192"),
        ("2035/4096", "1523/2048", "3033/4096")),
    3: (("-13/16384", "2035/8192", "8153/16384"),
        ("4083/8192", "3059/4096", "6105/8192")),
}


def verify(root):
    global PHASE_FILE
    PHASE_FILE = None
    phase_path = root / "periodic-curved-reader-phases.jsonl"
    with phase_path.open("x", encoding="utf-8"):
        pass
    PHASE_FILE = phase_path
    faulthandler.enable(all_threads=True)
    phase("before-import", executable=sys.executable, pid=os.getpid())
    import ngsolve as ng
    from netgen.meshing import IdentificationType

    phase("after-import", version=str(ng.__version__))
    phase("loaded-module-identity", modules=module_identity())
    ng.SetNumThreads(1)
    xyz = ng.CF((ng.x, ng.y, ng.z))
    # The first reference point is on vertices 1/2/3; the second is on
    # vertices 2/3/4. Both exercise their shared edge and q3 face bubbles.
    samples = [(0.25, 0.25, 0.5), (0.0, 0.25, 0.25)]
    for suffix, mapping in (("", (8, 7, 6, 5)), ("-mixed", (5, 7, 6, 8))):
      for order in (1, 2, 3):
        path = root / f"periodic-curved-q{order}{suffix}.vol"
        case = f"q{order}{suffix}"
        phase("before-load", case=case, sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        verify_raw_exterior(path)
        phase("raw-exterior-verified", case=case)
        mesh = ng.Mesh(str(path))
        phase("after-load", case=case)
        actual = {(int(a.nr), int(b.nr), int(nr))
                  for a, b, nr in mesh.ngmesh.GetIdentifications()}
        expected = {(i + 1, mapping[i], 1) for i in range(4)}
        if actual != expected or mesh.ngmesh.GetNrIdentifications() != 1:
            raise AssertionError(f"actual IDs={actual}, expected={expected}")
        if mesh.ngmesh.GetIdentificationType(1) != IdentificationType.PERIODIC:
            raise AssertionError("Expected PERIODIC")
        elements = list(mesh.Elements(ng.VOL))
        if len(elements) != 2:
            raise AssertionError(f"actual tets={len(elements)}, expected=2")
        for index, point in enumerate(samples):
            phase("before-reference-point", case=case, sample=index, reference=point)
            ip = next(iter(ng.IntegrationRule(points=[point], weights=[1.0])))
            phase("after-reference-point", case=case, sample=index)
            phase("before-map", case=case, sample=index, side="master", element=int(elements[0].nr))
            a = tuple(xyz(mesh.GetTrafo(elements[0])(ip)))
            phase("after-map", case=case, sample=index, side="master", xyz=a)
            phase("before-map", case=case, sample=index, side="slave", element=int(elements[1].nr))
            b = tuple(xyz(mesh.GetTrafo(elements[1])(ip)))
            phase("after-map", case=case, sample=index, side="slave", xyz=b)
            expected_a = tuple(float(Fraction(v)) for v in MASTER_COORDINATES[order][index])
            expected_b = (expected_a[0] + 3.0, expected_a[1], expected_a[2])
            if a != expected_a or b != expected_b:
                raise AssertionError(f"q{order}{suffix} point={point}: actual={(a,b)}, expected={(expected_a,expected_b)}")


if __name__ == "__main__":
    try:
        verify(Path(sys.argv[1]))
        phase("complete", modules=module_identity())
    except BaseException as error:
        phase("python-exception", error=repr(error))
        raise
