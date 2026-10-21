"""Do not accept NGSolve commonGPU host fallback as CUDA execution."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import pytest

spec = importlib.util.spec_from_file_location("gpu_adapter", Path(__file__).parents[1] / "src/radia_simulink/ngsolve_gpu.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_host_fallback_refused_before_device_request(monkeypatch):
    def forbidden():
        raise AssertionError("host device must not be accepted")
    monkeypatch.setattr("importlib.import_module", lambda name: SimpleNamespace(backend="host", GetGPUDevice=forbidden))
    with pytest.raises(RuntimeError, match="actual backend is 'host'"):
        module.cuda_device()


def test_cuda_fp64_required(monkeypatch):
    device = SimpleNamespace(has_float64=False)
    monkeypatch.setattr("importlib.import_module", lambda name: SimpleNamespace(backend="cuda", GetGPUDevice=lambda: device))
    with pytest.raises(RuntimeError, match="FP64"):
        module.cuda_device()
    device.has_float64 = True
    assert module.cuda_device() is device
