"""Phase F-1 smoke test: import + version + precision info + Newton kernel."""
import math


def test_import_and_version():
    import radia_vim
    assert radia_vim.version() == "0.1.0"
    assert "Newton-kernel" in radia_vim.hello()


def test_precision_info():
    from radia_vim import precision_info
    info = precision_info()
    # double is always available
    assert info.with_double is True
    # quad and mpfr are conditional
    print(repr(info))


def test_newton_kernel_double_unit_distance():
    from radia_vim import newton_kernel_double
    val = newton_kernel_double([0, 0, 0], [1, 0, 0])
    assert abs(val - 1.0) < 1e-15


def test_newton_kernel_double_known_value():
    from radia_vim import newton_kernel_double
    # 1 / sqrt(3) for r = (1,1,1)
    val = newton_kernel_double([0, 0, 0], [1, 1, 1])
    assert abs(val - 1.0 / math.sqrt(3.0)) < 1e-15


def test_newton_kernel_quad_when_available():
    try:
        from radia_vim import newton_kernel_quad
    except ImportError:
        return  # quad not compiled in
    s = newton_kernel_quad([0, 0, 0], [1, 0, 0])
    assert isinstance(s, str)
    assert s.startswith("1") and len(s) > 20


def test_newton_kernel_mpfr_when_available():
    try:
        from radia_vim import newton_kernel_mpfr
    except ImportError:
        return  # mpfr not compiled in
    s = newton_kernel_mpfr([0, 0, 0], [1, 0, 0], digits=60)
    assert isinstance(s, str)
    assert s.startswith("1") and len(s) >= 50


if __name__ == "__main__":
    test_import_and_version()
    test_precision_info()
    test_newton_kernel_double_unit_distance()
    test_newton_kernel_double_known_value()
    test_newton_kernel_quad_when_available()
    test_newton_kernel_mpfr_when_available()
    print("All smoke tests passed.")
