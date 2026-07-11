import copy
import math

from ltspice_converter.bandwidth_gate import measure_bandwidth_crossing_gate
from ltspice_converter.mcp_server import measure_bandwidth_crossing_gate as mcp_gate


def good():
    frequencies = [10.0 + index * 10.0 for index in range(100)]
    center = 500.0
    width = 200.0
    magnitudes = [10.0 / math.sqrt(1.0 + ((frequency - center) / width) ** 4) for frequency in frequencies]
    response_db = [20.0 * math.log10(value) for value in magnitudes]
    peak_db = max(response_db)
    threshold = peak_db - 20.0 * math.log10(math.sqrt(2.0))

    def crossing(rising):
        inside = [i for i, value in enumerate(response_db) if value >= threshold]
        left, right = (inside[0] - 1, inside[0]) if rising else (inside[-1], inside[-1] + 1)
        fraction = (threshold - response_db[left]) / (response_db[right] - response_db[left])
        return frequencies[left] + fraction * (frequencies[right] - frequencies[left])

    lower, upper = crossing(True), crossing(False)
    return {
        "frequency_hz": frequencies,
        "magnitude": magnitudes,
        "measured_peak_db": peak_db,
        "measured_lower_3db_hz": lower,
        "measured_upper_3db_hz": upper,
        "measured_bandwidth_hz": upper - lower,
    }


def test_accepts_two_sided_bandwidth_replay():
    result = measure_bandwidth_crossing_gate(good())
    assert result["status"] == "ok"
    assert result["checks"]["peak_lies_between_crossings"] is True
    assert mcp_gate(good())["status"] == "ok"


def test_rejects_stale_measure_crossings_and_bandwidth():
    payload = copy.deepcopy(good())
    payload["measured_upper_3db_hz"] *= 1.2
    payload["measured_bandwidth_hz"] *= 0.7
    result = measure_bandwidth_crossing_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["upper_crossing_matches_samples"] is False
    assert result["checks"]["bandwidth_matches_samples"] is False


def test_rejects_nonmonotone_frequency_grid():
    payload = good()
    payload["frequency_hz"][20], payload["frequency_hz"][21] = (
        payload["frequency_hz"][21], payload["frequency_hz"][20]
    )
    result = measure_bandwidth_crossing_gate(payload)
    assert result["status"] == "needs_attention"
    assert result["checks"]["frequency_strictly_increases"] is False
