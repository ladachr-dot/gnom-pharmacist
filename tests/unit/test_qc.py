# tests/unit/test_qc.py
from gpharm_core.qc import _classify_status, QCStatus

def test_classify_pass():
    assert _classify_status(95.0) == QCStatus.PASS

def test_classify_warn():
    assert _classify_status(85.0) == QCStatus.WARN
    assert _classify_status(70.0) == QCStatus.WARN

def test_classify_fail():
    assert _classify_status(69.9) == QCStatus.FAIL
    assert _classify_status(0.0)  == QCStatus.FAIL

def test_classify_boundary():
    assert _classify_status(90.0) == QCStatus.PASS
    assert _classify_status(89.9) == QCStatus.WARN