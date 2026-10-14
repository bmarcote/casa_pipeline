"""Tests for Source / SourceSet role handling and phase-referencing."""
import pytest

from vlbipy.errors import SourceNotFoundError
from vlbipy.models import SourceType
from vlbipy.sources import SourceSet


def _cfg(**roles):
    return {k: v for k, v in roles.items()}


def test_roles_and_accessors():
    ss = SourceSet.from_config(_cfg(targets="3C286", phase_calibrators=["J1048+7143"],
                                    fringe_finders=["3C345"], check_sources=["JCHK"]))
    assert ss.target.name == "3C286"
    assert [s.name for s in ss.phase_calibrators] == ["J1048+7143"]
    assert [s.name for s in ss.fringe_finders] == ["3C345"]
    assert [s.name for s in ss.check_sources] == ["JCHK"]
    assert {s.name for s in ss.calibrators} == {"J1048+7143", "3C345"}
    assert ss["3C345"].source_type is SourceType.FRINGE_FINDER
    assert "3C286" in ss and len(ss) == 4


def test_ambiguous_target_raises():
    ss = SourceSet.from_config(_cfg(targets=["A", "B"]))
    with pytest.raises(SourceNotFoundError):
        _ = ss.target
    assert [s.name for s in ss.targets] == ["A", "B"]


def test_missing_target_raises():
    ss = SourceSet.from_config(_cfg(fringe_finders=["3C345"]))
    with pytest.raises(SourceNotFoundError):
        _ = ss.target


def test_unknown_name_raises():
    ss = SourceSet.from_config(_cfg(targets="3C286"))
    with pytest.raises(SourceNotFoundError):
        _ = ss["nope"]


def test_phase_referencing_explicit_mapping():
    ss = SourceSet.from_config(
        _cfg(targets=["T1"], phase_calibrators=["C1", "C2"], fringe_finders=["FF"]),
        phaseref_cfg={"T1": ["C2", "C1"]})
    cals = ss.calibrators_for("T1")
    assert [c.name for c in cals] == ["C2", "C1"]


def test_phase_referencing_fallback_order():
    # No mapping, but fringe finder present -> fall back to fringe finders.
    ss = SourceSet.from_config(_cfg(targets=["T1"], phase_calibrators=["C1"], fringe_finders=["FF"]))
    assert [c.name for c in ss.calibrators_for("T1")] == ["FF"]
    # No fringe finder -> fall back to phase calibrators.
    ss2 = SourceSet.from_config(_cfg(targets=["T1"], phase_calibrators=["C1"]))
    assert [c.name for c in ss2.calibrators_for("T1")] == ["C1"]
    # Nothing -> fall back to the target itself.
    ss3 = SourceSet.from_config(_cfg(targets=["T1"]))
    assert [c.name for c in ss3.calibrators_for("T1")] == ["T1"]


# -- role inference from the data ---------------------------------------------------------------

def _phaseref_metadata(cycle_scans=30, finder_scans=2, cal_duration=60.0, target_duration=180.0):
    """Metadata of a classic phase-referenced run: a fringe-finder snapshot then cal/target cycles."""
    from vlbipy.models import ObsMetadata, Scan

    scans, time, number = [], 0.0, 1
    for _ in range(finder_scans):
        scans.append(Scan(scan_number=number, source="FF", time_start=time, time_end=time + 300.0,
                          antennas=["EF", "WB", "JB"]))
        time, number = time + 400.0, number + 1
    for _ in range(cycle_scans):
        for source, duration in (("CAL", cal_duration), ("TGT", target_duration)):
            scans.append(Scan(scan_number=number, source=source, time_start=time,
                              time_end=time + duration, antennas=["EF", "WB", "JB"]))
            time, number = time + duration + 20.0, number + 1
    return ObsMetadata(project_code="p", scans=scans, source_names=["FF", "CAL", "TGT"],
                       source_coords={"FF": (180.0, 40.0), "CAL": (200.0, 41.0), "TGT": (200.5, 41.2)})


def _survey_for(metadata, levels: dict[str, float]):
    """One-polarization survey giving every antenna of every scan its source's SNR level."""
    from vlbipy.models import ScanSNRSurvey

    antennas = ["EF", "WB", "JB"]
    rows = [[float("nan")] + [levels.get(s.source, 0.0)] * (len(antennas) - 1) for s in metadata.scans]
    return ScanSNRSurvey(project_code="p", scan_numbers=[s.scan_number for s in metadata.scans],
                         scan_sources=[s.source for s in metadata.scans], antennas=antennas,
                         snr={"RR": rows}, refant="EF")


def test_inference_finds_the_fringe_finder_phase_calibrator_and_target():
    from vlbipy.sources import infer_source_roles

    metadata = _phaseref_metadata()
    survey = _survey_for(metadata, {"FF": 900.0, "CAL": 50.0, "TGT": 30.0})
    roles, phaseref = infer_source_roles(metadata, survey)
    assert roles == {"FF": SourceType.FRINGE_FINDER, "CAL": SourceType.PHASE_CALIBRATOR,
                     "TGT": SourceType.TARGET}
    assert phaseref == {"TGT": ["CAL"]}


def test_a_target_brighter_than_its_calibrator_is_still_the_target():
    # Brightness cannot decide a pair: the scan pattern (who opens the cycle, who is visited
    # more often, whose scans are shorter) does.
    from vlbipy.sources import infer_source_roles

    metadata = _phaseref_metadata()
    roles, _ = infer_source_roles(metadata, _survey_for(metadata, {"FF": 900.0, "CAL": 20.0, "TGT": 400.0}))
    assert roles["CAL"] is SourceType.PHASE_CALIBRATOR and roles["TGT"] is SourceType.TARGET


def test_a_declared_target_settles_which_member_of_the_pair_is_the_calibrator():
    from vlbipy.sources import infer_source_roles

    metadata = _phaseref_metadata(cal_duration=180.0, target_duration=60.0)   # pattern now misleads
    declared = {"CAL": SourceType.TARGET}
    roles, phaseref = infer_source_roles(metadata, _survey_for(metadata, {"FF": 900.0}), declared=declared)
    assert "CAL" not in roles                                   # declared roles are never returned
    assert roles["TGT"] is SourceType.PHASE_CALIBRATOR
    assert phaseref == {"CAL": ["TGT"]}


def test_sources_far_apart_are_not_treated_as_a_phase_reference_pair():
    from vlbipy.sources import infer_source_roles

    metadata = _phaseref_metadata()
    metadata.source_coords["TGT"] = (20.0, -40.0)               # ~100 deg from CAL
    roles, phaseref = infer_source_roles(metadata, _survey_for(metadata, {"FF": 900.0, "CAL": 50.0}))
    assert roles["CAL"] is SourceType.TARGET and roles["TGT"] is SourceType.TARGET
    assert phaseref == {}


def test_inference_works_without_a_survey():
    from vlbipy.sources import infer_source_roles

    roles, phaseref = infer_source_roles(_phaseref_metadata(), None)
    assert roles == {"FF": SourceType.FRINGE_FINDER, "CAL": SourceType.PHASE_CALIBRATOR,
                     "TGT": SourceType.TARGET}
    assert phaseref == {"TGT": ["CAL"]}


def test_an_undetected_snapshot_source_is_not_made_a_fringe_finder():
    from vlbipy.sources import infer_source_roles

    metadata = _phaseref_metadata()
    roles, _ = infer_source_roles(metadata, _survey_for(metadata, {"CAL": 50.0, "TGT": 30.0}))
    assert roles["FF"] is SourceType.OTHER


def test_a_snapshot_run_makes_its_best_source_the_fringe_finder_and_the_rest_targets():
    from vlbipy.models import ObsMetadata, Scan
    from vlbipy.sources import infer_source_roles

    scans = [Scan(scan_number=i + 1, source=name, time_start=i * 400.0, time_end=i * 400.0 + 300.0,
                  antennas=["EF", "WB", "JB"]) for i, name in enumerate(["A", "B", "C"])]
    metadata = ObsMetadata(project_code="p", scans=scans, source_names=["A", "B", "C"])
    roles, phaseref = infer_source_roles(metadata, _survey_for(metadata, {"A": 10.0, "B": 900.0, "C": 20.0}))
    assert roles == {"B": SourceType.FRINGE_FINDER, "A": SourceType.TARGET, "C": SourceType.TARGET}
    assert phaseref == {}


def test_a_single_source_observation_is_its_own_target():
    from vlbipy.models import ObsMetadata, Scan
    from vlbipy.sources import infer_source_roles

    scans = [Scan(scan_number=i + 1, source="ONLY", time_start=i * 400.0, time_end=i * 400.0 + 300.0,
                  antennas=["EF", "WB"]) for i in range(8)]
    metadata = ObsMetadata(project_code="p", scans=scans, source_names=["ONLY"])
    roles, _ = infer_source_roles(metadata, None)
    assert roles == {"ONLY": SourceType.TARGET}


def test_adding_an_inferred_source_never_overrides_a_declared_one():
    ss = SourceSet.from_config(_cfg(targets="TGT"))
    from vlbipy.sources import Source

    assert ss.add(Source(name="CAL", source_type=SourceType.PHASE_CALIBRATOR)) is True
    assert ss.add(Source(name="TGT", source_type=SourceType.PHASE_CALIBRATOR)) is False
    assert ss["TGT"].source_type is SourceType.TARGET
    assert [s.name for s in ss.phase_calibrators] == ["CAL"]


def test_pipeline_infers_roles_for_undeclared_sources_and_fills_the_calibrator_field():
    from vlbipy import VLBIObs

    obs = VLBIObs(project="RSM07", network="EVN", backend="dummy")
    obs.import_data()
    observation = obs["RSM07"]
    assert observation.sources.names == [] and observation.calibrator_field == ""
    roles = observation.calibrate.infer_source_roles()
    assert roles and observation.calibrator_field == ",".join(observation.metadata.source_names)
    assert observation.calibrate.infer_source_roles() == {}      # nothing left to infer


def test_declared_roles_stop_the_inference_from_running():
    from vlbipy import VLBIObs

    obs = VLBIObs(project="RSM07", network="EVN", backend="dummy", target="UNKNOWN")
    obs.import_data()
    observation = obs["RSM07"]
    assert observation.calibrate.infer_source_roles() == {}
    assert observation.sources["UNKNOWN"].source_type is SourceType.TARGET


def test_inferred_roles_are_persisted_and_reused_by_a_later_run(tmp_path, monkeypatch):
    """The inference costs a fringe fit per scan, so a resumed run must not repeat it."""
    from vlbipy import VLBIObs

    work = tmp_path / "proj"
    obs = VLBIObs(project="RSM07", network="EVN", backend="casa", work_dir=str(work))
    observation = obs["RSM07"]
    observation._metadata = _phaseref_metadata()
    observation._metadata.project_code = "RSM07"
    surveys = []

    def fake_scan_snr(*args, **kwargs):
        surveys.append(kwargs.get("field", ""))
        return _survey_for(observation._metadata, {"FF": 900.0, "CAL": 50.0, "TGT": 30.0})

    monkeypatch.setattr(observation.calibrate, "scan_snr", fake_scan_snr)
    roles = observation.calibrate.infer_source_roles()
    assert {n: r.value for n, r in roles.items()} == {"FF": "fringefinder", "CAL": "phasecal",
                                                      "TGT": "target"}
    assert len(surveys) == 1 and (work / ".source_roles.json").is_file()

    # A fresh Observation over the same work_dir reuses the stored answer: no second survey.
    again = VLBIObs(project="RSM07", network="EVN", backend="casa", work_dir=str(work))["RSM07"]
    again._metadata = observation._metadata
    monkeypatch.setattr(again.calibrate, "scan_snr", fake_scan_snr)
    reused = again.calibrate.infer_source_roles()
    assert {n: r.value for n, r in reused.items()} == {n: r.value for n, r in roles.items()}
    assert len(surveys) == 1
    assert again.sources.phase_referencing == {"TGT": ["CAL"]}
