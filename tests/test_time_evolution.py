"""Tests for the raw time-evolution figures: the plan, the reader, the plot and the namespace wiring."""
from __future__ import annotations

import math
import re
from types import SimpleNamespace

import numpy as np
import pytest

from vlbipy import VLBIObs
from vlbipy.models import Antenna, FreqSetup, ObsMetadata, Scan, ScanSNRSurvey, Stokes
from vlbipy.selection import (FRINGE_COHERENCE, MIN_COHERENCE_PAIRS, NOISE_COHERENCE, format_runs, lag_coherence,
                              plan_time_evolution, plan_time_evolution_subbands, reference_antenna_gaps)

NAN = float("nan")
#: How the dashboard reads the reference antenna and the subband back from a file name.
DASHBOARD_NAME = re.compile(r"\.time_evolution\.ref([^.]+)\.spw(\d+)\.png$")


def metadata(subbands: dict, absent: dict = None, n_scans: int = 12, n_subbands: int = 4, unobserved=()) -> ObsMetadata:
    """Small observation: ``subbands`` maps antenna -> recorded subbands, ``absent`` antenna -> scans it misses."""
    absent = absent or {}
    antennas = {name: Antenna(name=name, observed=True, subbands=tuple(spws)) for name, spws in subbands.items()}
    antennas.update({name: Antenna(name=name, observed=False) for name in unobserved})
    scans = [Scan(scan_number=k, source="CAL" if k % 2 else "TGT", time_start=5.0e9 + 600.0 * k,
                  time_end=5.0e9 + 600.0 * k + 118.0, integration_time=2.0,
                  antennas=[name for name in subbands if k not in absent.get(name, ())])
             for k in range(1, n_scans + 1)]
    return ObsMetadata(project_code="te", antennas=antennas, scans=scans, source_names=["CAL", "TGT"],
                       time_range=(scans[0].time_start, scans[-1].time_end),
                       freq_setup=FreqSetup(ref_freq=5.0e9, n_subbands=n_subbands, n_channels=16, channel_width=5.0e5,
                                            polarizations=[Stokes.RR, Stokes.RL, Stokes.LR, Stokes.LL]))


def survey(snr: dict, refant: str = "EF") -> ScanSNRSurvey:
    """One-scan, one-hand survey: ``snr`` maps antenna -> fringe SNR (nan = not detected)."""
    return ScanSNRSurvey(project_code="te", scan_numbers=[1], scan_sources=["CAL"], antennas=list(snr),
                         snr={"RR": [list(snr.values())]}, refant=refant)


FULL = (0, 1, 2, 3)


# -- the plan -------------------------------------------------------------------------------------

def test_format_runs_writes_ranges_and_counts_what_it_leaves_out():
    assert format_runs([7, 2, 3, 4, 3]) == "2-4, 7"
    assert format_runs([]) == ""
    assert format_runs([1, 3, 5, 7, 9, 10, 11]) == "1, 3, 5, 7 and 3 more"


def test_one_figure_when_every_antenna_shares_a_subband_and_the_reference_antenna_is_always_there():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL, "CM": (1, 2), "DE": (1, 2)}, unobserved=("TR",))
    plan = plan_time_evolution(meta, None, "EF")
    assert [(entry["refant"], entry["subband"]) for entry in plan] == [("EF", 1)]      # 1 and 2 tie: lower index
    assert plan[0]["reason"] == "reference antenna; subband 1 is recorded by all 5 antennas with data"
    assert plan[0]["covers"] == []


def test_the_subband_nearest_the_band_centre_wins_a_tie():
    meta = metadata({"EF": range(8), "WB": range(8)}, n_subbands=8)
    assert plan_time_evolution_subbands(meta)[0][0] == 3


def test_partial_band_antennas_force_a_second_subband():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": (0, 1), "CM": (2, 3), "DE": (2, 3)})
    subbands = plan_time_evolution_subbands(meta)
    assert [subband for subband, _ in subbands] == [2, 1]                  # 0 and 1 both add JB: nearer the centre
    assert subbands[0][1] == "subband 2 is recorded by 4 of the 5 antennas with data, the most of any subband"
    assert subbands[1][1] == "subband 1 because JB did not record subband 2"
    plan = plan_time_evolution(meta, None, "EF")
    assert [(entry["refant"], entry["subband"]) for entry in plan] == [("EF", 2), ("EF", 1)]


def test_only_subbands_of_the_reference_antenna_are_considered():
    meta = metadata({"EF": (2, 3), "WB": FULL, "JB": FULL, "MC": (0, 1), "CM": (0, 1)})
    assert [subband for subband, _ in plan_time_evolution_subbands(meta)] == [1, 2]
    subbands = plan_time_evolution_subbands(meta, reference="EF")
    assert subbands == [(2, "subband 2 is recorded by 3 of the 5 antennas with data, the most of any subband "
                            "EF recorded")]
    assert plan_time_evolution_subbands(meta, 1, reference="WB", wanted=["MC", "CM"])[0][0] == 1


def test_subbands_stop_at_the_limit_and_when_nothing_is_gained():
    meta = metadata({"EF": (0,), "WB": (0,), "JB": (1,), "MC": (1,), "CM": (2,), "DE": (2,), "KN": (3,), "PI": (3,)})
    assert len(plan_time_evolution_subbands(meta, max_subbands=3)) == 3
    assert len(plan_time_evolution_subbands(metadata({"EF": FULL, "WB": FULL}))) == 1


def test_a_reference_antenna_missing_scans_gets_a_second_one_chosen_by_fringe_snr():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL, "MC": FULL}, absent={"EF": (4, 5, 6, 9)})
    plan = plan_time_evolution(meta, survey({"EF": NAN, "WB": 20.0, "JB": 90.0, "MC": 50.0}), "EF")
    assert [(entry["refant"], entry["subband"]) for entry in plan] == [("EF", 1), ("JB", 1)]
    assert plan[1]["covers"] == [4, 5, 6, 9]
    assert plan[1]["reason"].startswith("EF is absent from scans 4-6, 9; subband 1")


def test_an_undetected_antenna_is_not_made_a_reference_antenna():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL}, absent={"EF": (3,)})
    plan = plan_time_evolution(meta, survey({"EF": NAN, "WB": NAN, "JB": 15.0}), "EF")
    assert [entry["refant"] for entry in plan] == ["EF", "JB"]
    plan = plan_time_evolution(meta, survey({"EF": NAN, "WB": NAN, "JB": NAN}), "EF")
    assert [entry["refant"] for entry in plan] == ["EF"]                    # nobody the survey can vouch for


def test_the_reference_antenna_of_the_survey_is_the_first_candidate_when_another_one_was_asked_for():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL}, absent={"WB": (1, 2)})
    plan = plan_time_evolution(meta, survey({"EF": NAN, "WB": 30.0, "JB": 80.0}, refant="EF"), "WB")
    assert [entry["refant"] for entry in plan] == ["WB", "EF"]


def test_without_a_survey_the_antenna_covering_the_most_missing_scans_is_added():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL, "MC": FULL, "O8": FULL},
                    absent={"EF": (2, 3, 4, 5), "WB": (3, 4, 5), "JB": (4, 5), "O8": (2,)})
    plan = plan_time_evolution(meta, None, "EF")
    assert [entry["refant"] for entry in plan] == ["EF", "MC"]             # MC is in all four, O8 in three
    assert plan[1]["covers"] == [2, 3, 4, 5]


def test_a_scan_with_a_single_antenna_needs_no_reference_antenna():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL}, absent={"EF": (2, 3), "WB": (2,), "JB": (2, 3)})
    assert [entry["refant"] for entry in plan_time_evolution(meta, None, "EF")] == ["EF"]


def test_several_reference_antennas_are_added_until_the_limit():
    # O8 is in every scan (so each has a baseline) but the survey never detected it: it cannot be a reference.
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL, "MC": FULL, "O8": FULL},
                    absent={"EF": (1, 2, 3, 4), "WB": (3, 4, 5), "JB": (1, 2, 6), "MC": (1, 2, 3, 4, 5, 6, 7, 8)})
    ranked = survey({"EF": NAN, "WB": 90.0, "JB": 50.0, "MC": 10.0, "O8": NAN})
    plan = plan_time_evolution(meta, ranked, "EF")
    assert [entry["refant"] for entry in plan] == ["EF", "WB", "JB"]
    assert plan[1]["covers"] == [1, 2] and plan[2]["covers"] == [3, 4]
    assert plan[2]["reason"].startswith("EF is absent from scans 3-4 (not covered by WB either)")
    assert [entry["refant"] for entry in plan_time_evolution(meta, ranked, "EF", max_refants=2)] == ["EF", "WB"]


def test_antennas_sharing_no_subband_with_the_reference_antenna_get_a_figure_on_another_one():
    meta = metadata({"EF": (2, 3), "WB": FULL, "JB": FULL, "MC": (0, 1), "CM": (0, 1)})
    plan = plan_time_evolution(meta, None, "EF")
    assert [(entry["refant"], entry["subband"]) for entry in plan] == [("EF", 2), ("WB", 1)]
    assert plan[1]["reason"] == "MC, CM share no subband with EF; subband 1 on WB" and plan[1]["covers"] == []
    ranked = survey({"EF": NAN, "WB": 20.0, "JB": 70.0, "MC": 30.0, "CM": 5.0})
    assert plan_time_evolution(meta, ranked, "EF")[1]["refant"] == "JB"    # the best detected one that can show them


def test_an_added_reference_antenna_without_the_chosen_subband_gets_one_it_recorded():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": (0, 1), "CM": (2, 3), "DE": (2, 3)}, absent={"EF": (5, 6)})
    ranked = survey({"EF": NAN, "WB": 20.0, "JB": 90.0, "CM": 10.0, "DE": 10.0})
    plan = plan_time_evolution(meta, ranked, "EF", max_subbands=1)
    assert [(entry["refant"], entry["subband"]) for entry in plan] == [("EF", 2), ("JB", 1)]
    assert plan[1]["reason"] == "EF is absent from scans 5-6; subband 1 because JB did not record subband 2"


def test_scans_without_data_or_without_fringes_count_as_not_covered():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL})
    ranked = survey({"EF": NAN, "WB": 40.0, "JB": 60.0})
    assert len(plan_time_evolution(meta, ranked, "EF")) == 1
    plan = plan_time_evolution(meta, ranked, "EF", without_data={"EF": [2, 3]}, without_fringes={"EF": [8]})
    assert [entry["refant"] for entry in plan] == ["EF", "JB"]
    assert plan[1]["covers"] == [2, 3, 8]
    assert plan[1]["reason"].startswith("EF has no unflagged data in scans 2-3, shows no fringes in scan 8;")


def test_a_reference_antenna_chain_is_reduced_to_its_first_member_with_data_and_blank_picks_one():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL}, unobserved=("TR",))
    assert plan_time_evolution(meta, None, "TR,WB,EF")[0]["refant"] == "WB"
    assert plan_time_evolution(meta, None, "")[0]["refant"] == "EF"         # the built-in priority list


def test_nothing_is_planned_without_a_baseline():
    assert plan_time_evolution(metadata({"EF": FULL}), None, "EF") == []
    assert plan_time_evolution(metadata({"EF": (0,), "WB": (1,)}), None, "EF") == []


# -- what the data says about the reference antenna -------------------------------------------------

def fringes(n: int, rate: float = 0.03, amplitude: float = 3.0, seed: int = 0) -> np.ndarray:
    """``n`` samples of a signal turning ``rate`` cycles per sample, in unit complex noise."""
    rng = np.random.default_rng(seed)
    return amplitude * np.exp(2j * np.pi * rate * np.arange(n)) + rng.normal(size=n) + 1j * rng.normal(size=n)


def noise(n: int, seed: int = 0) -> np.ndarray:
    """``n`` samples of unit complex noise."""
    rng = np.random.default_rng(seed)
    return rng.normal(size=n) + 1j * rng.normal(size=n)


def test_lag_coherence_tells_a_turning_fringe_from_noise():
    turning = fringes(60)
    assert lag_coherence(turning) > FRINGE_COHERENCE
    assert lag_coherence(fringes(60, rate=0.2)) > FRINGE_COHERENCE       # the rate does not matter
    assert lag_coherence(noise(60)) < NOISE_COHERENCE
    # A plain vector average over the scan cancels the same detection: that is what the lag avoids.
    assert abs(np.mean(turning)) / np.mean(np.abs(turning)) < NOISE_COHERENCE


def test_lag_coherence_needs_enough_samples_and_skips_missing_ones():
    assert math.isnan(lag_coherence(fringes(MIN_COHERENCE_PAIRS)))
    assert math.isnan(lag_coherence(np.zeros(60)))
    gappy = fringes(80)
    gappy[::4] = np.nan
    assert lag_coherence(gappy) > 0.5


def series_for(meta: ObsMetadata, refant: str, per_scan: dict, others=("WB", "JB")) -> dict:
    """A reader-like series: ``per_scan`` maps scan number -> 'fringes', 'noise' or None (no data)."""
    baselines, with_data = {}, set()
    for index, name in enumerate(others):
        times, scans, vis = [], [], []
        for scan in meta.scans:
            kind = per_scan.get(scan.scan_number, "fringes")
            if kind is None:
                continue
            seed = scan.scan_number + index
            samples = fringes(59, seed=seed) if kind == "fringes" else noise(59, seed=scan.scan_number)
            times.extend(scan.time_start + 2.0 * np.arange(59))
            scans.extend([scan.scan_number] * 59)
            vis.extend(np.stack([samples, samples], axis=1))
            with_data.add(scan.scan_number)
        baselines[name] = {"times": np.asarray(times), "scans": np.asarray(scans, dtype=np.int32),
                           "vis": np.asarray(vis).reshape(-1, 2)}
    return {"refant": refant, "subband": 1, "column": "data", "polarizations": ["RR", "LL"], "antennas": list(others),
            "baselines": baselines, "notes": {}, "scan_sources": {s.scan_number: s.source for s in meta.scans},
            "sources": list(meta.source_names), "scans_with_data": sorted(with_data), "channels": [6, 9],
            "n_channels": 16, "frequency_ghz": 5.0, "bandwidth_mhz": 2.0, "time_range": list(meta.time_range),
            "n_samples": sum(len(b["times"]) for b in baselines.values())}


def test_gaps_are_the_scans_without_data_and_the_scans_without_fringes_on_a_source_that_is_seen():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL}, absent={"EF": (12,)})
    gaps = reference_antenna_gaps(series_for(meta, "EF", {4: None, 5: None, 7: "noise", 12: None}), meta)
    assert gaps == {"no_data": [4, 5], "no_fringes": [7]}          # scan 12: EF is not listed there at all


def test_a_source_nobody_fringes_on_says_nothing_about_the_reference_antenna():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL})
    target_scans = {scan.scan_number: "noise" for scan in meta.scans if scan.source == "TGT"}
    assert reference_antenna_gaps(series_for(meta, "EF", target_scans), meta) == {"no_data": [], "no_fringes": []}


def test_a_reference_antenna_with_no_series_at_all_has_no_data_anywhere():
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL})
    empty = series_for(meta, "EF", {scan.scan_number: None for scan in meta.scans})
    assert reference_antenna_gaps(empty, meta)["no_data"] == list(range(1, 13))


# -- the reader -------------------------------------------------------------------------------------

def test_central_channel_range_matches_the_casa_selection_string():
    from vlbipy.backends.casa import central_channel_range, central_channel_selection
    assert central_channel_range(64, 0.25) == (24, 39)
    assert central_channel_range(32, 0.8) == (3, 28) and central_channel_selection(32, 0.8) == "*:3~28"
    assert central_channel_range(32, 1.0) == (0, 31) and central_channel_range(2, 0.5) == (0, 1)
    assert central_channel_range(64, 0.0) == (31, 32)                    # never an empty window


def test_channels_are_vector_averaged_over_what_is_unflagged_and_recorded():
    from vlbipy.backends.casa import average_recorded_channels
    values = np.ones((4, 4, 3), dtype=complex)                          # (correlation, channel, row)
    values[0, :, 0] = [1, 1j, -1, -1j]                                  # RR of row 0 cancels: a vector average
    values[3, :, 1] = 0.0                                               # LL of row 1 was never recorded
    flags = np.zeros((4, 4, 3), dtype=bool)
    flags[0, :, 2] = True                                               # RR of row 2 is flagged
    flags[3, :2, 2] = True
    values[3, 2:, 2] = 2.0
    averaged = average_recorded_channels(values, flags, [0, 3])
    assert averaged.shape == (3, 2) and averaged.dtype == np.complex64
    assert abs(averaged[0, 0]) < 1e-6 and averaged[0, 1] == 1.0
    assert averaged[1, 0] == 1.0 and np.isnan(averaged[1, 1])
    assert np.isnan(averaged[2, 0]) and averaged[2, 1] == 2.0


class FakeMsTool:
    """Stand-in for ``casatools.ms``: answers ``getdata`` from scripted rows, honouring the selections."""

    def __init__(self, rows: dict) -> None:
        self.rows, self.opened, self.selections, self.closed = rows, None, [], False
        self.channels, self.selection = (0, 0), {}

    def open(self, path, nomodify=True):
        self.opened = (path, nomodify)
        return True

    def selectinit(self, datadescid=0):
        self.spw = datadescid
        return True

    def selectchannel(self, nchan, start, width, inc):
        self.channels = (start, start + nchan)
        return width == 1 and inc == 1

    def msselect(self, selection):
        self.selection = selection
        self.selections.append(dict(selection, spw=self.spw))
        return True

    def getdata(self, items):
        wanted = [int(number) for number in self.selection["scan"].split(",")]
        keep = np.isin(self.rows["scan_number"], wanted)
        if not keep.any():
            return {}
        first, stop = self.channels
        record = {key: self.rows[key][keep] for key in ("antenna1", "antenna2", "time", "scan_number")}
        record["data"] = self.rows["data"][:, first:stop, keep]
        record["flag"] = self.rows["flag"][:, first:stop, keep]
        return record

    def reset(self):
        self.selection = {}

    def close(self):
        self.closed = True


def scripted_rows(meta: ObsMetadata, refant: str) -> dict:
    """Rows of every baseline to ``refant``, two integrations per scan, value = (1 + 0.1j) x scan number."""
    names = list(meta.antennas)
    antenna1, antenna2, times, scans = [], [], [], []
    for scan in meta.scans:
        for other in scan.antennas:
            if other == refant or refant not in scan.antennas:
                continue
            pair = sorted((names.index(refant), names.index(other)))
            for step in range(2):
                antenna1.append(pair[0])
                antenna2.append(pair[1])
                times.append(scan.time_start + 2.0 * step)
                scans.append(scan.scan_number)
    scans = np.asarray(scans)
    data = np.ones((4, 16, scans.size), dtype=complex) * (1 + 0.1j) * scans
    return {"antenna1": np.asarray(antenna1), "antenna2": np.asarray(antenna2), "time": np.asarray(times),
            "scan_number": scans, "data": data, "flag": np.zeros(data.shape, dtype=bool)}


def casa_reader(rows: dict):
    """A ``CasaDataOps`` whose backend hands out a :class:`FakeMsTool` over ``rows`` (no CASA needed)."""
    from vlbipy.backends.casa import CasaDataOps
    tool = FakeMsTool(rows)
    ops = CasaDataOps.__new__(CasaDataOps)
    ops._backend = SimpleNamespace(tools=SimpleNamespace(ms=lambda: tool), ms_path=lambda code: f"/nowhere/{code}.ms")
    return ops, tool


def test_reader_selects_the_central_channels_of_the_baselines_to_the_reference_antenna_read_only():
    meta = metadata({"WB": FULL, "EF": FULL, "JB": FULL, "CM": (2, 3)}, absent={"JB": (1, 2)}, n_scans=25,
                    unobserved=("TR",))
    rows = scripted_rows(meta, "EF")
    names = list(meta.antennas)
    rows["data"][:, :, rows["antenna2"] == names.index("CM")] = 0.0          # CM never recorded subband 1
    rows["flag"][:, :, (rows["scan_number"] == 3) & (rows["antenna2"] == names.index("JB"))] = True
    ops, tool = casa_reader(rows)
    series = ops.read_time_evolution("te", refant="EF", subband=1, metadata=meta)

    assert tool.opened == ("/nowhere/te.ms", True) and tool.closed
    assert tool.channels == (6, 10) and series["channels"] == [6, 9]      # the central quarter of 16 channels
    assert all(selection["baseline"] == "EF&*" and selection["spw"] == 1 for selection in tool.selections)
    assert [len(selection["scan"].split(",")) for selection in tool.selections] == [10, 10, 5]   # bounded reads
    assert series["antennas"] == ["WB", "JB", "CM"]                       # table order; TR never observed
    assert series["polarizations"] == ["RR", "LL"] and series["sources"] == ["CAL", "TGT"]
    assert series["notes"] == {"CM": "subband not recorded"}
    assert len(series["baselines"]["CM"]["times"]) == 0
    assert series["scans_with_data"] == list(range(1, 26))
    jb = series["baselines"]["JB"]
    assert sorted(set(jb["scans"].tolist())) == list(range(4, 26))         # absent from 1-2, flagged in 3
    assert np.all(np.diff(jb["times"]) >= 0) and jb["vis"].shape == (44, 2)
    assert series["n_samples"] == 50 + 44 and series["frequency_ghz"] > 0 and series["bandwidth_mhz"] == 2.0


def test_reader_turns_every_baseline_so_that_the_reference_antenna_leads():
    meta = metadata({"WB": FULL, "EF": FULL, "JB": FULL}, n_scans=2)
    ops, _ = casa_reader(scripted_rows(meta, "EF"))
    series = ops.read_time_evolution("te", refant="EF", subband=0, metadata=meta)
    # Rows are stored lower antenna first: WB-EF is conjugated into EF-WB, EF-JB is kept as it is.
    assert np.allclose(series["baselines"]["WB"]["vis"][0], np.conj(1 + 0.1j))
    assert np.allclose(series["baselines"]["JB"]["vis"][0], 1 + 0.1j)


def test_reader_names_what_it_cannot_read():
    from vlbipy.errors import BackendError
    meta = metadata({"WB": FULL, "EF": FULL}, n_scans=2)
    ops, _ = casa_reader(scripted_rows(meta, "EF"))
    with pytest.raises(BackendError, match="not in the array"):
        ops.read_time_evolution("te", refant="XX", subband=0, metadata=meta)
    with pytest.raises(BackendError, match="subband 9 does not exist"):
        ops.read_time_evolution("te", refant="EF", subband=9, metadata=meta)


def test_reader_reports_no_unflagged_data_for_a_fully_flagged_antenna():
    meta = metadata({"WB": FULL, "EF": FULL, "JB": FULL}, n_scans=2)
    rows = scripted_rows(meta, "EF")
    rows["flag"][:, :, rows["antenna2"] == list(meta.antennas).index("JB")] = True
    series = casa_reader(rows)[0].read_time_evolution("te", refant="EF", subband=0, metadata=meta)
    assert series["notes"] == {"JB": "no unflagged data"} and series["scans_with_data"] == [1, 2]


# -- the figure -------------------------------------------------------------------------------------

def test_file_name_is_the_one_the_dashboard_reads_back():
    from vlbipy.dashboard import time_evolution_key
    from vlbipy.plotting import time_evolution_filename
    name = time_evolution_filename("rsm07", "EF", 1)
    assert name == "rsm07.raw.time_evolution.refEF.spw1.png"
    assert DASHBOARD_NAME.search(name).groups() == ("EF", "1")
    assert time_evolution_key(SimpleNamespace(name=name)) == ("EF", 1)
    assert DASHBOARD_NAME.search(time_evolution_filename("x", "A.1", 12)).groups() == ("A_1", "12")


def test_plot_writes_one_png_and_keeps_a_strip_for_a_baseline_without_data(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt
    from vlbipy.plotting import SOURCE_COLORS, plot_time_evolution
    meta = metadata({"EF": FULL, "WB": FULL, "JB": FULL, "CM": (2, 3)})
    series = series_for(meta, "EF", {}, others=("WB", "JB"))
    series["antennas"] = ["WB", "JB", "CM"]
    series["baselines"]["CM"] = {"times": np.zeros(0), "scans": np.zeros(0, dtype=np.int32),
                                 "vis": np.zeros((0, 2), dtype=np.complex64)}
    series["baselines"]["JB"]["vis"][:, 1] = np.nan                       # a dead polarization
    series["notes"] = {"CM": "subband not recorded"}
    figures = []
    monkeypatch.setattr(plt, "close", figures.append)
    path = plot_time_evolution(series, tmp_path, "te", source_roles={"CAL": "phasecal", "TGT": "target"},
                               reason="reference antenna; subband 1 is recorded by all 4 antennas with data")
    assert path == tmp_path / "te.raw.time_evolution.refEF.spw1.png" and path.stat().st_size > 10_000
    figure = figures[0]
    titles = [axis.get_title(loc="left") for axis in figure.axes]
    assert [t for t in titles if t] == ["EF-WB  RR", "EF-WB  LL", "EF-JB  RR", "EF-JB  LL"]
    texts = [text.get_text() for axis in figure.axes for text in axis.texts]
    assert "EF-CM" in texts and "subband not recorded" in texts and "no unflagged LL data" in texts
    legend = [text.get_text() for text in figure.legends[0].get_texts()]
    assert legend == ["CAL (phase calibrator)", "TGT (target)"]
    assert "reference antenna; subband 1" in figure._suptitle.get_text()
    # One colour per source, the same in every panel; the two hands of a baseline share the amplitude scale.
    drawn = {line.get_color() for axis in figure.axes for line in axis.lines}
    assert drawn == set(SOURCE_COLORS[:2])
    amplitude_axes = [axis for axis in figure.axes if axis.get_title(loc="left").startswith("EF-JB")]
    assert amplitude_axes[0].get_ylim() == amplitude_axes[1].get_ylim()
    plt.close("all")


def test_plot_handles_one_hand_and_a_figure_with_nothing_at_all(tmp_path):
    from vlbipy.plotting import plot_time_evolution
    meta = metadata({"EF": FULL, "WB": FULL})
    series = series_for(meta, "EF", {}, others=("WB",))
    series["polarizations"] = ["RR"]
    series["baselines"]["WB"]["vis"] = series["baselines"]["WB"]["vis"][:, :1]
    assert plot_time_evolution(series, tmp_path, "one").is_file()
    empty = series_for(meta, "EF", {scan.scan_number: None for scan in meta.scans}, others=("WB",))
    empty["baselines"]["WB"]["vis"] = np.zeros((0, 2), dtype=np.complex64)
    assert plot_time_evolution(empty, tmp_path, "none").is_file()        # only the strip, and it says why
    nothing = dict(empty, antennas=[], baselines={})
    assert not plot_time_evolution(nothing, tmp_path, "void").exists()


# -- the namespace ------------------------------------------------------------------------------------

def make_obs(code: str = "te01"):
    """A dummy-backend observation with a fringe finder, a phase calibrator and a target, imported."""
    obs = VLBIObs(code, network="EVN", backend="dummy", target="3C395", phasecal="J1848+3219", fringe_finder="3C345")
    obs.import_data()
    return obs


def test_namespace_makes_one_figure_on_the_reference_antenna_by_default():
    obs = make_obs()
    observation = obs["te01"]
    paths = observation.plot.time_evolution()
    refant = str(observation.refant).split(",")[0]
    assert len(paths) == 1 and "/plots/raw/" in paths[0]
    assert DASHBOARD_NAME.search(paths[0]).groups() == (refant, "3")          # 8 subbands: 3 is next to the centre
    assert paths[0].endswith(f"te01.raw.time_evolution.ref{refant}.spw3.png")


def test_raw_diagnostics_include_the_time_evolution_and_the_calibrated_ones_do_not():
    obs = make_obs("te02")
    assert sum(".time_evolution." in path for path in obs.plot.diagnostics(column="data")) == 1
    assert not any(".time_evolution." in path for path in obs.plot.diagnostics(column="corrected"))


def test_backends_declare_the_operation_and_the_dummy_series_has_the_documented_keys():
    from vlbipy.backends import get_backend
    from vlbipy.backends.base import DataOps
    backend = get_backend("dummy")
    assert backend.supports("plot", "time_evolution") and backend.supports("data", "read_time_evolution")
    meta = backend.data.get_metadata("te03", ["A", "B"], "EVN")
    series = backend.data.read_time_evolution("te03", refant="EF", subband=2, metadata=meta)
    documented = re.findall(r"``(\w+)``", DataOps.read_time_evolution.__doc__.split("Returns")[1])
    assert set(documented) <= set(series), sorted(set(documented) - set(series))
    assert series["scans_with_data"] == [scan.scan_number for scan in meta.scans]
    assert reference_antenna_gaps(series, meta) == {"no_data": [], "no_fringes": []}


def test_namespace_adds_a_reference_antenna_when_the_data_show_a_gap_and_drops_one_that_adds_nothing():
    obs = make_obs("te04")
    observation = obs["te04"]
    data = observation._backend.data
    refant = str(observation.refant).split(",")[0]
    real = type(data).read_time_evolution
    lost = [scan.scan_number for scan in observation.metadata.scans][:2]

    def flagged_at_the_start(project_code, *, refant, subband, **kwargs):
        series = real(data, project_code, refant=refant, subband=subband, **kwargs)
        if refant == primary:
            series["scans_with_data"] = [n for n in series["scans_with_data"] if n not in lost]
        return series

    primary = refant
    data.read_time_evolution = flagged_at_the_start
    paths = observation.plot.time_evolution()
    assert len(paths) == 2 and f".ref{refant}." in paths[0] and f".ref{refant}." not in paths[1]

    def flagged_for_everyone(project_code, *, refant, subband, **kwargs):
        series = real(data, project_code, refant=refant, subband=subband, **kwargs)
        series["scans_with_data"] = [n for n in series["scans_with_data"] if n not in lost]
        return series

    data.read_time_evolution = flagged_for_everyone
    paths = observation.plot.time_evolution()
    assert len(paths) == 1 and f".ref{refant}." in paths[0]                # nobody has those scans: no extra figure


def test_resumed_run_only_adds_the_figures_when_there_is_data_and_none_exists():
    obs = make_obs("te05")
    assert obs["te05"].plot._lacks_time_evolution() is False                # the dummy backend has no data on disk
