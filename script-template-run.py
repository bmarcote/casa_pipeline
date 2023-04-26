import glob
import logging
import casa_pipeline as capi
# from importlib import reload

params = {'projectcode': 'rg013a',
          'obsdate': '221016',
          'password': '80f5fab7fcff',
          # 'target': ['GRB221009A'],
          # 'phaseref': ['J1908+2028', 'J1905+1943'],
          'phaseref': ['GRB221009A'],
          'target': ['J1908+2028', 'J1905+1943'],
          'fringefinder': ['J1800+3848', 'J1925+2106'],
          'reference_antenna': ['T6', 'MC'],  # can be a list of antennas too
          # 'calibration': {'sbd_timerange': '10:01:30~10:03:30'}}
          'calibration': {'sbd_timerange': '2022/10/16/12:00:00~2022/10/16/12:03:00'}}


obs = capi.Project(params['projectcode'], params=params, logging_level=logging.DEBUG)

# Importing FITS-IDI
if not obs.ms.msfile.exists():
    if len(glob.glob(f"{params['projectcode']}_1_1.IDI*")) == 0:
        obs.importdata.evn_download(obsdate=params['obsdate'], username=params['projectcode'].lower(),
                                    password=params['password'])

    obs.importdata.evn_fitsidi(f"{params['projectcode']}_1_1.IDI*", delete=False)
else:
    # obs.calibrate.clearcal()
    obs.calibrate.callib.remove_entry('all')
    # obs.flag.flagdata(mode='unflag')
    pass

# Manual config (should be in the params auto)
# obs.ms.refant = params['reference_antenna']
obs.store()


# Params to check
# obs.ms.sources
# obs.ms.sources.targets
# obs.ms.sources.fringe_finders
# obs.ms.antennas
# obs.ms.antennas.observed
# obs.ms.timerange
# obs.ms.freqsetup
obs.summary(outfile=f"results/summary-{obs.projectname}.txt")


# Initial flagging
# obs.flag.apply_flags_from_file()
# obs.flag.edge_channels(0.05)  # 0.1 (10% edges) by default
# obs.flag.quack("WB", 4.0)  # 4 s on WB data
# obs.flag.quack("SR", 4.0)  # 4 s on WB data
# obs.flag.aoflagger()

# Initial calibration
obs.calibrate.a_priori_calibration()

# ACCOR, ACSCL.
# obs.calibrate.accor()

# obs.calibrate.ionospheric_corrections()
# obs.calibrate.main_calibration(fixed_fringe_callib=False, dispersive=True)
obs.calibrate.main_calibration(fixed_fringe_callib=False, dispersive=False)
obs.store()
obs.calibrate.apply_calibration()

obs.flag.tfcrop()

splits = []
splits.append(obs.ms.split(sources=obs.sources.phase_calibrators.names + obs.sources.targets.names))
# by default sources=obs.sources.names

for src in obs.sources.phase_calibrators.names + obs.sources.targets.names:
    splits[-1][src].export_uvfits()

# obs.calibrate.recalibration(dispersive=True)
obs.calibrate.recalibration(dispersive=False)
obs.store()
obs.calibrate.apply_calibration()

# obs.flag.tfcrop()

splits.append(obs.ms.split(sources=obs.sources.phase_calibrators.names + obs.sources.targets.names))

for src in obs.sources.phase_calibrators.names + obs.sources.targets.names:
    splits[-1][src].export_uvfits()

