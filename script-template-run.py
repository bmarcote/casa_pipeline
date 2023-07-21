import sys
import glob
import logging
from rich import print as rprint
import casa_pipeline as capi
# from importlib import reload

params = {'projectname': 'ek051c',
          'obsdate': '220817',
          'sci_package': 'AIPS',
          'password': 'RQxHE91GWrMI',
          'sources': {'target': ['R190520_D', 'J1603-1007'],
                      'phaseref': ['J1605-1139'],
                      'fringefinder': ['J1642-0621', 'J1927+6117']},
          'reference_antenna': ['EF', 'TR'],  # can be a list of antennas too
          'calibration': {'sbd_timerange': '2022/08/17/19:31:30~2022/08/17/19:33:00',
                          'do_ionospheric': True}}


obs = capi.Project(params['projectname'], params=params, logging_level=logging.DEBUG,
                   sci_package=params['sci_package'])

# Importing FITS-IDI
if len(glob.glob(f"{obs.projectname}_1_1.IDI*")) == 0:
    obs.importdata.evn_download(obsdate=params['obsdate'], username=params['projectname'].lower(),
                                password=params['password'])

if not obs.msfile.exists():
    if obs.sci_package == 'AIPS':
        uvfits = obs.calibrate.aips.a_priori_calibration()
        obs.importdata.import_uvfits(uvfitsfile=uvfits, delete=True)
    else:
        obs.importdata.evn_fitsidi(f"{obs.projectname}_1_1.IDI*", delete=False)
else:
    obs.flag.flagdata(mode='unflag')
    if not obs.sci_package == 'AIPS':
        obs.calibrate.clearcal()
        obs.calibrate.callib.remove_entry('all')
        obs.flag.apply_flags_from_file()
        obs.flag.flagdata(autocorr=True)


obs.store()

# params to check
# obs.ms.sources
# obs.ms.sources.targets
# obs.ms.sources.fringe_finders
# obs.ms.antennas
# obs.ms.antennas.observed
# obs.ms.timerange
# obs.ms.freqsetup
obs.summary(outfile=f"summary-{obs.projectname}.txt")

# obs.plot.tplot()


# initial flagging
obs.flag.edge_channels(0.05)  # 0.1 (10% edges) by default
# obs.flag.quack("wb", 10.0)  # 10 s on wb data
obs.flag.aoflagger()
obs.flag.tfcrop()


splits = []
if obs.sci_package == 'AIPS':
    obs.export_uvfits(outfitsfilename=f"{obs.projectname}.2.uvfits", overwrite=True)
    splits.append(obs.calibrate.aips.main_calibration(uvfits=f"{obs.projectname}.2.uvfits", avgchan=False))

    for a_src in splits[-1]:
        splits[-1][a_src].flag.aoflagger()
        splits[-1][a_src].flag.tfcrop(ntime='30min')
        # todo: something here that euqals amplitudes accross different subbands. accor?
        splits[-1][a_src].split(inplace=True, datacolumn='data', chanbin=-1)
        splits[-1][a_src].export_uvfits(outfitsfilename=f"{splits[-1][a_src].projectname}.split.2.uvfits",
                                        overwrite=True)

else:
    # initial calibration
    obs.calibrate.a_priori_calibration()
    # accor, acscl.
    # obs.calibrate.accor()
    # obs.calibrate.ionospheric_corrections()
    obs.calibrate.main_calibration(fixed_fringe_callib=False, dispersive=True, parang=True)
    obs.calibrate.apply_calibration()
    obs.flag.tfcrop()
    obs.calibrate.recalibration(dispersive=False, parang=False)
    obs.calibrate.apply_calibration()
    splits.append(obs.split(sources=obs.sources.phase_calibrators.names + obs.sources.targets.names))
    # by default sources=obs.sources.names

    for src in obs.sources.phase_calibrators.names + obs.sources.targets.names:
        splits[-1][src].export_uvfits()

    obs.store()


rprint("\n\n[bold green]Stopping here so you can have fun with difmap[/bold green]")
sys.exit(0)


###### SELF-CALIBRATION AFTER DIFMAP
# obs.calibrate.aips.selfcal_with_difmapimage('ek051c.J1605-1139.difmap.2.fits')
obs.calibrate.aips.selfcal_with_difmapimage('EM161AB.J1605-1139.difmap_noTr.fits')


# Apply to the burst data
obs.calibrate.aips.transfer_calibration(f"{obs.projectname}_2_1.IDI", source="R190520_D")








