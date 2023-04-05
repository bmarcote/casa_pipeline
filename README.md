

# CAPI:  CASA Pipeline for radio data

**_NOTE THAT THIS CODE IS UNDER HEAVY DEVELOPMENT_**


The CASA Pipeline (CAPI) is a [CASA](http://casa.nrao.edu)-based higher-level layer to be executed on top of CASA (either within an interactive session or directly from the terminal) to help with the calibration of radio data.

This is a quite personal project which will be completely focused on the calibration of data from the European VLBI Network ([EVN](https://evlbi.org)) and from the Giant Metrewave Radio Telescope [GMRT](http://gmrt.ncra.tifr.res.in).


**CAPI has the following goals:**
- Simplify the user interaction layer with CASA to make it more dynamic and avoid command repetitions within a single session.
- Use the default and more optimized approach for data reduction for the different facilities, instead of a fully manual step-by-step approach.
- Still leave room for an interactive, manual improvement within the execution of CAPI.
- Get better diagnostic plots from the data and the calibration tables.



## Usage

The idea of CAPI is to allow the user to do something like the following for a standard data reduction:

```python
import casa_pipeline as capi

# This can be given as a dict or by using the template input file.
params = {'target': ['target1', 'target2'],
          'phaseref': ['pcal1', 'pcal1'],
          'fringefinder': ['ff'],
          'reference_antenna': 'EF',   # can be a list of antennas too
          'calibration': {'sbd_timerange': '2023/02/19/04:53:00~2023/02/19/04:55:30'}}

obs = capi.Project('obsname', params=params)

# Importing FITS-IDI
if not obs.ms.msfile.exists():
    obs.importdata.evn_fitsidi('obsname_1_1.IDI*', delete=False)
else:
    obs.calibrate.clearcal()
    obs.calibrate.callib.remove_entry('all')
    obs.flag.flagdata(mode='unflag')

# If you want to see some metadata from the observation:
obs.summary(outfile=f"results/summary-{obs.projectname}.md")
# Or indivually:
obs.ms.sources
obs.ms.sources.targets
obs.ms.sources.fringe_finders
obs.ms.antennas
obs.ms.antennas.observed  # will tell you which antennas have data
obs.ms.timerange
obs.ms.freqsetup


# # Initial flagging
obs.flag.apply_flags_from_file()
obs.flag.edge_channels()  # 0.1 (10% of the edges) by default
obs.flag.quack("WB", 4.0)  # 4 s on WB data
obs.flag.aoflagger()

# Initial calibration
obs.calibrate.a_priori_calibration()
# obs.calibrate.ionospheric_corrections()
obs.calibrate.main_calibration()
obs.calibrate.apply_calibration()

obs.flag.tfcrop()

obs.calibrate.backup()
obs.calibrate.recalibration()
obs.calibrate.apply_calibration()

obs.flag.tfcrop()

splits1 = obs.ms.split() # by default means sources=obs.sources.names

# Both split() and export_uvfits() prepare the data as Difmap-compatible
for src in obs.sources.phase_calibrators.names + obs.sources.targets.names:
    splits1[src].export_uvfits()


# steps for selfcal will follow
```

As you can infer, CAPI basically will call the standard functions in CASA. However, it adds the object-oriented higher layer that saves you from typing most of parts that are already known by the program.
For instance, instead of calling to:
```python
bandpass(msfile, **kwargs)
```
you will be able to run it as:
```python
obs.calibrate.bandpass(**kwargs)
```
This will make use of all the already-known information if not given by the kwargs. For example, it will use the MS file stored in `obs`, and it will use the fringe finder sources, as defined in `obs`, without need to type them again.





## External Dependencies

While CAPI is based primarily on CASA, it does not use exclusively such package. Instead, it uses different alternative programs which produce better and faster solutions for some specific parts.

- [jplotter](https://github.com/haavee/jiveplot): used to create plots from the data.
- [AOFlagger](https://aoflagger.readthedocs.io/en/latest/introduction.html): fast and accurate flagger to recognize radio frequency interference (RFI) in the data.
- [WSClean](https://wsclean.readthedocs.io/en/latest/index.html): fast and optimized tool to create cleaned maps.


## Python Dependencies

CAPI assumes that you already have [CASA](http://casa.nrao.edu) installed in your system. Independently of the type of installation (either the standard or modular), you will need to install some extra Python dependencies in the environment that CASA uses.

For a standard (self-contained) installation of CASA, this means that you will need to run `<path-to-casa-folder>/bin/python3 -m pip install -r requirements.txt` (probably as `sudo` if you have CASA in a protected folder). For a modular installation, just run the same `pip` command with the Python environment that CASA uses.

Then you can install CAPI with a similar procedure. After downloading it from GitHub, just place it under a directory included in the `$PYTHONPATH`.
CAPI will then be already visible for CASA as `import casa_pipeline as capi`.



## Limitations and assumptions

- CAPI is still under heave development and most of the functionality is still not implemented.
- CAPI is meant to be used for standard continuum, phase referencing, EVN/LBA and GMRT data, with no subarrays or multi-band data. These are thus the data types where CAPI will be tested.



## Acknowledgements

CAPI is not an original idea, but it only takes a different approach on existing programs to run such calibration.
In particular, CAPI has been motivated with code from:
- [rPicard](https://bitbucket.org/M_Janssen/Picard): Radboud Pipeline for the Calibration of high Angular Resolution Data, from M. Janssen.
- [VLBI_Pipeline](https://github.com/jradcliffe5/VLBI_pipeline): Generic CASA-based VLBI pipeline for use on job cluster, from J. Radcliffe.
- [CAPTURE-CASA6](https://github.com/ruta-k/CAPTURE-CASA6): CAsa Pipeline-cum-Toolkit for Upgrade GMRT data REduction, from K. Ruta.
- [CASA eMERLIN](https://github.com/e-merlin/eMERLIN_CASA_pipeline): Pipeline to calibrate data from the e-MERLIN array, by J. Moldón.

It also integrates the code from the [casa-vlbi](https://github.com/jive-vlbi/casa-vlbi) (from JIVE) repository to assist with EVN data compatibility with CASA.


