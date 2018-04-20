# EVN Pipeline calibration with CASA


## VLA Pipeline

### Description in the website

- Loads the data into from the archival format (Science Data Model-Binary Data Format [SDM-BDF]) to a CASA MeasurementSet (MS), applies Hanning smoothing, and obtains information about the observing set-up from the MS;
- Applies online flags and other deterministic flags (shadowed data, edge channels of sub-bands, etc.);
- Prepares models for primary flux density calibrators;
- Derives pre-determined calibrations (antenna position corrections, gain curves, atmospheric opacity corrections, requantizer gains, etc.);
- Iteratively determines initial delay and bandpass calibrations, including flagging of RFI (Radio Frequency Interference) and some automated identification of system problems;
- Derives initial gain calibration, and derives the spectral index of the bandpass calibrator;
- RFI flagging is done on data with the initial calibration applied;
- Derives final delay, bandpass, and gain/phase calibrations, and applies them to the data;
- Runs the RFI flagging algorithm on the target data;
- Calculates data weights based on the inverse square RMS noise of the MS;
- Creates diagnostic images of calibrators.
 

### Flow chart (VLA CASA pipeline)

- load main ms
- extract info about contents

- hanning smooth the data

- apply online flags
- deterministic flags

- set flux cal models (setjy)

- derive prior cals (antenna pos. errors, opacity, gains,...)

- test delay & BP cals

- falg bad spws

- rfalg on calibrated delay & BP cals

- semifinal delay & BP cals (not normalized, spectral index of BP not yet determined)

- split out calibrators with spw avg into calibrators.ms
- do a test gaincal to determine short and long solints

- make amp gain table for flux density bootstrapping
- Do flux density bootstrapping (derive spectral index of calibrators and re-inserting into main MS)

- make final cal tables (delay, BP, amp, phase)
- redo split, fluxboot, etc.

- apply final calibrations to main MS

- run rflag on all fields including target

- run statwt to set weights using rms per spw

- make some diagnostic plots

- prepare list of images to make

- make images of calibrators per spw

- science quality imaging



## Picart Pipeline (poster from Michael Janssen)

- Automatic flagging heuristics based on autocorrelations

- Correct field rotation angle
- Scalar bandpass
- a-priori calibration
- correct for digital sampling (ACCOR)

- manual phase calibration (instrumental phase and delay per spw using bright calibrator)

- global fringe-fit

- complex bandpass

- polarization calibration (RL delay and phase offsets, D-terms)

- calbirated dataset




## Steps to implement in the EVN CASA Pipeline

### Inspection

plotants(vis= , figfile=)


### Ionosphere corrections

- Import the TEC image
from recipes import tec_maps
tec_image, tec_rms_image, tec_graph = tec_maps.create(vis=msdata, doplot=True)

gencal(vis=msdata, caltype='tecim', caltable='tecim.cal', infile=tec_image)


## eMERLIN CASA pipeline

### How to run the pipeline

casa -c code_pipeline.py -i <input_file>

or for a parallelized version:

mpicasa -n <num_cores> casa -c code_pipeline.py -i <input_file>

run_in_casa = True
pipeline_path = '/path/to/pipeline_path/'   # You need to define this variable explicitly
execfile(pipeline_path + 'eMERLIN_CASA_pipeline.py')
inputs, msinfo = run_pipeline(inputs_path=<input file>)





## Current CASA Data Reduction for EVN data


in plotms:  overwrite=True, showgui=False,...

