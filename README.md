# EVN CASA Pipeline

Pipeline to analyze radio data from the European VLBI Network (EVN) using the NRAO CASA software package.


This pipeline is in a very early design stage, so any question, comment, feature requests, and suggestions, please be sent to marcote@jive.eu.

The Pipeline will be executed by loading a input file that will contain all the relevant information. All paramenters from each step can be overwritten if desired, making the pipeline quite flexible for any purpose.




## Pipeline Workflow

This is the current workflow of the EVN CASA Pipeline and all steps that will do.

We assume that data are in a Measurement Set (MS) format. Steps in italic will be optional (either because they are done internally at JIVE or depending on the datasets).

- _Include Tsys info to MS_.
- _Create gain curve corrections_.
- MS inspection (metadata, plots, observing information).
- Apply a-priori amplitude calibration.
- _Ionospheric corrections_.
- _Instrumental delay_.
- Global fringe.
- Bandpass correction.
- Split.
- Imaging/self-calibration for the calibrators.
- _Apply corrections to all antennas/sources_ (how to do it in CASA? Supported?).
- Imaging target.
- Produce final uvdata files.
- Produce final plots.

The configuration file (`config.inp`) will allow the users to overwrite any parameter used in any of the CASA tasks. The pipeline will thus be completely flexible in terms of calibration steps.


## Class Diagram

You can find an UML Class Diagram of the project in *diagram-uml.dia* created with the [Gnome Dia](https://wiki.gnome.org/Apps/Dia) program. A workflow of the pipeline steps can also be found in *diagram-workflow.dia*.



## File structure during pipeline execution

When running the EVN CASA Pipeline several files will be read and other ones will be created.


The basic files required to run the pipeline will be:
- *config.inp*: input file containing all commands to run the pipeline.
- *MS dataset*: (specified in the `config.inp`) MS file with the raw EVN data.

The pipeline will create several files and directories (under the output directory specified in `config.inp`):

- calibration_tables/: directory where all calibration tables will lie.
- cycle_n/: (with n=0,1,2,..) directories where output from the different calibration cycles will be stored.







