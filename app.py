#! /usr/bin/env python3
"""
"""
import os
import sys
import glob
import time
import argparse
import traceback
from pathlib import Path
import datetime as dt
import json
from evn_pipeline import project
from evn_pipeline import preprocess
from evn_pipeline import flagging
from evn_pipeline import calibration
# from .evn_pipeline import


# Rename the file to __main__.py. Then it can be executed by python -m evn_postprocess

# {path-to-casa}/bin/mpicasa -n *N* {path-to-casa}/bin/casa    {usual options if needed}

__version__ = 0.1
__prog__ = 'casa_pipeline.py'
usage = "%(prog)s  [-h]  <experiment_name>  <support_scientist>\n"
description = """XXXXX. #TODO
"""

def main():
    all_steps = {'import_fits': [preprocess.append_tsys, preprocess.generate_gain_curve, preprocess.convert_uvflg,
                                 preprocess.importfitsidi2ms],
                 'summary': [preprocess.get_ms_info, preprocess.listobs],
                 'apriori_cal': [calibration.a_priori_calibration, flagging.apply_apriori_flagging,
                                 flagging.edge_channels], # add aoflagger
                 'ionos': [calibration.ionospheric_corrections],
                 'main_cal': [calibration.main_calibration, calibration.recalibration],
                 'apply_cal': [calibration.apply_calibration],
                 # 'plots': [plotting.plotcals, plotting.plotdata],
                 'split': [calibration.split, calibration.export_uvfits]
                }
    parser = argparse.ArgumentParser(description=description, prog=__prog__, usage=usage)
    parser.add_argument('input_file', type=str, help=f"Path to the input file to be used to run {__prog__}.")
    parser.add_argument('-f', '--first', type=str, default=None, help="First step in the data reduction to run.")
    parser.add_argument('-l', '--last', type=str, default=None, help="Last step in the data reduction to run.")
    parser.add_argument('-i', '--ignore', default=False, action='store_true',
                        help="In case a previous run was conducted, will not resume from the last step that run there.")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    args = parser.parse_args()

    inputs = Path(args.input_file)
    if not inputs.exists():
        print(f"\n\nThe input file {inputs.name} has not been found.")
        sys.exit(1)

    inp_pars = json.load(open(inputs, 'r'))
    the_project = project.Project(name=inp_pars['globals']['project_code'], targets=inp_pars['globals']['targets'],
                          bandpass_calibrators=inp_pars['globals']['bandpass_calibrators'],
                          phase_calibrators=inp_pars['globals']['phase_calibrators'],
                          reference_antennas=inp_pars['globals']['reference_antennas'],
                          instr_delay_timerange=inp_pars['globals']['instr_delay_timerange'])
    if inp_pars['globals']['fits_idi_files'] is not None:
        the_project.idi_files = sorted(glob.glob(inp_pars['globals']['fits_idi_files']))
    else:
        the_project.idi_files = []

    # If there is already a project stored with this name, then it will retrieve the last step.
    # It will still update the parameterr from the input file as those could have been updated
    # TODO: this is a loop that can be completely removed!
    if the_project.local_copy.exists():
        the_project_old = the_project.load()
        the_project_old.targets = inp_pars['globals']['targets']
        the_project_old.bandpass_calibrators = inp_pars['globals']['bandpass_calibrators']
        the_project_old.phase_calibrators = inp_pars['globals']['phase_calibrators']
        the_project_old.refants = inp_pars['globals']['reference_antennas']
        the_project_old.instr_delay_timerange = inp_pars['globals']['instr_delay_timerange']
        the_project = the_project_old


    step_keys = list(all_steps.keys())

    first_step = step_keys.index(args.first) if (args.first is not None) else 0
    last_step = step_keys.index(args.last) if (args.last is not None) else len(step_keys) - 1

    if args.ignore or (the_project.last_step is None):
        the_steps = step_keys[first_step:last_step+1]
    else:
        assert step_keys.index(the_project.last_step) <= last_step, "The pipeline already run up to a step farther" \
                f" than {args.last}. Specify a different one or run with the '--ignore' flag."
        first_step = max(first_step, step_keys.index(the_project.last_step))
        the_steps = step_keys[first_step:last_step+1]
        print(f"Starting the data reduction at the step {the_steps[0]} to {the_steps[-1]}.")
        time.sleep(2)

    try:
        for a_step in the_steps:
            for a_task in all_steps[a_step]:
                a_task(the_project)

            the_project.last_step = a_step

        print('\nThe data reduction has finished properly.')
    except :
        for a_logpath in glob.glob('casa-*.log'):
            a_logfile = Path(a_logpath)
            a_logfile.replace(the_project.logdir / a_logfile.name)

        traceback.print_exc()
        sys.exit(1)
    finally:
        the_project.store()

if __name__ == '__main__':
    main()

