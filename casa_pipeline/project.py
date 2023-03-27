#!/usr/bin/env python3
"""Defines a VLBI experiment with all the relevant metadata required during the post-processing.
The metadata is obtained from the MS itself.
"""
import pickle
import logging
import datetime as dt
from pathlib import Path
from typing import Optional, Union # Iterable, NoReturn, List, Tuple
# import blessed
# tomli was introduced in the standard library as tomllib in Python 3.11
try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib

from . import obsdata
from . import calibration
from . import flagging
from . import plotting
from . import imaging

class Project(object):
    """Defines a (VLBI) observation with all relevant metadata.
    """
    @property
    def projectname(self) -> str:
        """Name of the project associated to the observation.
        """
        return self._projectname

    @property
    def cwd(self) -> Path:
        """Returns the working directory for the project.
        """
        return self._cwd

    @property
    def observatory(self) -> str:
        """Returns the name of the observatory that conducted the observations.
        """
        return self._observatory

    @observatory.setter
    def observatory(self, observatory_name: str):
        assert isinstance(observatory_name, str)
        self._observatory = observatory_name

    @property
    def epoch(self) -> dt.time:
        return self._ms.epoch

    @property
    def timerange(self) -> tuple:
        return self._ms.timerange

    @property
    def antennas(self) -> obsdata.Antennas:
        return self._ms.antennas

    @property
    def refant(self) -> list:
        return self._ms.refant

    @refant.setter
    def refant(self, new_refant: Union[str, list]):
        """Defines the antenna to use as reference during calibration.
        It can be either a str (with a single antenna, or comma-separated antennas),
        or a list with the antenna(s) to use as reference.
        """
        if isinstance(new_refant, str):
            self._ms.refant = [ant.strip() for ant in new_refant.split(',')]
        elif isinstance(new_refant, list):
            self._ms.refant = new_refant

    @property
    def freqsetup(self) -> obsdata.FreqSetup:
        return self._ms.freqsetup

    @property
    def sources(self) -> obsdata.Sources:
        return self._ms.sources

    @property
    def logdir(self) -> Path:
        """Returns the directory where the log outputs from the process are stored.
        """
        return self._logdir

    @property
    def caldir(self) -> Path:
        """Returns the directory where the calibration table files are stored.
        """
        return self._caldir

    @property
    def outdir(self) -> Path:
        """Returns the output directory where images, plots and other final products are stored.
        """
        return self._outdir

    @property
    def params(self) -> dict:
        return self._args

    @params.setter
    def params(self, new_params: dict):
        assert isinstance(new_params, dict), \
               f"{new_params} should be a dict but is {type(new_params)}"
        self._args = new_params

    @property
    def last_step(self):
        """Returns the last post-processing step that did run properly in a tentative previous run.
        """
        return self._last_step

    @last_step.setter
    def last_step(self, last_step):
        self._last_step = last_step

    @property
    def logger(self) -> logging.Logger:
        return self._logger

    @property
    def ms(self):
        return self._ms

    @property
    def calibrate(self):
        return self._calibration

    @property
    def flag(self):
        return self._flagging

    @property
    def plot(self):
        return self._plotting

    @property
    def image(self):
        return self._imaging

    @property
    def importdata(self):
        return self._importing

    def __init__(self, projectname: str, observatory: Optional[str] = '',
                 params: Union[dict, str, Path] = None,
                 cwd: Optional[Union[Path, str]] = None, logging_level=logging.INFO):
        """Initializes an EVN experiment with the given name.

        Inputs:
            projectname : str
               The name of the project (case insensitive).

            observatory : str  [default = '']
                Name of the observatory associated to this project.
                If not specified, it will be read from the MS when exists.

            params : dict/Path/str  [optional]
                Default parameters read from the pipeline input file and are loaded here.
                If it is a Path or str, it is assumed to be the file name containing the
                input parameters.

            cwd : Path or str  [optional. Default = $PWD]
                The default directory where the pipeline will run, search for and create the files.
        """
        assert (projectname != '') and (isinstance(projectname, str)), \
               "The project name needs to be a non-empty string."
        self._projectname = projectname
        self._observatory = observatory
        logging.basicConfig(level=logging_level, format="%(asctime)s %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M")
        self._logger = logging.getLogger(self.projectname)
        self.logger.debug(f"Initializing Project '{self.projectname}' from observatory " \
                          f"'{self.observatory}'.")
        if cwd is None:
            self._cwd = Path('./')
        elif isinstance(cwd, str):
            self._cwd = Path(cwd)
        elif isinstance(cwd, Path):
            self._cwd = cwd
        else:
            TypeError(f"The working directory ({cwd}) should be either None, Path or str.")

        self._logdir = self.cwd / 'log'
        self._caldir = self.cwd / 'caltables'
        self._outdir = self.cwd / 'results'
        self._last_step = None
        self._args = {}

        for a_dir in (self.cwd, self.logdir, self.caldir, self.outdir):
            a_dir.mkdir(exist_ok=True)


        self._local_copy = self.cwd / Path(f"{self.projectname.lower()}.obj")
        if self.exists_local_copy():
            self.load(self._local_copy)

        if params is not None:
            if isinstance(params, dict):
                self._args = params
            elif isinstance(params, Path) or isinstance(params, str):
                with open(params, "rb") as fp:
                    self._args = tomllib.load(fp)

        if observatory != '':
            self._observatory = observatory
            self.logger.debug(f"Observatory changed to {self.observatory}.")

        if not self.exists_local_copy():
            self.logger.debug("Project Initializing from scratch.")
            self._ms = obsdata.Ms(projectname, observatory=self.observatory, cwd=self.cwd,
                                  params=self._args, logger=self._logger)
            # TODO: differentiate as function of the observatory
            self._calibration = calibration.Calibration(self._ms, self.caldir)
            self._flagging = flagging.Flagging(self._ms)
            self._plotting = plotting.Plotting(self._ms)
            self._imaging = imaging.Imaging(self._ms)
            self._importing = obsdata.Importing(self._ms)
        elif params is not None:
            self._calibration._ms._params = self._args
            self._flagging._ms._params = self._args
            self._plotting._ms._params = self._args
            self._imaging._ms._params = self._args
            self._importing._ms._params = self._args

        # self._last_step = None
        self.store()


    def exists_local_copy(self):
        """Checks if there is a local copy of the Experiment object stored in a local file.
        """
        return self._local_copy.exists()


    def store(self, path : Optional[Union[Path, str]] = None):
        """Stores the current Experiment into a file in the indicated path. If not provided,
        it will be '.{projectname.lower()}.obj' where exp is the name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'wb') as file:
            pickle.dump(self, file)

        self.logger.debug(f"Local copy of {self.projectname} stored at {self._local_copy}.")


    def load(self, path: Optional[Union[Path, str]] = None) -> bool:
        """Loads the current Experiment that was stored in a file in the indicated path.
        If path is None, it assumes the standard path of '.{exp}.obj', where exp is the
        name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'rb') as file:
            obj = pickle.load(file)

        self._observatory = obj._observatory
        self._logger = obj._logger
        self._args = obj._args
        self._last_step = obj._last_step
        self._ms = obj._ms
        self._calibration = obj._calibration
        self._flagging = obj._flagging
        self._plotting = obj._plotting
        self._imaging = obj._imaging
        self._importing = obj._importing
        self.logger.info(f"Loaded Project object from the stored local copy at {self._local_copy}.")
        return True

    def summary(self, outfile: Union[str, Path, None] = None):
    #     TODO:  this is commented out because I have issues installing blessed within the CASA
    #     environment in my desktop
    #     term = blessed.Terminal(force_styling=True)
        s_file = []
    #     with term.fullscreen(), term.cbreak():
    #         s = term.white_on_red(term.center(term.bold(f"--- {self.projectname} ---")))
        s_file += [f"# {self.projectname}"]
    #         s += f"{term.normal}\n\n{term.normal}"
    #         s += term.bold_green('General Information\n')
    #         s += ["## General Information"]
        if self.observatory != '':
    #             s += term.bright_black('Observatory: ') + f"{self.observatory}\n"
            s_file += [f"Observatory: {self.observatory}"]
    #
    #         # TODO: if it is EVN, put a direct link to the archive (same function as in evn_postprocess)
    #         s += term.bright_black('') + f"\n"
    #         s_file += [f""]
    #         s += term.bright_black('') + f"\n"
    #         s_file += [f""]
    #         s += term.bold_green('SOURCES\n')
        s_file += [f"Central frequency: {self.freqsetup.frequency:.2}"]
        s_file += [f"With a bandwith of {self.freqsetup.bandwidth.to(u.MHz)} divided in " \
                   f"{self.freqsetup.n_subbands} x {self.freqsetup.bandwidth_per_subband.to(u.MHz)} " \
                   f"subbands."]
        s_file += [f"{self.freqsetup.channels} spectral channels per subband.\n"]

        s_file += ["## Sources"]
    #         s += term.bright_black('Fringe finders: ') + \
    #              f"{', '.join([s.name for s in self.sources.fringe_finders])}\n"
        s_file += ["Fringe finders: " \
                   f"{', '.join([s.name for s in self.sources.fringe_finders])}"]
    #         s += term.bright_black('Phase calibrators: ') + \
    #              f"{', '.join([s.name for s in self.sources.phase_calibrators])}\n"
        s_file += ["Phase calibrators: " \
                   f"{', '.join([s.name for s in self.sources.phase_calibrators])}"]
    #         s += term.bright_black('Target sources: ') + \
    #              f"{', '.join([s.name for s in self.sources.targets])}\n\n"
        s_file += [f"Target sources: {', '.join([s.name for s in self.sources.targets])}\n"]
        longest_src_name = max([len(s) for s in self.sources.names])
        for src in self.sources:
    #             s += term.bright_black(f"{src.name}: ") + f"src.coordinates.\n"
            s_file += [f"{src.name}:{' '*(longest_src_name-len(src.name))} {src.coordinates.to_string('hmsdms')}"]
    #
    #
        s_file += ["\n## Antennas"]
        s_file += ["   Did Observe?  Subbands"]
        for ant in self.antennas:
            s_file += [f"{ant.name} {'yes' if ant.observed else 'no '} " \
                       f"{' '*(3*(ant.subbands[0]))}{ant.subbands}"]

        if outfile is not None:
            with open(outfile, 'w') as fout:
                fout.write('\n'.join(s_file))

        print('\n'.join(s_file))

        # s_file += ["## Data files"]
    #         s += term.bright_black('') + f"\n"
    #         s += term.bright_black('') + f"\n"
    #         s += term.bright_black('') + f"\n"
    #         s += term.bright_black('') + f"\n"
    #         s += term.bright_black('') + f"\n"
    #         s += term.bright_black('') + f"\n"
    #         s += term.bright_black('') + f"\n"


    def __repr__(self, *args, **kwargs):
        rep = super().__repr__(*args, **kwargs)
        rep.replace("object", f"object ({self.projectname})")
        return rep

    def __str__(self):
        return f"<Project {self.projectname}>"
