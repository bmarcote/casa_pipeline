#!/usr/bin/env python3
"""Defines a VLBI experiment with all the relevant metadata required during the post-processing.
The metadata is obtained from the MS itself.
"""
import pickle
from pathlib import Path
from typing import Optional, Union # Iterable, NoReturn, List, Tuple
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

    def __init__(self, projectname: str, observatory: Optional[str] = '', params: dict = None,
                 cwd: Optional[Union[Path, str]] = None):
        """Initializes an EVN experiment with the given name.

        Inputs:
            projectname : str
               The name of the project (case insensitive).

            observatory : str  [default = '']
                Name of the observatory associated to this project.
                If not specified, it will be read from the MS when exists.

            params : dict  [optional]
                Default parameters read from the pipeline input file and are loaded here.

            cwd : Path or str  [optional. Default = $PWD]
                The default directory where the pipeline will run, search for and create the files.
        """
        assert (projectname != '') and (isinstance(projectname, str)), \
            "The project name needs to be a non-empty string."
        self._projectname = projectname
        self._observatory = observatory
        # logpath = Path("./logs")
        # logpath.mkdir(parents=True, exist_ok=True)
        # self._logs = {'dir': logpath, 'file': Path("./processing.log)")}
        # if 'cwd' in params:
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

        for a_dir in (self.cwd, self.logdir, self.caldir, self.outdir):
            a_dir.mkdir(exist_ok=True)

        self._local_copy = self.cwd / Path(f"{self.projectname.lower()}.obj")
        if self.exists_local_copy():
            self = self.load(self._local_copy)

        # self._last_step = None
        if params is None:
            self._args = {}
        else:
            self._args = params

        self._ms = obsdata.Ms(projectname, cwd)
        # TODO: differentiate as function of the observatory
        self._calibration = calibration.Calibration(self._ms)
        self._flagging = flagging.Flagging(self._ms)
        self._plotting = plotting.Plotting(self._ms)
        self._imaging = imaging.Imaging(self._ms)
        self._importing = obsdata.Importing(self._ms)

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

    def load(self, path: Optional[Union[Path, str]] = None):
        """Loads the current Experiment that was stored in a file in the indicated path.
        If path is None, it assumes the standard path of '.{exp}.obj', where exp is the
        name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'rb') as file:
            obj = pickle.load(file)

        return obj

    def __repr__(self, *args, **kwargs):
        rep = super().__repr__(*args, **kwargs)
        rep.replace("object", f"object ({self.projectname})")
        return rep

    def __str__(self):
        return f"<Project {self.projectname}>"
