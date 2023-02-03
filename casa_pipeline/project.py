#!/usr/bin/env python3
"""Defines a VLBI experiment with all the relevant metadata required during the post-processing.
The metadata is obtained from the MS itself.
"""
from __future__ import annotations
import os
import sys
import pickle
import glob
import json
import subprocess
import datetime as dt
import numpy as np
from pathlib import Path
from typing import Optional, Iterable, NoReturn, List, Union, Tuple
from dataclasses import dataclass
from collections import defaultdict
from natsort import natsort_keygen
from pyrap import tables as pt
import casatasks
from enum import Enum
from astropy import units as u
from rich import print as rprint
from rich import progress
# import blessed
from astropy import coordinates as coord
from casa_pipeline.casavlbitools import fitsidi
from . import calibration
from . import flagging
from . import imaging
from . import plotting


def chunkert(counter: int, max_length: int, increment: int) -> Tuple[int, int]:
    """Silly function to select a subset of an interval in
       [counter, counter + increment] : 0 < counter < max_length.

    Yields the tuple (counter, + interval_increment) : interval_increment = min(increment, max_length - counter))
    """
    while counter < max_length:
        this_increment = min(increment, max_length - counter)
        yield (counter, this_increment)
        counter += this_increment

def percentage(x, y):
    """Returns the percentage value of  100 * x / y.
    """
    return (x / y)*100.0


class SourceType(Enum):
    """Type of source (target, calibrator, fringefinder, or other)
    """
    target = 0
    calibrator = 1
    fringefinder = 2
    other = 3


class Source(object):
    """Defines a source by name, type (i.e. target, reference, fringefinder, other)
    and if it must be protected or not (password required to get its data).
    """
    @property
    def name(self) -> str:
        """Name of the source.
        """
        return self._name

    @property
    def coordinates(self) -> coord.SkyCoord:
        """Coordinates of the source, as astropy.coordinate.SkyCoord type.
        """
        return self._coordinates

    @property
    def type(self) -> SourceType:
        """Type of the source (SorceType object).
        """
        return self._type

    @property
    def model(self) -> Path:
        """Source model file (generated after imaging or other tasks and can be used for calibration).
        """
        return self._model

    @model.setter
    def model(self, new_model: Union[Path, str]):
        if isinstance(new_model, str):
            new_model = Path(new_model)

        self._model = new_model

    def __init__(self, name: str, sourcetype: SourceType, coordinates: coord.SkyCoord,
                 model: Optional[Union[Path, str]] = None):
        assert isinstance(name, str), f"The name of the source must be a string (currrently {name})"
        assert isinstance(sourcetype, SourceType), \
               f"The name of the source must be a SourceType object (currrently {sourcetype})"
        assert isinstance(coordinates, coord.SkyCoord), "The coordinates of the source must " \
               f"be an astropy.coordinates.SkyCoord object (currently {coordinates})"
        self._name = name
        self._type = sourcetype
        self._coordinates = coordinates
        self.model = model

    def __iter__(self):
        for key in ('name', 'type', 'protected'):
            yield key, getattr(self, key)

    def __str__(self):
        return f"Source({self.name}, {self.type.name}, at {self.coordinates.to_string('hmsdms')})"

    def __repr__(self):
        return self.__str__()

    def json(self):
        """Returns a dict with all attributes of the object.
        I define this method to use instead of .__dict__ as the later only reporst
        the internal variables (e.g. _username instead of username) and I want a better
        human-readable output.
        """
        d = dict()
        for key, val in self.__iter__():
            if isinstance(val, SourceType):
                d[key] = val.name
            elif isinstance(val, coord.SkyCoord):
                d[key] = val.to_string('hmsdms')
            else:
                d[key] = val

        return d


class Sources(object):
    """Defines a collection of sources.
    """
    def __init__(self, sources: Optional[Iterable[Source]] = None):
        if sources is not None:
            self._sources = sources[:]
        else:
            self._sources = []

    def append(self, new_source: Source) -> NoReturn:
        """Adds a new source to the collection of sources.

        Note that this function and "add" are completely equivalent.
        """
        assert isinstance(new_source, Source)
        self._sources.append(new_source)

    def add(self, new_source: Source) -> NoReturn:
        """Adds a new source to the collection of sources.

        Note that this function and "add" are completely equivalent.
        """
        self.append(new_source)

    @property
    def names(self) -> List[str]:
        """Returns the names of all sources included in Sources.
        """
        return [s.name for s in self._sources]

    @property
    def types(self) -> List[SourceType]:
        """Returns the types of all sources in Sources.
        """
        return [s.type for s in self._sources]

    @property
    def coordinates(self) -> List[coord.SkyCoord]:
        """Returns the coordinates of all sources in Sources.
        """
        return [s.coordinates for s in self._sources]

    @property
    def targets(self) -> Sources:
        """Returns the target sources.
        """
        return Sources([s for s in self._sources if s.type == SourceType.target])

    @property
    def calibrators(self) -> Sources:
        """Returns the phase calibrator sources.
        """
        return Sources([s for s in self._sources if s.type == SourceType.calibrator])

    @property
    def other_sources(self) -> Sources:
        """Returns the sources catalogued as "others".
        """
        return Sources([s for s in self._sources if s.type == SourceType.other])

    @property
    def fringe_finders(self) -> Sources:
        """Returns the fringe-finder sources.
        """
        return Sources([s for s in self._sources if s.type == SourceType.fringefinder])

    def __len__(self):
        return len(self._sources)

    def __getitem__(self, key):
        return self._sources[self.names.index(key)]

    def __delitem__(self, key):
        return self._sources.remove(self.names.index(key))

    def __iter__(self):
        return self._sources.__iter__()

    def __reversed__(self):
        return self._sources[::-1]

    def __contains__(self, key):
        return key in self.names

    def __str__(self):
        return f"Sources([{', '.join(self.names)}])"

    def __repr__(self):
        return f"Sources([{', '.join(self.names)}])"

    def json(self):
        """Returns a dict with all attributes of the object.
        I define this method to use instead of .__dict__ as the later only reporst
        the internal variables (e.g. _username instead of username) and I want a better
        human-readable output.
        """
        d = dict()
        for ant in self.__iter__():
            d['Antenna'] = ant.__dict__

        return d


@dataclass
class Antenna:
    """Defines an antenna.
    It has three parameters:
        name : str
            Name of the antenna
        observed : bool
            If the antenna has observed (has no-null data).
        subbands : tuple
            Tuple with the subbands numbers where the antenna observed.
            It may be all subbands covered in the observation or a subset of them.
    """
    name: str
    observed: bool = True
    subbands: tuple = ()


class Antennas(object):
    """List of antennas (Antenna class)
    """
    def __init__(self, antennas: Optional[Iterable[Antenna]] = None):
        if antennas is not None:
            self._antennas = antennas[:]
        else:
            self._antennas = []

    def append(self, new_antenna: Antenna) -> NoReturn:
        """Adds a new antenna to the collection of antennas.

        Note that this function and "add" are completely equivalent.
        """
        self.add(new_antenna)

    def add(self, new_antenna: Antenna) -> NoReturn:
        """Adds a new antenna to the collection of antennas.

        Note that this function and "add" are completely equivalent.
        """
        assert isinstance(new_antenna, Antenna)
        self._antennas.append(new_antenna)

    @property
    def names(self) -> List[str]:
        """Returns the name of the antenna
        """
        return [a.name for a in self._antennas]

    @property
    def observed(self) -> List[tuple]:
        """Returns the antenna names that have data (have observed) in the observation
        """
        return [a.name for a in self._antennas if a.observed]

    @property
    def subbands(self):
        """Returns the subbands that each antenna observed.
        """
        return [a.subbands for a in self._antennas if a.observed]

    def __len__(self):
        return len(self._antennas)

    def __getitem__(self, key):
        return self._antennas[self.names.index(key)]

    def __delitem__(self, key):
        return self._antennas.remove(self.names.index(key))

    def __iter__(self):
        return self._antennas.__iter__()

    def __reversed__(self):
        return self._antennas[::-1]

    def __contains__(self, key):
        return key in self.names

    def __str__(self):
        return f"Antennas([{', '.join([ant if (ant in self.observed) else '['+ant+']' for ant in self.names])}])"

    def __repr__(self):
        return f"Antennas([{', '.join([ant if ant in self.observed else '['+ant+']' for ant in self.names])}])"

    def json(self):
        """Returns a dict with all attributes of the object.
        I define this method to use instead of .__dict__ as the later only reporst
        the internal variables (e.g. _username instead of username) and I want a better
        human-readable output.
        """
        d = dict()
        for ant in self.__iter__():
            d['Antenna'] = ant.__dict__

        return d


class FreqSetup(object):
    """Defines the frequency setup of a given observation with the following data:
        - n_subbands :  int
            Number of subbands.
        - channels : int
            Number of channels per subband.
        - frequencies : array-like
            Reference frequency for each channel and subband (NxM array, with N
            number of subbands, and M number of channels per subband).
        - central_freq : astropy.units.Quantity
            Central frequency of the observation.
        - bandwidths : astropy.units.Quantity or float
            Total bandwidth for each subband.
    """
    @property
    def n_subbands(self) -> int:
        """Returns the number of subbands (spectral windows) of the observation.
        """
        return self._n_subbands

    @property
    def channels(self) -> int:
        """Returns the number of channels per subband of the observation.
        Assumes that all subbands have the same number of channels.
        """
        return self._channels

    @property
    def frequencies(self) -> np.ndarray:
        """Returns the frequencies of each channel and subband in a NxM array, with
        N the number of subbands, and M the number of channels per subband.
        """
        return self._freqs

    @property
    def frequency(self) -> u.Quantity:
        """Returns the central frequency of the observation, in a astropy.units.Quantity format.
        """
        return (np.mean(self.frequencies)*u.Hz).to(u.GHz)

    @property
    def bandwidth_per_subband(self) -> u.Quantity:
        """Returns the total bandwidth per subband
        """
        return self._bandwidth_spw

    @property
    def bandwidth(self) -> u.Quantity:
        """Returns the total bandwidth of the observation
        """
        return self._bandwidth_spw * self.n_subbands

    @property
    def frequency_range(self) -> tuple[u.Quantity, u.Quantity]:
        """Returns the lowest and highest frequency observed.
        """
        return (self.frequency - self.bandwidth/2, self.frequency + self.bandwidth/2)

    def __init__(self, channels: int, frequencies, bandwidth_spw: Union[float, u.Quantity]):
        """Inputs:
            - chans : int
                Number of channels per subband.
            - freqs : array-like
                Reference frequency for each channel and subband (NxM array, M number
                of channels per subband.
            - bandwidth_spw : float or astropy.units.Quantity
                Total bandwidth for each subband. If not units are provided, Hz are assumed.
        """
        self._n_subbands = frequencies.shape[0]
        assert isinstance(channels, (int, np.int32, np.int64)), \
            f"Chans {channels} is not an int as expected (found type {type(channels)})."
        assert isinstance(bandwidth_spw, float) or isinstance(bandwidth_spw, u.Quantity), \
            f"Bandwidth {bandwidth_spw} is not a float or Quantity as expected (found {type(bandwidth_spw)})."
        assert frequencies.shape == (self._n_subbands, channels)
        self._channels = int(channels)
        self._freqs = np.copy(frequencies)
        if isinstance(bandwidth_spw, float):
            self._bandwidth_spw = bandwidth_spw*1e-9*u.GHz
        else:
            self._bandwidth_spw = bandwidth_spw

    def __iter__(self):
        for key in ('n_subbands', 'channels', 'bandwidth_spw', 'frequencies'):
            yield key, getattr(self, key)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"{self.frequency} central frequency.\n{self.n_subbands} x {self.bandwidth_per_subband} " \
               f"subbands. {self.channels} channels per subband."

    def json(self):
        """Returns a dict with all attributes of the object.
        I define this method to use instead of .__dict__ as the later only reporst
        the internal variables (e.g. _username instead of username) and I want a better
        human-readable output.
        """
        d = dict()
        for key, val in self.__iter__():
            if isinstance(val, u.Quantity):
                d[key] = val.to(u.Hz).value
            elif isinstance(val, np.ndarray):
                d[key] = list(val)
            else:
                d[key] = val

        return d



class Project(object):
    """Defines a (VLBI) observation with all relevant metadata.
    """
    @property
    def projectname(self) -> str:
        """Name of the project associated to the observation.
        """
        return self._expname

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
    def cwd(self) -> Path:
        """Returns the working directory for the project.
        """
        return self._cwd

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
        assert isinstance(new_params, dict), f"{new_params} should be a dict. Instead is {type(new_params)}"
        self._args = new_params

    @property
    def last_step(self):
        """Returns the last post-processing step that did run properly in a tentative previous run.
        """
        return self._last_step

    @last_step.setter
    def last_step(self, last_step):
        self._last_step = last_step

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
        assert projectname != '', "The projectname cannot be an empty string"
        self._expname = projectname
        self.observatory = observatory
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
            TypeError(f"The working directory ({cwd}) should be either None, a Path type, or a str.")

        self._logdir = self.cwd / 'log'
        self._caldir = self.cwd / 'caltables'
        self._outdir = self.cwd / 'results'

        for a_dir in (self.cwd, self.logdir, self.caldir, self.outdir):
            a_dir.mkdir(exist_ok=True)

        self._local_copy = self.cwd / Path(f"{self.projectname.lower()}.obj")
        # self._last_step = None
        # self._args = params


    def exists_local_copy(self):
        """Checks if there is a local copy of the Experiment object stored in a local file.
        """
        return self._local_copy.exists()

    def store(self, path : Optional[Union[Path, str]] = None):
        """Stores the current Experiment into a file in the indicated path. If not provided,
        it will be '.{expname.lower()}.obj' where exp is the name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'wb') as file:
            pickle.dump(self, file)

    def load(self, path: Optional[Union[Path, str]] = None):
        """Loads the current Experiment that was stored in a file in the indicated path. If path is None,
        it assumes the standard path of '.{exp}.obj' where exp is the name of the experiment.
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



class Ms(Project):
    """Defines one MS object for a given observation.
    It contains all relevant information that is MS-depended, e.g. associated MS file,
    frequency setup, sources, etc. But also objects associated to the calibration tasks that
    can run in the MS.
    """
    @property
    def msfile(self) -> Path:
        """Returns the name of the MS file (libpath.Path object) associated to this correlator pass.
        """
        return self._msfile

    @msfile.setter
    def msfile(self, new_msfile: Union[Path, str]):
        if isinstance(new_msfile, Path):
            self._msfile = new_msfile
        elif isinstance(new_msfile, str):
            self._msfile = Path(new_msfile)
        else:
            TypeError(f"msfile() expected a str or Path type. Found {new_msfile} "
                      f"({type(new_msfile)}).")

    @property
    def uvfitsfile(self) -> Path:
        """Returns the UV FITS file associated to this observation, if any.
        """
        return self._uvfitsfile

    @uvfitsfile.setter
    def uvfitsfile(self, new_uvfitsfile: Union[Path, str]):
        if isinstance(new_uvfitsfile, Path):
            self._uvfitsfile = new_uvfitsfile
        elif isinstance(new_uvfitsfile, str):
            self._uvfitsfile = Path(new_uvfitsfile)
        else:
            TypeError(f"uvfitsfile() expected a str or Path type. Found {new_uvfitsfile} "
                      f"({type(new_uvfitsfile)}).")

    @property
    def sources(self) -> Sources:
        """List of sources present in this correlator pass.
        """
        return self._sources

    @sources.setter
    def sources(self, list_of_sources: Sources):
        assert isinstance(list_of_sources, Sources)
        self._sources = list_of_sources

    @property
    def antennas(self) -> Antennas:
        """List of antennas available in the experiment.
        """
        return self._antennas

    @antennas.setter
    def antennas(self, new_antennas: Antennas):
        isinstance(new_antennas, Antennas)
        self._antennas = new_antennas

    @property
    def refant(self) -> list:
        return self._refant

    @refant.setter
    def refant(self, new_refant):
        """Defines the antenna to use as reference during calibration.
        It can be either a str (with a single antenna, or comma-separated antennas), or a list with the
        antenna(s) to use as reference.
        """
        if isinstance(new_refant, str):
            new_refant = [ant.strip() for ant in new_refant.split(',')]
        elif not isinstance(new_refant, list):
            raise TypeError(f"The antenna introduced must be either a  str or list")

        for antenna in new_refant:
            assert antenna in self.antennas.names, f"Antenna {antenna} for refant is not in the array."
        self._refant = list(new_refant)

    @property
    def freqsetup(self) -> FreqSetup:
        """Returns the frequency setup of the observation (FreqSetup object).
        """
        return self._freqsetup

    @freqsetup.setter
    def freqsetup(self, freqsetup: FreqSetup):
        """Sets the frequency setup for the given correlator pass.
        """
        assert isinstance(freqsetup, FreqSetup)
        self._freqsetup = freqsetup

    @property
    def epoch(self) -> dt.date:
        """Epoch at which the project was observed (starting date).
        """
        return self._obsdate

    @epoch.setter
    def epoch(self, obsdate: dt.date):
        assert isinstance(obsdate, dt.date), "obsdate needs to be of datetime.date type " \
                                             f"(currently {obsdate} is {type(obsdate)})."
        self._obsdate = obsdate

    @property
    def timerange(self) -> tuple:
        """Returns a tuple with the start and end time of the observation in datetime format.
        """
        return self._startime, self._endtime

    @timerange.setter
    def timerange(self, times: Tuple[dt.datetime]):
        """Start and end time of the observation in datetime format.
        Input:
            - times : tuple of datetime
                Tupple with (startime, endtime), each of them in datetime format.
        """
        starttime, endtime = times
        assert isinstance(starttime, dt.datetime)
        assert isinstance(endtime, dt.datetime)
        self._startime = starttime
        self._endtime = endtime
        self.epoch = starttime.date()

    def __init__(self, projectname: str, observatory: Optional[str] = '', params: Optional[dict] = None,
                 cwd: Optional[Union[Path, str]] = None):
        super().__init__(projectname, observatory, params, cwd)
        self._msfile = self.cwd / Path(f"{projectname.lower()}.ms")
        self._uvfitsfile = self.cwd / self.msfile.name.replace('.ms', '.uvfits')
        self._freqsetup = None
        self._sources = Sources()
        self._antennas = Antennas()
        self._refant = list()
        if self.msfile.exists():
            try:
                self.get_metadata_from_ms()
            except KeyError as e:
                self.listobs()
                rprint("\n\n[center][bold red] --- No source information is provided --- [/bold red][/center]\n")
                rprint("listobs() has run on the associated MS and you should now provide the inforamation "
                       "about the sources (which ones are targets, calibrators, etc.) in order to continue.\n"
                       "You can add such information in the input file.")
                raise KeyError(e)

    def initialize(self):
        # TODO: differenciate between observatories
        self.calibrate = calibration.Calibration(self)
        # if self.observatory.upper() in ('EVN', 'EEVN', 'E-EVN', 'VLBI', 'VLBA', 'LBA'):
        #     self.calibrate = vlbicalibration.Calibration(self)
        # elif self.observatory.upper() is 'GMRT':
        #     self.calibrate = gmrtcalibration.Calibration(self)

        # self.flag = flagging.Flagging(self)
        # self.plot = plotting.Plotting(self)
        # self.image = imaging.Imaging(self)
        pass

    def exists(self):
        """Returns if the associated MS file exists.
        """
        return self.msfile.exists()

    def get_metadata_from_ms(self):
        """Recovers all useful metadata from the MS file and stores it in the object.

        -- May Raise
        It may raise KeyError if the type of sources (target, phaseref, fringefinder) are not specified
        in the self.params dict.
        """
        try:
            with pt.table(self.msfile.name, readonly=True, ack=False) as ms:
                with pt.table(ms.getkeyword('ANTENNA'), readonly=True, ack=False) as ms_ant:
                    antenna_col = ms_ant.getcol('NAME')
                    for ant_name in antenna_col:
                        ant = Antenna(name=ant_name, observed=True)
                        self.antennas.add(ant)

                with pt.table(ms.getkeyword('DATA_DESCRIPTION'), readonly=True, ack=False) as ms_spws:
                    spw_names = ms_spws.getcol('SPECTRAL_WINDOW_ID')

                with pt.table(ms.getkeyword('OBSERVATION'), readonly=True, ack=False) as ms_obs:
                    telescope_name = ms_obs.getcol('TELESCOPE_NAME')
                    if self.observatory is None:
                        self.observatory = ','.join(telescope_name)
                    elif self.observatory != telescope_name[0]:
                        rprint("[orange]WARNING: the observatory name in MS does not match the one "
                               "provided in the project[/orange]")

                ant_subband = defaultdict(set)
                print('\nReading the MS to find which antennas actually observed...')
                with progress.Progress() as progress_bar:
                    task = progress_bar.add_task("[yellow]Reading MS...", total=len(ms))
                    for (start, nrow) in chunkert(0, len(ms), 5000):
                        ants1 = ms.getcol('ANTENNA1', startrow=start, nrow=nrow)
                        ants2 = ms.getcol('ANTENNA2', startrow=start, nrow=nrow)
                        spws = ms.getcol('DATA_DESC_ID', startrow=start, nrow=nrow)
                        msdata = ms.getcol('DATA', startrow=start, nrow=nrow)

                        for ant_i,antenna_name in enumerate(antenna_col):
                            for spw in spw_names:
                                cond = np.where((ants1 == ant_i) & (ants2 == ant_i) & (spws == spw))
                                if not (abs(msdata[cond]) < 1e-5).all():
                                    ant_subband[antenna_name].add(spw)

                        progress_bar.update(task, advance=nrow)

                for antenna_name in self.antennas.names:
                    self.antennas[antenna_name].subbands = tuple(ant_subband[antenna_name])
                    self.antennas[antenna_name].observed = len(self.antennas[antenna_name].subbands) > 0

                with pt.table(ms.getkeyword('FIELD'), readonly=True, ack=False) as ms_field:
                    src_names = ms_field.getcol('NAME')
                    src_coords = ms_field.getcol('PHASE_DIR')
                    for a_name, a_coord in zip(src_names, src_coords):
                        if a_name in self.params['target']:
                            a_type = SourceType.target
                        elif a_name in self.params['phaseref']:
                            a_type = SourceType.calibrator
                        elif a_name in self.params['fringefinder']:
                            a_type = SourceType.fringefinder
                        else:
                            a_type = SourceType.other

                        self.sources.append(Source(a_name, a_type, coord.SkyCoord(*a_coord[0], unit=(u.rad, u.rad))))

                with pt.table(ms.getkeyword('OBSERVATION'), readonly=True, ack=False) as ms_obs:
                    self.timerange = dt.datetime(1858, 11, 17, 0, 0, 2) + \
                         ms_obs.getcol('TIME_RANGE')[0]*dt.timedelta(seconds=1)

                with pt.table(ms.getkeyword('SPECTRAL_WINDOW'), readonly=True, ack=False) as ms_spw:
                    self.freqsetup = FreqSetup(ms_spw.getcol('NUM_CHAN')[0],
                                               ms_spw.getcol('CHAN_FREQ'),
                                               ms_spw.getcol('TOTAL_BANDWIDTH')[0])
        except RuntimeError:
            print(f"WARNING: {self.msfile} not found.")


    def import_lba_fits(self, fitsfile: str):
        # TODO: LBA import in CASA?
        raise NotImplementedError

    def import_evn_fitsidi(self, fitsidifiles: Union[list, str], ignore_antab=False, ignore_uvflg=False):
        """Imports the provided FITS-IDI files from an EVN observation into a single MS file.
        If checks if the FITS-IDI files already contain the Tsys and GC tables. Otherwise, it
        will first append such information and then import the files.

        If a .uvflg file exists (AIPS-format flag file), it will convert it to a CASA-compatible
        .flag file.

        Inputs
            fitsidifiles : list or str
                If str, it will retrieve all files that match the given path/name.
                If list, it must be an ordered list of all FITS-IDI files that will be imported.
            ignore_antab : False
                If the FITS-IDI files should be imported into a MS without caring about the ANTAB
                information (not recommended).
            ignore_uvflg : False
                If the FITS-IDI files should be imported into a MS without caring about the .uvflg
                file (AIPS-format a-priori flag table) information (not recommended).
        """
        # First run the Tsys and GC appending to the FITS-IDI. These functions will do nothing
        # if the information is already there.
        if isinstance(fitsidifiles, str):
            fitsidifiles = sorted(glob.glob(fitsidifiles), key=natsort_keygen())
        elif not isinstance(fitsidifiles, list):
            raise TypeError(f"The argument fitsidifiles should be either a list or str. "
                            f"Currently is {fitsidifiles} (type {type(fitsidifiles)})")

        for a_fitsidi in fitsidifiles:
            if not os.path.isfile(a_fitsidi):
                raise FileNotFoundError(f"The file {a_fitsidi} could not be found.")

        antabfile = self.cwd / Path(f"{self.projectname.lower()}.antab")
        uvflgfile = self.cwd / Path(f"{self.projectname.lower()}.uvflg")
        if not ignore_antab:
            assert antabfile.exists()
            fitsidi.append_tsys(str(antabfile), fitsidifiles)
            fitsidi.append_gc(str(antabfile), fitsidifiles[0])

        if not ignore_uvflg:
            assert uvflgfile.exists()
            fitsidi.convert_flags(infile=str(self.cwd / Path(f"{self.projectname.lower()}.uvflg")),
                                  idifiles=fitsidifiles,
                                  outfile=str(self.cwd / Path(f"{self.projectname.lower()}.flag")))

        casatasks.importfitsidi(vis=str(self.msfile), fitsidifile=fitsidifiles, constobsid=True,
                                scanreindexgap_s=8.0, specframe='GEO')
        self.get_metadata_from_ms()

    def listobs(self, listfile: Union[Path, str] = None):
        listfile = self.logdir / f"{self.projectname}-listobs.log" if listfile is None else listfile
        casatasks.listobs(self.msfile.name, listfile=str(listfile))
