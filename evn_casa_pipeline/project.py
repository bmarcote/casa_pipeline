#!/usr/bin/env python3
"""Defines a VLBI experiment with all the relevant metadata required during the post-processing.
The metadata is obtained from the MS itself.
"""
import os
import sys
import pickle
import glob
import json
import subprocess
import datetime as dt
import numpy as np
from pathlib import Path
from dataclasses import dataclass
from collections import defaultdict
from natsort import natsort_keygen
from pyrap import tables as pt
import casatasks
from enum import Enum
from astropy import units as u
from rich import print as rprint
# import blessed
from astropy import coordinates as coord
from casavlbitools import fitsidi
from . import calibration
from . import flagging
from . import imaging


def chunkert(f, l, cs):
    """Silly function to select a subset of an interval (f, f+cs) : 0 < f < l.
    """
    while f<l:
        n = min(cs, l-f)
        yield (f, n)
        f = f + n

percent = lambda x, y: (float(x)/float(y))*100.0


def cli_progress_bar(current_val, end_val, bar_length=40):
    """Creates a progress bar in the terminal.
    """
    percent = current_val/end_val
    hashes = '#'*int(round(percent*bar_length))
    spaces = ' '*(bar_length-len(hashes))
    sys.stdout.write("\rProgress: [{0}] {1}%".format(hashes+spaces, int(round(percent*100))))
    sys.stdout.flush()


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

    def __init__(self, name: str, sourcetype: SourceType, coordinates: coord.SkyCoord):
        assert isinstance(name, str), f"The name of the source must be a string (currrently {name})"
        assert isinstance(sourcetype, SourceType), \
               f"The name of the source must be a SourceType object (currrently {sourcetype})"
        assert isinstance(coordinates, coord.SkyCoord), "The coordinates of the source must " \
               f"be an astropy.coordinates.SkyCoord object (currently {coordinates})"
        self._name = name
        self._type = sourcetype
        self._coordinates = coordinates

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
    def __init__(self, sources=None):
        if sources is not None:
            self._sources = sources[:]
        else:
            self._sources = []

    def add(self, new_source):
        self.append(new_source)

    def append(self, new_source):
        assert isinstance(new_source, Source)
        self._sources.append(new_source)

    @property
    def names(self):
        """Returns the names of all sources included in Sources.
        """
        return [s.name for s in self._sources]

    @property
    def types(self):
        """Returns the types of all sources in Sources.
        """
        return [s.type for s in self._sources]

    @property
    def coordinates(self):
        """Returns the coordinates of all sources in Sources.
        """
        return [s.coordinates for s in self._sources]

    @property
    def targets(self):
        """Returns the target sources.
        """
        return Sources([s for s in self._sources if s.type == SourceType.target])

    @property
    def calibrators(self):
        """Returns the phase calibrator sources.
        """
        return Sources([s for s in self._sources if s.type == SourceType.calibrator])

    @property
    def other_sources(self):
        """Returns the sources catalogued as "others".
        """
        return Sources([s for s in self._sources if s.type == SourceType.other])

    @property
    def fringe_finders(self):
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
    """Defines an antenna in the array.
    """
    name: str
    observed: bool = False
    subbands: tuple = ()


class Antennas(object):
    """List of antennas (Antenna class)
    """
    def __init__(self, antennas=None):
        if antennas is not None:
            self._antennas = antennas[:]
        else:
            self._antennas = []

    def add(self, new_antenna: Antenna):
        """Adds a new antenna to the list of antennas.
        """
        assert isinstance(new_antenna, Antenna)
        self._antennas.append(new_antenna)

    @property
    def names(self):
        """Returns the name of the antenna
        """
        return [a.name for a in self._antennas]

    @property
    def observed(self):
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
    def frequency_range(self) -> tuple:
        """Returns the lowest and highest frequency observed.
        """
        return (self.frequency - self.bandwidth/2, self.frequency + self.bandwidth/2)

    def __init__(self, chans: int, freqs, bandwidth_spw):
        """Inputs:
            - chans : int
                Number of channels per subband.
            - freqs : array-like
                Reference frequency for each channel and subband (NxM array, M number
                of channels per subband.
            - bandwidth_spw : float or astropy.units.Quantity
                Total bandwidth for each subband. If not units are provided, Hz are assumed.
        """
        self._n_subbands = freqs.shape[0]
        assert isinstance(chans, (int, np.int32, np.int64)), \
            f"Chans {chans} is not an int as expected (found type {type(chans)})."
        assert isinstance(bandwidth_spw, float) or isinstance(bandwidth_spw, u.Quantity), \
            f"Bandwidth {bandwidth_spw} is not a float or Quantity as expected (found {type(bandwidth_spw)})."
        assert freqs.shape == (self._n_subbands, chans)
        self._channels = int(chans)
        self._freqs = np.copy(freqs)
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
    """Defines and EVN experiment with all relevant metadata.
    """
    @property
    def projectname(self) -> str:
        """Name of the EVN experiment, in upper case.
        """
        return self._expname

    @property
    def epoch(self) -> dt.date:
        """Epoch at which the EVN experiment was observed (starting date).
        """
        return self._obsdate

    @epoch.setter
    def epoch(self, obsdate: dt.date):
        assert isinstance(obsdate, dt.date), "obsdate needs to be of datetime.date type " \
                                             f"(currently {obsdate} is {type(obsdate)})."
        self._obsdate = obsdate

    @property
    def timerange(self):
        """Start and end time of the observation in datetime format.
        """
        return self._startime, self._endtime

    @timerange.setter
    def timerange(self, times):
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

    def __init__(self, projectname: str, params: dict):
        """Initializes an EVN experiment with the given name.

        Inputs:
        - expname : str
               The name of the experiment (case insensitive).
        """
        self._expname = projectname
        self._obsdate = None
        self._refant = []
        self._startime = None
        self._endtime = None
        # logpath = Path("./logs")
        # logpath.mkdir(parents=True, exist_ok=True)
        # self._logs = {'dir': logpath, 'file': Path("./processing.log)")}
        self._local_copy = Path(f"./{self._expname.lower()}.obj")
        self._last_step = None
        self._args = params

    # @property
    # def logfile(self):
    #     """Returns a dict with the logs, with two keys:
    #     - 'dir': the directory where individual log files can be stored (by default 'logs/')
    #     - 'file': the 'processing.log' file which stores all steps that run during the post-processing
    #               of the experiment.
    #     """
    #     return self._logs

    # def log(self, entry, timestamp=False):
    #     """Writes into the processing.log file a new entry.
    #     """
    #     if timestamp:
    #         cmd = f"# {dt.datetime.today().strftime('%d-%m-%Y %H:%M')}\n{entry}\n"
    #     else:
    #         cmd = f"{entry}\n"
    #
    #     with open(self.logfile['file'], 'a') as logfile:
    #         logfile.write(cmd)


    def exists_local_copy(self):
        """Checks if there is a local copy of the Experiment object stored in a local file.
        """
        return self._local_copy.exists()

    def store(self, path=None):
        """Stores the current Experiment into a file in the indicated path. If not provided,
        it will be '.{expname.lower()}.obj' where exp is the name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'wb') as file:
            pickle.dump(self, file)

    def store_json(self, path=None):
        """Stores the current Experiment into a JSON file.
        If path not prvided, it will be '{expname.lower()}.json'.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'wb') as file:
            json.dump(self.json(), file, cls=ExpJsonEncoder, indent=4)

    def load(self, path=None):
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

    # def __iter__(self):
    #     for key in ('expname', 'eEVNname', 'piname', 'email', 'supsci', 'obsdate', 'obsdatetime',
    #                 'timerange', 'sources', 'sources_stdplot', 'antennas', 'refant', 'credentials',
    #                 'cwd', 'logfile', 'vix', 'expsum', 'special_params',
    #                 'last_step', 'gui', 'correlator_passes'):
    #         yield key, getattr(self, key)

    def json(self):
        """Returns a dict with all attributes of the object.
        I define this method to use instead of .__dict__ as the later only reporst
        the internal variables (e.g. _username instead of username) and I want a better
        human-readable output.
        """
        d = dict()
        for key, val in self.__iter__():
            if hasattr(val, 'json'):
                d[key] = val.json()
            elif isinstance(val, Path):
                d[key] = val.name
            elif isinstance(val, dt.datetime):
                d[key] = val.strftime('%Y-%m-%d')
            elif isinstance(val, dt.date):
                d[key] = val.strftime('%Y-%m-%d')
            elif isinstance(val, list) and (len(val) > 0) and hasattr(val[0], 'json'):
                d[key] = [v.json() for v in val]
            elif isinstance(val, tuple) and (len(val) > 0) and isinstance(val[0], dt.datetime):
                d[key] = [v.strftime('%Y-%m-%d %H:%M:%S') for v in val]
            elif isinstance(val, dict):
                d[key] = {}
                for k, v in val:
                    if hasattr(v, 'json'):
                        d[key][k] = v.json()
                    elif hasattr(v, 'name'):
                        d[key][k] = v.name
                    else:
                        d[key][k] = v
            else:
                d[key] = val

        return d


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
    def msfile(self, new_msfile):
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
    def uvfitsfile(self, new_uvfitsfile):
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
    def inputs(self) -> dict:
        """Returns a dictionary with the inputs that the user set in the input file to run the pipeline.
        """
        return self._inputs

    def __init__(self, projectname: str, inputs: dict = {}):
        super().__init__(projectname, inputs)
        self._msfile = Path(f"{projectname.lower()}.ms")
        self._uvfitsfile = str(self.msfile).replace('.ms', '.uvfits') if '.ms' in str(self.msfile) \
                           else f"{str(self.msfile)}.uvfits"
        self._freqsetup = None
        self._sources = Sources()
        self._antennas = Antennas()
        self._refant = list()
        if self.msfile.exists():
            self.get_metadata_from_ms()

        self.calibrate = calibration.Calibration(self)
        self.flag = flagging.Flagging(self)
        self.images = imaging.Imaging(self)

    def exists(self):
        """Returns if the associated MS file exists.
        """
        return self.msfile.exists()

    def get_metadata_from_ms(self):
        """Recovers all useful metadata from the MS file and stores it in the object.
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

                ant_subband = defaultdict(set)
                print('\nReading the MS to find the missing antennas...')
                for (start, nrow) in chunkert(0, len(ms), 5000):
                    ants1 = ms.getcol('ANTENNA1', startrow=start, nrow=nrow)
                    ants2 = ms.getcol('ANTENNA2', startrow=start, nrow=nrow)
                    spws = ms.getcol('DATA_DESC_ID', startrow=start, nrow=nrow)
                    msdata = ms.getcol('DATA', startrow=start, nrow=nrow)
                    cli_progress_bar(start, len(ms), bar_length=40)

                    for ant_i,antenna_name in enumerate(antenna_col):
                        for spw in spw_names:
                            cond = np.where((ants1 == ant_i) & (ants2 == ant_i) & (spws == spw))
                            if not (abs(msdata[cond]) < 1e-5).all():
                                ant_subband[antenna_name].add(spw)

                for antenna_name in self.antennas.names:
                    self.antennas[antenna_name].subbands = tuple(ant_subband[antenna_name])
                    self.antennas[antenna_name].observed = len(self.antennas[antenna_name].subbands) > 0

                with pt.table(ms.getkeyword('FIELD'), readonly=True, ack=False) as ms_field:
                    src_names = ms_field.getcol('NAME')
                    src_coords = ms_field.getcol('PHASE_DIR')
                    for a_name, a_coord in zip(src_names, src_coords):
                        #TODO: check how I do specify this in the json file
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
        raise NotImplementedError

    def import_evn_fitsidi(self, fitsidifiles):
        """Imports the provided FITS-IDI files from an EVN observation into a single MS file.
        If checks if the FITS-IDI files already contain the Tsys and GC tables. Otherwise, it
        will first append such information and then import the files.

        Inputs
            fitsidifiles : list or str
                If str, it will retrieve all files that match the given path/name.
                If list, it must be an ordered list of all FITS-IDI files that will be imported.
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

        fitsidi.append_tsys(f"{self.projectname.lower()}.antab", fitsidifiles)
        fitsidi.append_gc(f"{self.projectname.lower()}.antab", fitsidifiles[0])
        casatasks.importfitsidi(vis=str(self.msfile), fitsidifile=fitsidifiles, constobsid=True,
                                scanreindexgap_s=8.0, specframe='GEO')
        self.get_metadata_from_ms()

    def __iter__(self):
        for key in ('msfile', 'fitsidifiles', 'sources', 'antennas', 'freqsetup'):
            yield key, getattr(self, key)

    def json(self):
        """Returns a dict with all attributes of the object.
        I define this method to use instead of .__dict__ as the later only reporst
        the internal variables (e.g. _username instead of username) and I want a better
        human-readable output.
        """
        d = dict()
        for key, val in self.__iter__():
            if hasattr(val, 'json'):
                d[key] = val.json()
            elif isinstance(val, Path):
                d[key] = val.name
            elif isinstance(val, list):
                d[key] = [v.json() for v in val]
            else:
                d[key] = val

        return d



class ExpJsonEncoder(json.JSONEncoder):
    """Encodes properly the Experiment class to be able to be written as a JSON format
    """
    def default(self, obj):
        if isinstance(obj, dt.datetime):
            return obj.strftime('%Y-%m-%d %H:%M:%S')
        elif isinstance(obj, dt.date):
            return obj.strftime('%Y-%m-%d')
        elif hasattr(obj, 'json'):
            return obj.json()
        elif isinstance(obj, Path):
            return obj.name
        elif isinstance(obj, np.ndarray):
            return list(obj)
        else:
            print(obj)
            return json.JSONEncoder.default(self, obj)


