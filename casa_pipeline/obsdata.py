#!/usr/bin/env python3
"""Defines all clases related to a VLBI experiment with all the relevant metadata.
The metadata is obtained from the MS itself.
"""
from __future__ import annotations
import os
import glob
import shutil
import logging
import datetime as dt
import numpy as np
from enum import IntEnum
from pathlib import Path
from concurrent import futures
from threading import Lock
from typing import Optional, Iterable, NoReturn, List, Union, Tuple
from dataclasses import dataclass
from collections import defaultdict
from natsort import natsort_keygen
import casatasks
from casatools import msmetadata as msmd
from casatools import table as tb
from enum import Enum
from astropy import units as u
from rich import print as rprint
from rich import progress
# import blessed
from astropy import coordinates as coord
from casa_pipeline.casavlbitools import fitsidi
from casa_pipeline.casa_pipeline import tools
# from .calibration import Calibration
# from . import flagging
# from . import imaging
# from . import plotting


# Because as far as I know CASA folks do not have a straight away to get this simple information
class Stokes(IntEnum):
    """The Stokes types defined as in the enum class from casacore code.
    """
    Undefined = 0 # Undefined value
    I = 1 # standard stokes parameters
    Q = 2
    U = 3
    V = 4
    RR = 5 # circular correlation products
    RL = 6
    LR = 7
    LL = 8
    XX = 9 # linear correlation products
    XY = 10
    YX = 11
    YY = 12
    RX = 13 # mixed correlation products
    RY = 14
    LX = 15
    LY = 16
    XR = 17
    XL = 18
    YR = 19
    YL = 20
    PP = 21 # general quasi-orthogonal correlation products
    PQ = 22
    QP = 23
    QQ = 24
    RCircular = 25 # single dish polarization types
    LCircular = 26
    Linear = 27
    Ptotal = 28 # Polarized intensity ((Q^2+U^2+V^2)^(1/2))
    Plinear = 29 #  Linearly Polarized intensity ((Q^2+U^2)^(1/2))
    PFtotal = 30 # Polarization Fraction (Ptotal/I)
    PFlinear = 31 # linear Polarization Fraction (Plinear/I)
    Pangle = 32 # linear polarization angle (0.5  arctan(U/Q)) (in radians)


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

    @type.setter
    def type(self, new_type: Union[str, SourceType]):
        """Sets the type of the given source.

        new_type : [str/SourceType]
            The new type for the given source. If string, it needs to match
            'target', 'phaseref', 'fringefinder' or 'other'.
            Otherwise, it is expected to be a 'SourceType' parameter.
        """
        if isinstance(new_type, SourceType):
            self._type =new_type
        elif isinstance(new_type, str):
            if new_type == 'target':
                self._type = SourceType.target
            elif new_type == 'phaseref':
                self._type = SourceType.calibrator
            elif new_type == 'fringefinder':
                self._type = SourceType.fringefinder
            elif new_type == 'other':
                self._type = SourceType.other
            else:
                raise ValueError(f"The new type given for the source {self.name} " \
                                 f"({new_type}) has an unexpected value, " \
                                 f"while 'target'/'phaseref'/'fringefinder'/'other' were expected.")
        else:
            raise ValueError(f"The new type given for the source {self.name} " \
                             f"({new_type}) has an unexpected type {type(new_type)}, " \
                             "while str/SourceType were expected.")

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
    def phase_calibrators(self) -> Sources:
        """Returns the phase calibrator sources.
        """
        return Sources([s for s in self._sources if s.type == SourceType.calibrator])

    @property
    def all_calibrators(self) -> Sources:
        """Returns all calibrator sources (including phase calibrators and fringe finders).
        """
        return Sources([s for s in self._sources if (s.type == SourceType.calibrator) or \
                       (s.type == SourceType.fringefinder)])

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
        return self.__repr__()

    def __repr__(self):
        return "Antennas:\n" \
           f"([{', '.join([ant if ant in self.observed else '['+ant+']' for ant in self.names])})"

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
    def frequency(self) -> u.Quantity:
        """Returns the central frequency of the observation, in a astropy.units.Quantity format.
        """
        return self._frequency.to(u.GHz)

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

    def __init__(self, channels: int, n_subbands: int, central_frequency: Union[float, u.Quantity],
                 bandwidth_spw: Union[float, u.Quantity]):
        """Inputs:
            - channels : int
                Number of channels per subband.
            - n_subbands : int
                Number of subbands.
            - central_frequency: float or astropy.units.Quantity
                Central frequency of the observation. If no units are provided, Hz are assumed.
            - bandwidth_spw : float or astropy.units.Quantity
                Total bandwidth for each subband. If no units are provided, Hz are assumed.
        """
        assert isinstance(n_subbands, (int, np.int32, np.int64)) and n_subbands > 0, \
            f"n_subbands {n_subbands} is not a positive integer as expected " \
            f"(found type {type(n_subbands)})."
        assert isinstance(channels, (int, np.int32, np.int64)), \
            f"Chans {channels} is not an integer as expected (found type {type(channels)})."
        assert isinstance(bandwidth_spw, float) or isinstance(bandwidth_spw, u.Quantity), \
            f"Bandwidth {bandwidth_spw} is not a float or Quantity as expected " \
            f"(found {type(bandwidth_spw)})."
        self._channels = int(channels)
        self._n_subbands = n_subbands
        if isinstance(bandwidth_spw, float):
            self._frequency = central_frequency*u.Hz
        else:
            self._frequency = central_frequency
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
        return f"{self.frequency:.3f} central frequency.\n{self.n_subbands} x " \
               f"{self.bandwidth_per_subband.to(u.MHz):.1f} " \
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


class Ms(object):
    """Defines one MS object for a given observation.
    It contains all relevant information that is MS-depended, e.g. associated MS file,
    frequency setup, sources, etc. But also objects associated to the calibration tasks that
    can run in the MS.
    """
    @property
    def prefixname(self) -> str:
        """Name of the project associated to the observation.
        """
        return self._prefixname

    @property
    def observatory(self) -> str:
        """Name of the observatory that conducted the observation.
        """
        return self._observatory_name

    @property
    def cwd(self) -> Path:
        """Returns the working directory for the project.
        """
        return self._cwd

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
    def refant(self, new_refant: Union[str, list]):
        """Defines the antenna to use as reference during calibration.
        It can be either a str (with a single antenna, or comma-separated antennas),
        or a list with the antenna(s) to use as reference.
        """
        if isinstance(new_refant, str):
            new_refant = [ant.strip() for ant in new_refant.split(',')]
        elif not isinstance(new_refant, list):
            raise TypeError("The antenna introduced must be either a  str or list")

        # As it may be set before the antenna information was retrieved from the MS
        if len(self.antennas) > 0:
            for antenna in new_refant:
                assert antenna in self.antennas.names, \
                       f"Antenna {antenna} for refant is not in the array."

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
    def params(self) -> dict:
        return self._params

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

    @property
    def splits(self):
        return self._splits

    @property
    def logger(self) -> logging.Logger:
        return self._logger

    def __init__(self, prefixname: str, observatory: Optional[str] = None,
                 cwd: Optional[Union[Path, str]] = None, params: dict = None,
                 logger: Optional[logging.Logger] = None):
        assert (prefixname != '') and (isinstance(prefixname, str)), \
               "The prefix name needs to be a non-empty string."
        self._prefixname = prefixname
        self._observatory_name = observatory if observatory is not None else ''
        if cwd is None:
            self._cwd = Path('./')
        elif isinstance(cwd, str):
            self._cwd = Path(cwd)
        elif isinstance(cwd, Path):
            self._cwd = cwd
        else:
            TypeError(f"The working directory ({cwd}) should be either None, a Path type, or str.")

        self._msfile = self.cwd / Path(f"{prefixname}.ms")
        self._uvfitsfile = self.cwd / self.msfile.name.replace('.ms', '.uvfits')
        self._freqsetup = None
        self._sources = Sources()
        self._antennas = Antennas()
        self._splits = defaultdict(list)
        self._params = params if params is not None else {}
        if logger is None:
            self._logger = logging.getLogger(f"Ms-{self.prefixname}")
        else:
            self._logger = logger

        if 'reference_antenna' in params:
            self.refant = params['reference_antenna']

        if self.msfile.exists():
            self.get_metadata_from_ms()
            self.listobs()

        self._logger.debug(f"MS object created with {self.__dict__}.")


    def exists(self):
        """Returns if the associated MS file exists.
        """
        return self.msfile.exists()

    def get_metadata_from_ms(self):
        """Recovers all useful metadata from the MS file and stores it in the object.

        -- May Raise
        It may raise KeyError if the type of sources (target, phaseref, fringefinder)
        are not specified in the self.params dict.
        """
        m = msmd(str(self.msfile))
        if not m.open(str(self.msfile)):
            return ValueError(f"The MS file {self.msfile} could not be openned.")

        try:
            antenna_names = m.antennanames()

            for ant_name in antenna_names:
                ant = Antenna(name=ant_name, observed=True)
                self.antennas.add(ant)

            spw_names = range(m.nspw())
            telescope_name = m.observatorynames()[0]
            if self.observatory == '':
                self._observatory_name = telescope_name
            elif self.observatory != telescope_name:
                rprint("[yellow]WARNING: the observatory name in MS does not match the one "
                       "provided in the project[/yellow]")

            src_names = m.fieldnames()
            src_coords = [m.phasecenter(s) for s in range(len(src_names)) ]
            for a_name, a_coord in zip(src_names, src_coords):
                try:
                    if a_name in self.params['target']:
                        a_type = SourceType.target
                    elif a_name in self.params['phaseref']:
                        a_type = SourceType.calibrator
                    elif a_name in self.params['fringefinder']:
                        a_type = SourceType.fringefinder
                    else:
                        a_type = SourceType.other
                except KeyError:
                    rprint("[bold yellow]-- No source type information has been " \
                           "found --[bold yellow]")
                    rprint("You better define manually or through the inputs which sources " \
                           "are target/calibrator/etc.")
                    a_type = SourceType.other

                self.sources.append(Source(a_name, a_type,
                                           coord.SkyCoord(ra=a_coord['m0']['value'],
                                                          dec=a_coord['m1']['value'],
                                                          unit=(a_coord['m0']['unit'],
                                                                a_coord['m1']['unit']),
                                                          equinox=a_coord['refer'])))

            timerange = m.timerangeforobs(0)
            self.timerange = [tools.mjd2date(timerange['begin']['m0']['value']),
                              tools.mjd2date(timerange['end']['m0']['value'])]

            mean_freq = (m.meanfreq(spw_names[-1]) + m.meanfreq(0)) / 2.0
            self.freqsetup = FreqSetup(m.nchan(0), m.nspw(), mean_freq,
                                       m.bandwidths()[0])

            nrows = int(m.nrows())

            # To be able  to get the parallel hands, either circular or linear
            corr_order = [Stokes(i) for i in m.corrtypesforpol(0)]
            corr_pos = []
            try:
                corr_pos.append(corr_order.index(Stokes.RR))
                corr_pos.append(corr_order.index(Stokes.LL))
            except ValueError:
                try:
                    corr_pos.append(corr_order.index(Stokes.XX))
                    corr_pos.append(corr_order.index(Stokes.YY))
                except ValueError:
                    rprint("[bold red]The associated MS does not have neither circular nor " \
                           "linear-based polarization information[/bold red]")
        finally:
            m.close()

        m = tb(str(self.msfile))
        try:
            if not m.open(str(self.msfile)):
                return ValueError(f"The MS file {self.msfile} could not be openned.")

            ant_subband = defaultdict(set)
            # TODO: parallelize the following. This approach blocks indefinitely. Because of the MS.
            # if test:
            #     ant_subband_mutex = Lock()
            #     rprint('\n[bold]Reading the MS to find which antennas actually observed.[/bold]')
            #     with futures.ProcessPoolExecutor(max_workers=8) as executor:
            #         workers = []
            #         with progress.Progress() as progress_bar:
            #             task = progress_bar.add_task("[yellow]Reading MS...", total=nrows)
            #             for (start, nrow) in tools.chunkert(0, nrows, 5000):
            #                 kwargs = {'mset': m, 'ant_subband_dict': ant_subband, 'startrow': start,
            #                           'nrow': nrow, 'antenna_names': antenna_names,
            #                           'spw_names': spw_names, 'corr_pos': corr_pos,
            #                           'mutex': ant_subband_mutex}
            #                 print(f"Launching {start}")
            #                 workers.append(executor.submit(self._get_spw_per_ant, **kwargs))
            #                 progress_bar.update(task, advance=nrow)
            # else:
            rprint('\n[bold]Reading the MS to find which antennas actually observed...[/bold]')
            with progress.Progress() as progress_bar:
                task = progress_bar.add_task("[yellow]Reading MS...", total=nrows)
                for (start, nrow) in tools.chunkert(0, nrows, 5000):
                    ants1 = m.getcol('ANTENNA1', startrow=start, nrow=nrow)
                    ants2 = m.getcol('ANTENNA2', startrow=start, nrow=nrow)
                    spws = m.getcol('DATA_DESC_ID', startrow=start, nrow=nrow)
                    msdata = m.getcol('DATA', startrow=start, nrow=nrow)

                    for ant_i,antenna_name in enumerate(antenna_names):
                        for spw in spw_names:
                            cond = np.where((ants1 == ant_i) & (ants2 == ant_i) & (spws == spw))
                            if not (abs(msdata[corr_pos][:, :, cond[0]]) < 1e-5).all():
                                ant_subband[antenna_name].add(spw)

                    progress_bar.update(task, advance=nrow)
        finally:
            m.close()

        for antenna_name in self.antennas.names:
            self.antennas[antenna_name].subbands = tuple(ant_subband[antenna_name])
            self.antennas[antenna_name].observed = len(self.antennas[antenna_name].subbands) > 0

        self.listobs()


    def _get_spw_per_ant(self, mset, ant_subband_dict, startrow, nrow, antenna_names, spw_names,
                         corr_pos, mutex):
        ants1 = mset.getcol('ANTENNA1', startrow=startrow, nrow=nrow)
        ants2 = mset.getcol('ANTENNA2', startrow=startrow, nrow=nrow)
        spws = mset.getcol('DATA_DESC_ID', startrow=startrow, nrow=nrow)
        msdata = mset.getcol('DATA', startrow=startrow, nrow=nrow)
        print(r"Got data from starrow: {startrow}")

        for ant_i,antenna_name in enumerate(antenna_names):
            for spw in spw_names:
                cond = np.where((ants1 == ant_i) & (ants2 == ant_i) & (spws == spw))
                if not (abs(msdata[corr_pos][:, :, cond[0]]) < 1e-5).all():
                    print(f"Acquiring mutex for {startrow}")
                    mutex.acquire()
                    print(f"Acquired mutex for {startrow}")
                    ant_subband_dict[antenna_name].add(spw)
                    mutex.release()
                    print(f"Released mutex for {startrow}")


    def antennas_in_scan(self, scanno: int) -> list:
        """Returns a list with all antennas that observed the given scan.
        """
        m = msmd(str(self.msfile))
        try:
            if not m.open(str(self.msfile)):
                return ValueError(f"The MS file {self.msfile} could not be openned.")

            antenna_names = m.antennanames()
            ant_ids = m.antennasforscan(scanno)
        finally:
            m.close()

        return [antenna_names[a] for a in ant_ids]


    def scans_with_source(self, source_name: str) -> list:
        """Returns a list with the numbers of the scans where the given source was observed.
        """
        m = msmd(str(self.msfile))
        if not m.open(str(self.msfile)):
            return ValueError(f"The MS file {self.msfile} could not be openned.")

        scans = msmd.scansforfield(source_name)
        m.close()
        return scans


    def listobs(self, listfile: Union[Path, str] = None, overwrite=True) -> dict:
        listfile = self.cwd / f"{self.prefixname}-listobs.log" if listfile is None else listfile
        return casatasks.listobs(self.msfile.name, listfile=str(listfile), overwrite=overwrite)


    def split(self, sources: Optional[Union[Iterable[str], str, None]] = None,
              keepflags=False, chanbin: int =-1, **kwargs):
        """Splits all the data from all calibrated sources.
        If sources is None, then all sources will be split.

        Params:
            - chanbin : int  (default = -1)
              Width (bin) of input channels to average to for an output channel.
              If -1, then if will create a single-output channel per spw.

        It returns a dict with all split source names as keys. The values are the new Ms objects.
        """
        splits = {}

        if isinstance(sources, str):
            sources = [source.strip() for source in sources.split(',')]

        if sources is not None:
            for source in sources:
                assert source in self.sources, \
                       f"The passed source {source} is not in the MS {self.msfile}."

        if chanbin == -1 or chanbin > 1:
            kwargs['chanaverage'] = True
            kwargs['chanbin'] = self.freqsetup.channels

        for a_source in self.sources.names if sources is None else sources:
            np.int = int
            np.float = float  # these two is because CASA uses deprecated types in newer numpy ver!
            casatasks.mstransform(vis=str(self.msfile),
                                  outputvis=f"{self.prefixname}.{a_source}.ms",
                                  field=a_source, keepflags=keepflags, **kwargs)
            splits[a_source] = Ms(f"{self.prefixname}.{a_source}", cwd=self.cwd,
                                  params=self._params, logger=self._logger)
            self.splits[a_source].append(splits[a_source])
        return splits


    def export_uvfits(self, outfitsfilename: Union[str, None] = None, overwrite=True,
                      datacolumn='corrected'):
        """Export the available MS file into a UVFITS file, per source.
        """
        if outfitsfilename is None:
            outfitsfilename = str(self.msfile).replace(".ms", ".uvfits")

        casatasks.exportuvfits(vis=str(self.msfile), fitsfile=outfitsfilename,
                               datacolumn=datacolumn, multisource=(len(self.sources) > 1),
                               combinespw=True, padwithflags=True, overwrite=overwrite)




class Importing(object):
    """Class that contains different static functions to import data from different observatories
    into a MS. Each function would conduct the necessary steps to get a properly prepared MS file.
    """

    def __init__(self, ms: Ms):
        self._ms = ms

    def lba_fits(self, fitsfile: str):
        # TODO: LBA import in CASA?
        raise NotImplementedError

    def evn_download(self, projectcode: str = None, username: str = None, password: str = None):
        """Downloads the data files associated with the given EVN project.
        It will download the associated FITS-IDI files, and the associated .antab and .uvflg files.

        Inputs
            projectcode : str  (default = None)
                Project code of the
        """
        # TODO: Do I need the observing epoch?
        pass


    def evn_fitsidi(self, fitsidifiles: Union[list, str], ignore_antab=False,
                           ignore_uvflg=False, delete=False, replace_tsys=False):
        """Imports the provided FITS-IDI files from an EVN observation into a single MS file.
        If checks if the FITS-IDI files already contain the Tsys and GC tables. Otherwise, it
        will first append such information and then import the files.

        If a .uvflg file exists (AIPS-format flag file), it will convert it to a CASA-compatible
        .flag file.

        Inputs
            fitsidifiles : list or str
                If str, it will retrieve all files that match the given path/name.
                If list, it must be an ordered list of all FITS-IDI files that will be imported.
            ignore_antab : bool (default False)
                If the FITS-IDI files should be imported into a MS without caring about the ANTAB
                information (not recommended).
            ignore_uvflg : bool (default False)
                If the FITS-IDI files should be imported into a MS without caring about the .uvflg
                file (AIPS-format a-priori flag table) information (not recommended).
            delete : bool (default False)
                Overwrites the MS if existing.
            replace_tsys : bool (default False)
                Overwrites the Tsys/GC information, if exists, in the FITS-IDI files.
        """
        if self._ms.msfile.exists():
            if delete:
                shutil.rmtree(self._ms.msfile)
            else:
                rprint(f"[bold red]The MS file {self._ms.msfile} already exists "
                       "and you told me to not remove it[/bold red]")
                raise FileExistsError

        # First run the Tsys and GC appending to the FITS-IDI. These functions will do nothing
        # if the information is already there.
        if isinstance(fitsidifiles, str):
            fitsidifiles = sorted(glob.glob(fitsidifiles), key=natsort_keygen())
        elif isinstance(fitsidifiles, list):
            fitsidifiles = sorted(fitsidifiles, key=natsort_keygen())
        else:
            raise TypeError("The argument fitsidifiles should be either a list or str. "
                            f"Currently is {fitsidifiles} (type {type(fitsidifiles)})")

        self._ms.logger.debug(f"FITS-IDI files to read: {', '.join(fitsidifiles)}.")
        if len(fitsidifiles) == 0:
                raise FileNotFoundError(f"No files matching {fitsidifiles} could be found.")

        for a_fitsidi in fitsidifiles:
            if not os.path.isfile(a_fitsidi):
                raise FileNotFoundError(f"The file {a_fitsidi} could not be found.")

        antabfile = self._ms.cwd / Path(f"{self._ms.prefixname.lower()}.antab")
        uvflgfile = self._ms.cwd / Path(f"{self._ms.prefixname.lower()}.uvflg")
        flagfile = self._ms.cwd / Path(f"{self._ms.prefixname.lower()}.flag")
        if not ignore_antab:
            assert antabfile.exists(), \
                   f"The associated file {antabfile} should exist but was not found."
            rprint("[bold]Appending the Tsys values to the FITS-IDI files[/bold]")
            try:
                fitsidi.append_tsys(str(antabfile), fitsidifiles, replace=replace_tsys)
                rprint("[bold]Appending the GC values to the FITS-IDI files[/bold]")
                fitsidi.append_gc(str(antabfile), fitsidifiles[0], replace=replace_tsys)
            except RuntimeError as e:
                rprint(f"[yellow]{e}[/yellow]")

        if (not ignore_uvflg) and (not flagfile.exists()):
            assert uvflgfile.exists(), \
                   f"The associated file {uvflgfile} should exist but was not found."
            rprint("[bold]Converting the AIPS-style flags into CASA-style flags[/bold]")
            rprint(f"[green]The file {self._ms.prefixname.lower()}.flag containing the "
                   "a-priori flagging has been created.[/green]")
            fitsidi.convert_flags(infile=str(uvflgfile), idifiles=fitsidifiles,
                                  outfile=str(flagfile))

        rprint("[bold]Importing the FITS-IDI files into a MS[/bold]")
        casatasks.importfitsidi(vis=str(self._ms.msfile), fitsidifile=fitsidifiles, constobsid=True,
                                scanreindexgap_s=8.0, specframe='GEO')
        rprint(f"[green]The file {self._ms.msfile} has been created[/green]")
        self._ms.get_metadata_from_ms()

