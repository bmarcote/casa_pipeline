import typing
import pickle
from pathlib import Path
from enum import Enum
import numpy as np
from astropy import units as u
from astropy import coordinates as coord
import casatools
import casatasks


class SourceType(Enum):
    target = 0
    calibrator = 1
    bandpass = 2
    other = 3


class Source(object):
    """Defines a source that has been observed.
    """
    @property
    def name(self) -> str:
        """All sources have a name.
        """
        return self._name

    @property
    def type(self) -> SourceType:
        return self._type

    @property
    def coordinates(self) -> coord.SkyCoord:
        return self._coordinates

    @property
    def source_image(self) -> Path:
        """Returns the path to the latest produced image of the source.
        """
        return self._source_image

    def __init__(self, name: str, coordinates: coord.SkyCoord, source_type: SourceType=None):
        assert isinstance(name, str), f"The name of the source must be a string (currrently {name})."
        assert isinstance(coordinates, coord.SkyCoord), \
               f"coordinates must be an astropy.coordinates.SkyCoord object (currrently {type(coordinates)})"
        if source_type is None:
            source_type = SourceType.other
        assert isinstance(source_type, SourceType), \
               f"The name of the source must be a SourceType object (currrently {type(source_type)})."
        self._name = name
        self._type = source_type
        self._coordinates = coordinates
        self._source_image = None


class Project(object):

    @property
    def project_name(self) -> str:
        """Name of the project.
        It is used to search for and write all corresponding files during the data reduction.
        """
        return self._name

    def is_source_in_project(self, a_source: str) -> bool:
        """Checks that the given source has been observed in this project (e.g. is in the MS source list).
        """
        return a_source in self.sources

    @property
    def targets(self) -> list:
        """List of target sources.
        The order in the list should match the order of the phase_calibrators associated (if any).
        """
        return self._targets

    @targets.setter
    def targets(self, targets: list):
        assert isinstance(targets, list), "Targets should be a list of sources."
        for target in targets:
            assert self.is_source_in_project(target), \
                   f"The target {target} was not observed. " \
                   f"The available sources are {', '.join(list(self.sources.keys()))}."

        self._targets = targets

    @property
    def phase_calibrators(self) -> list:
        """List of phase calibrators sources.
        The order in the list should match the order of the targets associated (if any).
        """
        return self._phase_calibrators

    @phase_calibrators.setter
    def phase_calibrators(self, phase_calibrators: list):
        assert isinstance(phase_calibrators, list) or (phase_calibrators is None)
        for a_calibrator in phase_calibrators:
            assert self.is_source_in_project(a_calibrator), \
                   f"The phase calibrator {a_calibrator} was not observed. " \
                   f"The available sources are {', '.join(list(self.sources.keys()))}."

        self._phase_calibrators = phase_calibrators

    @property
    def bandpass_calibrators(self) -> list:
        """List of fringe finder or bandpass calibrator sources.
        """
        return self._bandpass_calibrators

    @bandpass_calibrators.setter
    def bandpass_calibrators(self, bandpass_calibrators: list):
        assert isinstance(bandpass_calibrators, list)
        for a_calibrator in bandpass_calibrators:
            assert self.is_source_in_project(a_calibrator), \
                   f"The bandpass calibrator {a_calibrator} was not observed. " \
                   f"The available sources are {', '.join(list(self.sources.keys()))}."

        self._bandpass_calibrators = bandpass_calibrators

    @property
    def calibrators(self) -> list:
        if self.phase_calibrators is not None:
            return self.bandpass_calibrators + self.targets
        else:
            return self.bandpass_calibrators + self.phase_calibrators

    @property
    def sources(self) -> dict:
        """Dictionary with all sources in the project. The keys are the source names
        and the values are Source-type objects for each source.
        """
        if self._sources is None:
            self.get_ms_info()

        return self._sources

    @property
    def antennas(self) -> list:
        return self._antennas

    @property
    def refants(self) -> list:
        return self._refants

    @refants.setter
    def refants(self, refants: list):
        if self.msfile.exists():
            for refant in refants:
                assert refant in self.antennas, f"The ref. antenna {refant} not in the array for {self.project_name}. " \
                                                f"The available antennas are {', '.join(self.antennas)}."
        self._refants = refants


    @property
    def idi_files(self) -> list:
        """Returns a list with all FITS IDI files associated to the project.
        """
        return self._idi_files

    @idi_files.setter
    def idi_files(self, files: list):
        """Stores the FITS IDI files associated with the project.
        """
        assert isinstance(files, list)
        self._idi_files = files

    @property
    def msfile(self) -> Path:
        """Returns the Path to the MS file associated to the project.
        """
        return Path(f"{self.project_name}.ms")

    @property
    def logdir(self) -> Path:
        """Returns the Path directory where all log files should be stored.
        """
        return self._logdir

    @logdir.setter
    def logdir(self, path: Path):
        # safety check
        if isinstance(path, str):
            self._logdir = Path(path)
        else:
            assert isinstance(path, Path)
            self._logdir = path

        if not self._logdir.exists():
            Path.mkdir(self._logdir, exist_ok=True)
            print(f"The directory {self._logdir.name} has been created.")

    @property
    def caldir(self) -> Path:
        """Returns the Path directory where all calibration tables should be stored.
        """
        return self._caldir

    @caldir.setter
    def caldir(self, path: Path):
        # safety check
        if isinstance(path, str):
            self._caldir = Path(path)
        else:
            assert isinstance(path, Path)
            self._caldir = path

        if not self._caldir.exists():
            Path.mkdir(self._caldir, exist_ok=True)
            print(f"The directory {self._caldir.name} has been created.")

    @property
    def nchannels(self) -> int:
        """Returns the number of channels per subband.
        """
        if self._nchannels is None:
            self.get_ms_info()

        return self._nchannels

    @property
    def nsubbands(self) -> int:
        """Returns the number of subbands.
        """
        if self._nsubbands is None:
            self.get_ms_info()

        return self._nsubbands

    @property
    def bandwidth(self) -> u.Quantity:
        """Returns the total bandwith of the observation.
        """
        if self._bandwidth is None:
            self.get_ms_info()

        return self._bandwidth

    @property
    def frequency(self) -> u.Quantity:
        """Returns the central frequency of the observation.
        """
        if self._frequency is None:
            self.get_ms_info()

        return self._frequency

    @property
    def instr_delay_timerange(self) -> str:
        """Returns the timerange defined to run the instrumental delay calibration.
        """
        return self._sbd_timer

    @instr_delay_timerange.setter
    def instr_delay_timerange(self, new_timerange: str):
        self._sbd_timer = new_timerange

    @property
    def last_step(self) -> str:
        """Returns the last step conducted successfully in the data reduction, from a potential
        previous run. If this is the first time, it will return None.
        """
        return self._last_step

    @last_step.setter
    def last_step(self, the_last_step: str):
        self._last_step = the_last_step

    @property
    def local_copy(self) -> Path:
        """Returns the Path to the (binary) file where the project Python data is stored.
        """
        return self._local_copy

    def __init__(self, name: str, bandpass_calibrators: list, targets: list, reference_antennas: list,
                 phase_calibrators: typing.Optional[list]=None, logdir: str='./log',
                 caldir: str='./calibration_tables', instr_delay_timerange: str=None):
        """Creates a Project object.
        Inputs:
            name : str
                Name associated to the project.
                It should match the preffix from all input files and it will be used for the name of
                all created files.
        """
        self._name = name
        self._idi_files = None
        self.logdir = Path(logdir)
        self.caldir = Path(caldir)
        self._sources = None
        self._nchannels = None
        self._antennas = None
        self._targets = targets
        self._bandpass_calibrators = bandpass_calibrators
        self._phase_calibrators = phase_calibrators
        self._refants = reference_antennas
        self._sbd_timer = instr_delay_timerange
        self._last_step = None
        self._local_copy = self.logdir / f"{self._name.lower()}.pyobj"


    def get_ms_info(self):
        """Recovers all the metadata from the MS
        """
        tb = casatools.table()
        # ms = casatools.ms()
        tb.open(self.msfile)
        # getting the required keywords, otherwise I need to open the MS again
        tb_ant = tb.getkeyword('ANTENNA')
        tb_spw = tb.getkeyword('SPECTRAL_WINDOW')
        # tb_pol = tb.getkeyword('POLARIZATION')
        tb_src = tb.getkeyword('FIELD')
        tb.open(tb_ant.replace('Table: ', ''))
        self._antennas = list(tb.getcol('NAME'))
        tb.open(tb_spw.replace('Table: ', ''))
        self._nsubbands = len(tb.getcol('TOTAL_BANDWIDTH'))
        self._bandwidth = np.sum(tb.getcol('TOTAL_BANDWIDTH'))*u.Hz
        self._nchannels = int(tb.getcol('NUM_CHAN')[0])
        self._frequency = np.mean(tb.getcol('CHAN_FREQ').reshape(self.nchannels*self.nsubbands))*u.Hz
        tb.open(tb_src.replace('Table: ', ''))
        tmp_src = {}
        src_coords = tb.getcol('PHASE_DIR')
        for i, src_name in enumerate(tb.getcol('NAME')):
            if src_name in self.targets:
                src_type = SourceType.target
            elif src_name in self.phase_calibrators:
                src_type = SourceType.calibrator
            elif src_name in self.bandpass_calibrators:
                src_type = SourceType.bandpass
            else:
                src_type = SourceType.other

            tmp_src[src_name] = Source(name=src_name, coordinates=coord.SkyCoord(*src_coords[:, 0, i], unit=u.rad),
                                       source_type=src_type)

        self._sources = tmp_src

        # Save checks...
        for source_group in (self.targets, self.phase_calibrators, self.bandpass_calibrators):
            for a_source in source_group:
                assert a_source in self.sources, f"The source {a_source} defined in the input file was not observed. " \
                                                 f"The available sources are {', '.join(list(self.sources.keys()))}."

        for an_antenna in self.refants:
            assert an_antenna in self.antennas, f"The ref. antenna {an_antenna} defined in the input file is not " \
                                                f"in the array for {self.project_name}. " \
                                                f"The available antennas are {', '.join(self.antennas)}."


    def get_listobs(self):
        """Retrieves and stores the information about the scans in the project (runs listobs in CASA).
        If the file already exists, it will not be overwriten.
        """
        if not self.msfile.exists():
            print(f"The MS file {self.msfile.name} does not exist!")
            return False

        if (self.logdir / f"{self.project_name.lower()}.listobs").exists():
            print(f"{self.project_name.lower()}.listobs already exists and will not be overwriten.")
        else:
            casatasks.listobs(vis=str(self.msfile), listfile=str(self.logdir / f"{self.project_name.lower()}.listobs"))

    def store(self, path: Path=None):
        """Stores the current Experiment into a file in the indicated path. If not provided,
        it will be '.{expname.lower()}.obj' where exp is the name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'wb') as file:
            pickle.dump(self, file)

    def load(self, path: Path=None):
        """Loads the current Experiment that was stored in a file in the indicated path. If path is None,
        it assumes the standard path of '.{exp}.obj' where exp is the name of the experiment.
        """
        if path is not None:
            self._local_copy = path

        with open(self._local_copy, 'rb') as file:
            obj = pickle.load(file)

        return obj












