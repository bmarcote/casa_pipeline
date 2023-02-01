

from . import project
from . import tools


class Plotting(object):
    """Functions to produce plots from a Ms object.
    """
    def __init__(self, ms: project.Ms):
        self._ms = ms

    def jplotter(self, ):
