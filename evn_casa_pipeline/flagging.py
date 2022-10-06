from . import tools


class Flagging(object):
    def __init__(self, exp):
        self._exp = exp

    def aoflagger(self, strategy_file: str = None):
        """Runs the AOFlagger on the associated MS.
        If you have a costumized AOflagger strategy file, you can use it.
        """
        if strategy_file is None:
            tools.shell_command("aoflagger", self._exp.msfile.name)
        else:
            tools.shell_command("aoflagger", ["-strategy", strategy_file, self._exp.msfile.name])


