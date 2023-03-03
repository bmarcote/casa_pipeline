import sys
import subprocess
import datetime
from typing import Iterable, Union, Optional, Tuple


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


def shell_command(command: str, parameters: Optional[Union[str, Iterable[str]]] = None, shell: bool = True,
                  bufsize=-1, stdout=None, stderr=subprocess.STDOUT):
    """Runs the provided command in the shell with some arguments if necessary.
    Returns the output of the command, assuming a UTF-8 encoding, or raises ValueError
    if fails. Parameters must be either a single string or a list, if provided.
    """
    if isinstance(parameters, Iterable):
        full_shell_command = [command] + parameters
    else:
        full_shell_command = [command] if parameters is None else [command, parameters]

    print("\n\033[1m> " + f"{' '.join(full_shell_command)}" + "\033[0m")

    process = subprocess.Popen(' '.join(full_shell_command), shell=shell,
                               stdout=stdout, stderr=stderr, bufsize=bufsize)
    # process = subprocess.Popen(full_shell_command, shell=shell, stdout=subprocess.PIPE,
    # for line in process.stdout:
    #     print(line.decode('utf-8').replace('\n', ''))
    output_lines = []
    while process.poll() is None:
        if process.stdout is not None:
            out = process.stdout.readline().decode('utf-8')
            output_lines.append(out)
            sys.stdout.write(out)
            sys.stdout.flush()

    if (process.returncode != 0) and (process.returncode is not None):
        raise ValueError(f"Error code {process.returncode} when running {command} {parameters} in ccs.")

    return ' '.join(full_shell_command), ''.join(output_lines)

def mjd2date(mjd: float) -> datetime.datetime:
    """Returns the datetime for the given MJD date.
    """
    origin = datetime.datetime(1858, 11, 17)
    return origin + datetime.timedelta(mjd)


def date2mjd(date: datetime.datetime) -> float:
    """Returns the MJD day associated to the given datetime.
    """
    origin = datetime.datetime(1858, 11, 17)
    return  (date-origin).days + (date-origin).seconds/86400.0
