import sys

import subprocess




def shell_command(command, parameters=None, shell=True, bufsize=-1,
                  stdout=None, stderr=subprocess.STDOUT):
    """Runs the provided command in the shell with some arguments if necessary.
    Returns the output of the command, assuming a UTF-8 encoding, or raises ValueError
    if fails. Parameters must be either a single string or a list, if provided.
    """
    if isinstance(parameters, list):
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

