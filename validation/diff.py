import difflib
import sys

with open('k1000-serial.txt', 'r') as hosts0:
    with open('k1000_mpi.txt', 'r') as hosts1:
        diff = difflib.unified_diff(
            hosts0.readlines(),
            hosts1.readlines(),
            fromfile='hosts0',
            tofile='hosts1',
            n=0,
        )
        for line in diff:
            for prefix in ('---', '+++', '@@'):
                if line.startswith(prefix):
                    break
            else:
                sys.stdout.write(line[1:])
