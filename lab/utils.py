""" Definitions and functions for all the projects """

from inspect import getframeinfo, stack

import sys
import time

__all__ = ['time_stamp', 'caller_file_and_line', 'eprint', 'green']

_TIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_'))


# region functions for logging
def time_stamp():
    return _TIME_STAMP


def caller_file_and_line():
    """
    This function is used to trace.
    :return information for caller's filename and line number in the file
    """
    caller = getframeinfo(stack()[1][0])
    return "%s: %d" % (caller.filename, caller.lineno)


# ref: https://stackoverflow.com/questions/5574702/how-to-print-to-stderr-in-python
def eprint(*args, **kwargs):
    """
    print to the stderr
    """
    print(*args, file=sys.stderr, **kwargs)

# endregion


# region functions to add colors
def red(string, e=0):
    return '\033[%s31m%s\033[m' % ('' if e == 0 else '1;', string)


def green(string, e=0):
    return '\033[%s32m%s\033[m' % ('' if e == 0 else '1;', string)


def yellow(string, e=0):
    return '\033[%s33m%s\033[m' % ('' if e == 0 else '1;', string)


def blue(string, e=0):
    return '\033[%s34m%s\033[m' % ('' if e == 0 else '1;', string)


def magenta(string, e=0):
    return '\033[%s35m%s\033[m' % ('' if e == 0 else '1;', string)


def cyan(string, e=0):
    return '\033[%s36m%s\033[m' % ('' if e == 0 else '1;', string)


def white(string, e=0):
    return '\033[%s37m%s\033[m' % ('' if e == 0 else '1;', string)

# endregion
