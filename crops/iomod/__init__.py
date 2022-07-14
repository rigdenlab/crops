"""Input/output functions."""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

import os
import argparse
import logging


def check_path(path, typeofpath=None):
    """Return full path if correct.

    :param path: Path to file or directory.
    :type path: str
    :param typeofpath: The type of path, 'dir' or 'file'. All types are considered if not value is given, defaults to None.
    :type typeofpath: str, optional

    :raises ValueError: If given `typeofpath` is neither 'dir' nor 'file'.
    :raises argparse: If given `path` is not found.

    :return: Absolute path.
    :rtype: str

    """
    pathok = False
    if typeofpath == 'dir':
        path = os.path.abspath(path)
        if os.path.isdir(os.path.join(path, '')) is True:
            path = os.path.abspath(os.path.join(path, ''))
            pathok = True
    elif typeofpath == 'file':
        path = os.path.abspath(path)
        if os.path.isfile(path) is True:
            pathok = True
    elif typeofpath is None:
        if os.path.isdir(os.path.abspath(os.path.join(path, ''))) is True:
            path = os.path.abspath(os.path.join(path, ''))
            pathok = True
        else:
            path = os.path.abspath(path)
            if os.path.isfile(path) is True:
                pathok = True
    else:
        logging.critical("Input string 'typeofpath' should be either 'dir' or 'file'.")
        raise ValueError
    if pathok is True:
        return path
    else:
        logging.critical(f"readable_dir:{path} is not a valid path")
        raise argparse.ArgumentTypeError


def outpathgen(globaldir, subdir=None, filename=None, mksubdir=False):
    """Build a filepath in an existing directory and create subdirectory if requested.

    :param globaldir: General output dir.
    :type globaldir: str
    :param subdir: Additional subdirectory, defaults to None.
    :type subdir: str, optional
    :param filename: File name, defaults to None.
    :type filename: str, optional.
    :param mksubdir: Create directory if it does not exist, defaults to False.
    :type mksubdir: bool, optional

    :raises FileNotFoundError: If directory does not exist and mksubdir is False.

    :return: Output filepath.
    :rtype: str

    """
    newpath = check_path(globaldir, 'dir')

    if subdir is not None:
        newpath = os.path.join(newpath, subdir)

        if os.path.isdir(newpath) is False:
            if mksubdir is True:
                os.mkdir(newpath)
            else:
                logging.critical('Directory does not exist')
                raise FileNotFoundError
    if filename is not None:
        newpath = os.path.join(newpath, filename)

    return newpath
