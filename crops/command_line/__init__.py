"""This is CROPS: Cropping and Renumbering Operations for PDB structure and Sequence files"""

from crops.about import __prog__, __description__, __author__, __date__, __version__,__copyright__

import logging
import sys


def welcome():
    msg="** Running "+__prog__+" v."+__version__+' **\n'
    msg+="** "+__description__+'\n'
    msg+="** Developed by "+__author__+". Copyright: "+__copyright__+'.\n'
    return msg

def ok():
    msg="** "+__prog__+" finished **\n"
    return msg

def crops_logger(level="info"):
    """Logger setup.

    :param level: Console logging level, defaults to "info".
                      Options: [ notset | info | debug | warning | error | critical ]
    :type level: str, optional
    :param logfile: The path to a full file log, defaults to None
    :type logfile: str, optional
    :return: Configured logger.
    :rtype: :obj:`~logging.Logger`

    """
    class CropsFormatter(logging.Formatter):

        def __init__(self, fmt="%(levelname)s: %(message)s"):
            #logging.Formatter.__init__(self, fmt=fmt)
            super().__init__(fmt="%(levelname)s: %(message)s", datefmt=None, style='%')

        def format(self, record):
            # Save the original format configured by the user
            # when the logger formatter was instantiated
            format_orig = self._style._fmt

            # Replace the original format with one customized by logging level
            if record.levelno == logging.DEBUG:
                self._style._fmt ="%(levelname)s: [%(module)s] (%(lineno)d) %(message)s"
            elif record.levelno == logging.INFO:
                self._style._fmt = "%(message)s"
            else:
                self._style._fmt = "%(levelname)s: %(message)s"

            # Call the original formatter class to do the grunt work
            myformat = logging.Formatter.format(self, record)

            # Restore the original format configured by the user
            self._fmt = format_orig

            return myformat

    logging_levels = {
        "notset": logging.NOTSET,
        "info": logging.INFO,
        "debug": logging.DEBUG,
        "warning": logging.WARNING,
        "error": logging.ERROR,
        "critical": logging.CRITICAL
    }

    # Create logger and default settings
    logging.getLogger().setLevel(logging.DEBUG)

    # create console handler with a higher log level
    custom_fmt=CropsFormatter()
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging_levels.get(level, logging.INFO))
    ch.setFormatter(custom_fmt)
    logging.getLogger().addHandler(ch)

    logging.getLogger().debug("Console logger level: %s", logging_levels.get(level, logging.INFO))

    return logging.getLogger()