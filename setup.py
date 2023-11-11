"""Script to setup necessary dependencies for SWAMP and manage unit testing"""

import setuptools
from crops import __prog__, __description__, __version__
from crops import __author__, __date__, __copyright__

import glob
import logging
import os
import sys
from unittest import TestLoader, TextTestRunner, TestSuite
import argparse

from distutils.util import convert_path

ROOT_DIR = os.path.join(os.path.dirname(__file__), 'crops', '')
PACKAGES = ["core", "iomod", "elements"]
PYTHON_EXE = None
for arg in sys.argv:
    if arg[0:20] == "--script-python-path" and len(arg) == 20:
        option, value = arg, sys.argv[sys.argv.index(arg) + 1]
        PYTHON_EXE = value
    elif arg[0:20] == "--script-python-path" and arg[20] == "=":
        option, value = arg[:20], arg[21:]
        PYTHON_EXE = value

if not PYTHON_EXE:
    PYTHON_EXE = sys.executable

def parse_arguments():
    """Parse command line arguments"""

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='CROPS-SETUP: Setup necessary dependencies'
                                                 ' for CROPS and manage unit testing')
    parser.add_argument("--tests", action='store_true',
                        help='If set, run CROPS unit testing')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    args = parser.parse_args()

    return args


def scripts():
    # Modified from conkit.org
    extension = ".bat" if sys.platform.startswith("win") else ""
    header = "" if sys.platform.startswith("win") else "#!/bin/sh"
    bin_dir = "bin"
    command_dir = convert_path("crops/command_line")
    scripts = []
    for file in os.listdir(command_dir):
        if not file.startswith("_") and file.endswith(".py"):
            # Make sure we have a workable name
            f_name = os.path.basename(file).rsplit(".", 1)[0]
            for c in [".", "_"]:
                new_f_name = f_name.replace(c, "-")
            # Write the content of the script
            script = os.path.join(bin_dir, new_f_name + extension)
            with open(script, "w") as f_out:
                f_out.write(header + os.linesep)
                # BATCH file
                if sys.platform.startswith("win"):
                    string = "@{0} -m conkit.command_line.{1} %*"
                # BASH file
                else:
                    string = '{0} -m conkit.command_line.{1} "$@"'
                f_out.write(string.format(PYTHON_EXE, f_name) + os.linesep)
            os.chmod(script, 0o777)
            scripts.append(script)
    return scripts

class UnittestFramework(object):
    """Framework to run unit testing"""

    def run(self, buffer=False, pattern="test*.py", verbosity=2):
        """Main routine for running the test cases"""

        tests = self.find_tests(ROOT_DIR, pattern=pattern)
        if int(tests.countTestCases()) <= 0:
            msg = 'Could not find any tests to run in directory: {0}'.format(ROOT_DIR) + os.linesep
            sys.stderr.write(msg)
            sys.exit(1)
        logging.disable(logging.CRITICAL)
        result = TextTestRunner(verbosity=verbosity, buffer=buffer).run(tests)
        logging.disable(logging.NOTSET)
        if result.wasSuccessful():
            exit(0)
        else:
            exit(1)

    def find_tests(self, directory, pattern="test*.py"):
        """Load a unittest test suite"""

        search_pattern = os.path.join(directory, "*")
        cases = [os.path.basename(folder) for folder in glob.iglob(search_pattern)
                 if os.path.isdir(folder) and os.path.basename(folder) in PACKAGES]

        return self._load_suite(cases, pattern, directory)

    def _load_suite(self, cases, pattern, directory):
        suite = TestSuite()
        for case in cases:
            path = os.path.join(directory, case, "tests")
            try:
                _suite = TestLoader().discover(path, pattern=pattern, top_level_dir=directory)
                suite.addTests(_suite)
                del _suite
            except ImportError:
                print("*** not a package: {0} ***".format(path))
        return suite


if __name__ == "__main__":
    args = parse_arguments()

    if args.tests:
        test = UnittestFramework()
        test.run()
    else:
        with open("README.rst", "r") as fh:
            long_description = fh.read()

        with open("requirements.txt", "r") as fh:
            dependencies = fh.read().splitlines()

        setuptools.setup(
             name=__prog__,
             version=__version__,
             scripts=[os.path.join('bin', 'crops-renumber'),
                      os.path.join('bin', 'crops-cropseq'),
                      os.path.join('bin', 'crops-cropstr'),
                      os.path.join('bin', 'crops-splitseqs')] ,
             author=__author__,
             author_email="J.J.Burgos-Marmol@liverpool.ac.uk",
             description=__description__,
             long_description=long_description,
             long_description_content_type="text/x-rst",
             url="https://github.com/rigdenlab/crops",
             packages=setuptools.find_packages(),
             classifiers=[
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: BSD License",
                 "Operating System :: OS Independent",
                 "Programming Language :: Python",
                 "Programming Language :: Python :: 3.7",
                 "Programming Language :: Python :: 3.8",
                 "Programming Language :: Python :: 3.9",
                 "Programming Language :: Python :: 3.10",
                 "Programming Language :: Python :: 3.11",
                 "Topic :: Scientific/Engineering :: Bio-Informatics",
             ],
             install_requires=dependencies)
