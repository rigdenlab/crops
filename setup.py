# Minimal installation version. To be updated for cross-platform availability.
import setuptools
from crops import __prog__, __description__, __version__
from crops import __author__, __date__, __copyright__
import os

with open("README.rst", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as fh:
    dependencies = fh.read().splitlines()

setuptools.setup(
     name=__prog__,
     version=__version__,
     scripts=[os.path.join('bin', 'crops-renumber'),
              os.path.join('bin', 'crops-cropseq'),
              os.path.join('bin', 'crops-cropstr')] ,
     author=__author__,
     author_email="J.J.Burgos-Marmol@liverpool.ac.uk",
     description=__description__,
     long_description=long_description,
   long_description_content_type="text/x-rst",
     url="https://github.com/rigdenlab/crops",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: BSD License",
         "Operating System :: OS Independent",
     ],
     install_requires=dependencies
 )
