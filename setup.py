#!/bin/env python3

version="0.0.1"

import os
#os.system('rm -rf *.egg-info/')        # Make 'really clean'

# Prevent the setuptools_scm plugin from adding (only) the contents of the git repo to the tarball:
os.system('mv -f .git .git_temp')

with open("README.md", "r") as fh:
    long_description = fh.read()


from setuptools import setup
setup(
    name='elp-mpp02',
    description='Accurate Moon positions using the Lunar solution ELP/MPP02 in Python',
    author='Marc van der Sluys',
    url='https://github.com/MarcvdSluys/ELP-MPP02',
    
    packages=['elp_mpp02'],
    install_requires=['numpy','fortranformat'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    
    version=version,
    license='GPLv3+',
    keywords=['Moon','Astronomy','Ephemeris'],
    
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Astronomy",
    ]
)

# Put git repo back:
os.system('mv -f .git_temp .git')

# Do some basic checks:
print("\nPython source files included in tarball:")
os.system('tar tfz dist/elp-mpp02-'+version+'.tar.gz |grep -E "\.py"')
print()

os.system('twine check dist/elp-mpp02-'+version+'.tar.gz')
os.system('twine check dist/elp_mpp02-'+version+'-py3-none-any.whl')
print()

