#!/usr/bin/env python

# setup for FisherExact library packages

from __future__ import absolute_import
from __future__ import print_function
import os
import sys
from distutils import spawn
try:
    from setuptools import find_packages
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
    from setuptools import find_packages

from numpy.distutils.core import Extension as Ext
from numpy.distutils.core import setup

__version__ = "1.3"
__project__ = "FisherExact"
__author__ = "Emmanuel Noutahi"

sys.path.insert(0, os.path.realpath(os.path.join(os.path.dirname(__file__), "python")))


def configuration(top_path=''):
    fexact_sources = [
        'FisherExact/statlib/fexact.pyf',
        'FisherExact/statlib/FEXACT.F90',
        'FisherExact/statlib/prterr.f'
    ]

    fexact = Ext(name='FisherExact.statlib.fexact', sources=[
                 os.path.join(top_path, x) for x in fexact_sources])
    asa159 = Ext(name='FisherExact.statlib.asa159', sources=[
                 os.path.join(top_path, 'FisherExact/statlib', 'asa159.f90')])
    asa205 = Ext(name='FisherExact.statlib.asa205', sources=[
                 os.path.join(top_path, 'FisherExact/statlib', 'asa205.f90')])
    return [fexact, asa205, asa159]


def setup_package():
    if os.path.exists('README.md'):
        README = open('README.md').read()
    else:
        README = ""  # a placeholder, readme is generated on release
    print("\nVersion : %s\n" % __version__)

    fortran_extnsion = configuration()

    setup(
        name=__project__,
        version=__version__,
        maintainer=__author__,
        description="Fishe's Exact test for mxn contingency table",
        url='https://github.com/maclandrol/FisherExact',
        download_url='https://github.com/maclandrol/FisherExact/archive/%s.tar.gz' % __version__,
        author= __author__,
        author_email='fmr.noutahi@umontreal.ca',
        scripts= ['bin/fexact'],
        packages=find_packages(),
        entry_points={'console_scripts': []},
        keywords="statistic fisher independence test",

        long_description=(README + '\n'),
        license='GPL',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Natural Language :: English',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Operating System :: Microsoft :: Windows',
            'Operating System :: POSIX',
            'Operating System :: Unix',
            'Operating System :: MacOS',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],

        setup_requires=['numpy'],
        install_requires=[
            'future>=0.15.0',
            'numpy>=1.8.0',
            'scipy>=0.16.0',
        ],
        ext_modules=fortran_extnsion
    )

if __name__ == '__main__':
    setup_package()
