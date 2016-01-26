# try using distribute or setuptools or distutils.
try:
    import distribute_setup
    distribute_setup.use_setuptools()
except ImportError:
    pass

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


import sys
import re

# parse version from package/module without importing or evaluating the code
with open('thapbi_santi/__init__.py') as fh:
    for line in fh:
        m = re.search(r"^__version__ = '(?P<version>[^']+)'$", line)
        if m:
            version = m.group('version')
            break

if sys.version_info <= (2, 5):
    sys.stderr.write("ERROR: santi_script requires Python Version 2.6 " +
                     "or above...exiting.\n")
    sys.exit(1)

setup(
    name="santi_script",
    version=version,
    author="Leighton Pritchard",
    author_email="leighton.pritchard@hutton.ac.uk",
    description=''.join(["This script translates Santiago Garcia's "
                         "metabarcoding pipeline",
                         "into a Python script with associated modules,"
                         "designed for flexibility and further development"]),
    license="MIT",
    keywords="genome bioinformatics sequence sequencing metabarcoding",
    platforms="Posix; MacOS X",
    url="http://widdowquinn.github.io/THAPBI/santi_script",  # project home page
    download_url="https://github.com/widdowquinn/THAPBI/releases",
    scripts=['thapbi_santi_otus.py'],
    packages=['thapbi_santi'],
    install_requires=['biopython'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    )
