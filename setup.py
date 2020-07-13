from setuptools import setup
from os import path
from io import open

import nrpcalc

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='nrpcalc',
    
    # Link: https://www.python.org/dev/peps/pep-0440/#version-scheme
    version=nrpcalc.__version__,
    
    description=nrpcalc.__doc__,
    
    long_description=long_description,
    
    long_description_content_type='text/markdown/html',
    
    url='https://github.com/ayaanhossain/nrpcalc',
    
    author=nrpcalc.__authors__,
    
    # author_email='someone@somewhere.com',  # Optional
    
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        # These classifiers are *not* checked by 'pip install'. See instead
        # 'python_requires' below.
        'Programming Language :: Python :: 2.7',
    ],
    keywords=' '.join([
        'synthetic',
        'computational',
        'biology',
        'genetic',
        'parts',
        'calculator',
        'non-repetitive',
        'design',
        'discovery',
        'algorithm',
        'stable',
        'systems',
        'nrp',
        'repeats',
        'vertex',
        'cover',
        'path',
        'finding']),

    packages=['nrpcalc', 'nrpcalc.base'],

    python_requires='>=2.7, <3.0.*',

    install_requires=[
        'numpy>=1.16.6',
        'biopython<=1.76',
        'plyvel>=1.2.0',
        'scipy>=1.2.3',
        'networkx>=2.2'],

    project_urls={  # Optional
        'Bug Reports': 'https://github.com/ayaanhossain/nrpcalc/issues',
        'Source'     : 'https://github.com/ayaanhossain/nrpcalc/',
    },
)