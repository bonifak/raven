from distutils.core import setup, Extension
import os

# Replicating the methods used in the RAVEN Makefile to find CROW_DIR,
# If the Makefile changes to be more robust, so should this
# We should be doing a search for CROW, I would think, we should not force a
# directory structure
CURR_DIR = os.path.dirname(os.path.realpath(__file__))
HERD_TRUNK_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')
CROW_SUBMODULE = os.path.join(CURR_DIR,'crow')
if os.path.isfile(os.path.join(CROW_SUBMODULE,'Makefile')):
  CROW_DIR = CROW_SUBMODULE
else:
  CROW_DIR = os.path.join(HERD_TRUNK_DIR,'crow')

BOOST_INCLUDE_DIR = os.path.join(CROW_DIR,'contrib','include')
RAVEN_INCLUDE_DIR = os.path.join('include','contrib')

include_dirs=[RAVEN_INCLUDE_DIR,BOOST_INCLUDE_DIR]
swig_opts=['-c++','-I'+RAVEN_INCLUDE_DIR, '-I'+BOOST_INCLUDE_DIR]
setup(name='amsc',
      version='0.0',
      description='A library for computing the Approximate Morse-Smale Complex (AMSC)',
      ext_modules=[Extension('_amsc',['src/contrib/amsc.i',
                                      'src/contrib/UnionFind.cpp',
                                      'src/contrib/AMSC.cpp'],
                   include_dirs=include_dirs, swig_opts=swig_opts)])
