#!/usr/bin/env python
#-*- coding: utf-8 -*-

from distutils.core import setup
from distutils.extension import Extension
import glob, os, os.path
import util.templating

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

compiler = "cython" # "pyrex"
try:
    if compiler == "pyrex":
        from Pyrex.Distutils import build_ext
        has_pyrex = True
    elif compiler == "cython":
        from Cython.Distutils import build_ext
        has_pyrex = True
    else:
        raise ValueError("Invalid compiler '%s'." % compiler)
except ImportError:
    has_pyrex = False

if has_pyrex:
    source_ext = '.pyx'
    cmdclass    = {'build_ext': build_ext}
else:
    source_ext = '.c'
    cmdclass    = {}

# gather the release details
release = {}
execfile(os.path.join(os.path.dirname(__file__), 'pywt','release_details.py'), {}, release)

# tune the C compiler settings
extra_compile_args = ['-Wall', '-finline-limit=1', '-O2']
#extra_compile_args += ['-march=pentium4',  '-mtune=pentium4']
#extra_compile_args += ['-Wno-long-long', '-Wno-uninitialized', '-Wno-unused']

macros = [('PY_EXTENSION', None),
          #('OPT_UNROLL2', None), # enable some manual unroll loop optimizations
          #('OPT_UNROLL4', None)  # enable more manual unroll loop optimizations
         ]

dwt = Extension("pywt._pywt",
        sources = [(n + source_ext) for n in ['src/_pywt']] + ["src/common.c", "src/convolution.c", "src/wavelets.c", "src/wt.c"], 
        include_dirs = ['src'],
        library_dirs = [],
        runtime_library_dirs = [],
        libraries = [],
        define_macros = macros,
        extra_compile_args = extra_compile_args,
		extra_link_args = [],
		export_symbols = [],
)
 
ext_modules = [dwt]
packages =  ['pywt']
package_dir = {'pywt':'pywt'}

    
def do_setup(**extra_kwds):
    util.templating.expand_files('src/*.template', True)

    setup(
        name = release["name"],
        version = release["version"],
        description = release["description"],
        long_description = release["long_description"],
        author = release["author"],
        author_email = release["author_email"], 
        url = release["url"],
        download_url = release["download_url"],
        license = release["license"],
        keywords = release["keywords"],
        platforms = release["platforms"],
        classifiers = release["classifiers"],
        
        ext_modules = ext_modules,
        
        packages = packages,
        package_dir = package_dir,
        #script_args = ["build_ext"],
        
        cmdclass = cmdclass,
        **extra_kwds
    )

if __name__ == '__main__':
    do_setup()
