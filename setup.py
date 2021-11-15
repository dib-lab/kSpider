#!/usr/bin/env python

from distutils.command.build import build
from distutils.spawn import find_executable
import sys
import os
import subprocess
import errno

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python version >=3.6")

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')

try:
    from setuptools import setup, Extension, find_packages
except ImportError:
    from distutils.core import setup, Extension

try:
    with open('README.md') as f:
        readme = f.read()
except IOError:
    readme = ''

if os.path.islink("kSpider_BUILD_DIR"):
    os.unlink("kSpider_BUILD_DIR")

if os.path.exists("build/libkSpider.a"):
    os.symlink("build", "kSpider_BUILD_DIR")


def check_exist(dirs):
    all_exists = True
    not_found_files = list()
    for directory in dirs:
        if not (os.path.isdir(directory)):
            print(f"[ERROR] | DIR: {directory} does not exist.", file=sys.stderr)
            all_exists = False
            not_found_files.append(directory)

    if not all_exists:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), ",".join(not_found_files))


SOURCES = [
    'src/swig_interfaces/kSpider_internal.i',
    'lib/kProcessor/src/kDataFrames/kDataFrame.cpp',
    'lib/kProcessor/src/algorithms.cpp',
]

if not find_executable('swig'):
    sys.exit("Error:  Building this module requires 'swig' to be installed")

INCLUDES = [
    'include',
    'lib/argh',
    'lib/kProcessor/include/kProcessor',
    'lib/kProcessor/ThirdParty/MQF/include',
    'lib/kProcessor/ThirdParty/kmerDecoder/include',
    'lib/kProcessor/ThirdParty/kmerDecoder/lib/kseq/include',
    'lib/kProcessor/ThirdParty/sdsl-lite/include',
    'lib/kProcessor/ThirdParty/ntCard/include',
    'lib/kProcessor/ThirdParty/kmerDecoder/lib/parallel-hashmap',
]

check_exist(INCLUDES)

LINK_ARGS = [
    "-fopenmp",
    "-lgomp",
    "-lbz2",
    "-lz",
    "-ldl",
]

kSpider_BUILD_DIR_dir = "kSpider_BUILD_DIR"

LIBRARIES_DIRS = [
    f"{kSpider_BUILD_DIR_dir}/lib/kProcessor",
    f"{kSpider_BUILD_DIR_dir}",
    f"{kSpider_BUILD_DIR_dir}/lib/kProcessor/ThirdParty/MQF/src",
    "lib/kProcessor/ThirdParty/ntCard",
    f"{kSpider_BUILD_DIR_dir}/lib/kProcessor/ThirdParty/sdsl-lite/lib",
    f"{kSpider_BUILD_DIR_dir}/lib/kProcessor/ThirdParty/kmerDecoder",
    f"{kSpider_BUILD_DIR_dir}/lib/kProcessor/ThirdParty/MQF/ThirdParty/stxxl/lib",
]

check_exist(LIBRARIES_DIRS)

LIBRARIES = [
    'kProcessor',
    'kSpider',
    'sdsl',
    'MQF',
    'ntcard',
    'kmerDecoder',
    'stxxl_debug',
]

SWIG_OPTS = [
    '-c++',
    '-py3',
    '-keyword',
    '-outdir',
    './pykSpider/internal'
]


class CustomBuild(build):
    sub_commands = [
        ('build_ext', build.has_ext_modules),
        ('build_py', build.has_pure_modules),
        ('build_clib', build.has_c_libraries),
        ('build_scripts', build.has_scripts),
    ]


kSpider_module = Extension('_kSpider_internal',
                           library_dirs=LIBRARIES_DIRS,
                           libraries=LIBRARIES,
                           sources=SOURCES,
                           include_dirs=INCLUDES,
                           extra_link_args=LINK_ARGS,
                           extra_compile_args=["-O3", "-Ofast", "-std=c++17", "-fPIC"],
                           swig_opts=SWIG_OPTS,
                           )

classifiers = [
    "License :: OSI Approved :: Apache Software License",
    'Development Status :: 3 - Alpha',
    "Operating System :: POSIX :: Linux",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
]

from setuptools.command.build_py import build_py
class BuildPy(build_py):
    def run(self):
        self.run_command('build_ext')
        super(build_py, self).run()

setup(name='kSpider',
      version="2.0.2",
      description="""A simple yet powerful sequence clustering tool""",
      ext_modules=[kSpider_module],
      py_modules=['kSpider_internal'],
      packages=find_packages('pykSpider'),
      package_dir={'': 'pykSpider'},
      python_requires='>=3.6',
      cmdclass={'build_py': BuildPy},
      license='BSD 3-Clause',
      long_description_content_type='text/markdown',
      long_description=readme,
      classifiers=classifiers,
      install_requires=[
          'Click',
      ],
      include_package_data=True,
      entry_points='''
        [console_scripts]
        kSpider=kSpider2:cli
    ''',
      project_urls={
          'Bug Reports': 'https://github.com/dib-lab/kSpider/issues',
          'Source': 'https://github.com/dib-lab/kSpider/issues',
      },
      )

if os.path.exists("build/libkSpider.a") and os.path.islink("kSpider_BUILD_DIR"):
    os.unlink("kSpider_BUILD_DIR")
