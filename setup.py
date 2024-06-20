from __future__ import absolute_import, with_statement, print_function, division
from setuptools import setup, Extension, find_packages
import sys
import os
import platform
import numpy as np

source_dir = os.getenv('OSCODE_EIGEN_INCLUDE_DIR')
def readme(short=False):
    with open("README.rst") as f:
        if short:
            return f.readlines()[1].strip()
        else:
            return f.read()

extra_compile_args = []
if platform.system() == 'Windows':  # For Windows
    if "MSC" in platform.python_compiler():
        # Visual Studio (MSVC)
        extra_compile_args = ['/std:c++17']
    else:
        # assume MinGW or similar
        extra_compile_args = ['-std=c++17', '-Wall']
else:  # For Unix/Linux/MacOS
    extra_compile_args = ['-std=c++17', '-Wall']

pyoscode_module = Extension(
    name="_pyoscode",
    sources=["pyoscode/_pyoscode.cpp"],
    include_dirs=['include','pyoscode',np.get_include(), source_dir],
    depends=["pyoscode/_python.hpp", "pyoscode/_pyoscode.hpp"],
    extra_compile_args=extra_compile_args
    )

setup(
    name="pyoscode",
    version="1.3.0",
    description=readme(short=True),
    long_description=readme(),
    url="https://github.com/fruzsinaagocs/oscode",
    project_urls={"Documentation":"https://oscode.readthedocs.io"},
    author="Fruzsina Agocs",
    author_email="fa325@cam.ac.uk",
    packages=find_packages(),
    install_requires=["numpy"],
    extras_require={"examples:":["matplotlib", "scipy", "jupyter"],
    "docs":["sphinx","sphinx-rtd-theme","numpydoc"], "testing":["pytest"]},
    setup_requires=["pytest-runner","numpy"],
    tests_require=["pytest", "numpy", "scipy"],
    include_package_data=True,
    license="oscode",
    ext_modules=[pyoscode_module],
    headers=["pyoscode/_python.hpp", "pyoscode/_pyoscode.hpp"],
    include_dirs=[np.get_include()],
    keywords="PPS, cosmic inflation, cosmology, oscillatory, ODE",
    classifiers=[
                'Intended Audience :: Developers',
                'Intended Audience :: Science/Research',
                'Natural Language :: English',
                'Programming Language :: Python :: 3.4',
                'Programming Language :: Python :: 3.5',
                'Programming Language :: Python :: 3.6',
                'Programming Language :: Python :: 3.7',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Astronomy',
                'Topic :: Scientific/Engineering :: Physics',
                'Topic :: Scientific/Engineering :: Visualization',
                'Topic :: Scientific/Engineering :: Mathematics'
    ],
)

