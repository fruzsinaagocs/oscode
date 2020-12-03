from __future__ import absolute_import, with_statement, print_function, division
from setuptools import setup, Extension, find_packages
import os
import numpy.distutils.misc_util

def readme(short=False):
    with open("README.rst") as f:
        if short:
            return f.readlines()[1].strip()
        else:
            return f.read()

pyoscode_module = Extension(
    name="_pyoscode",
    sources=["pyoscode/_pyoscode.cpp"],
    include_dirs=['include','pyoscode',numpy.distutils.misc_util.get_numpy_include_dirs()],
    depends=["pyoscode/_python.hpp", "pyoscode/_pyoscode.hpp"],
    extra_compile_args=['-std=c++11','-Wall']
    )

setup(
    name="pyoscode",
    version="1.0.1",
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
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    keywords="PPS, cosmic inflation, cosmology, oscillatory, ODE",
    classifiers=[
                'Intended Audience :: Developers',
                'Intended Audience :: Science/Research',
                'Natural Language :: English',
                'Programming Language :: Python :: 2.7',
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

