from __future__ import absolute_import, with_statement, print_function, division
from setuptools import setup, Extension, find_packages
import os
import numpy as np

def readme(short=False):
    with open("README.rst") as f:
        if short:
            return f.readlines()[1].strip()
        else:
            return f.read()

pyoscode_module = Extension(
    name="_pyoscode",
    sources=["pyoscode/_pyoscode.cpp"],
    include_dirs=['include','pyoscode',np.get_include(),'deps/eigen'],
    depends=["pyoscode/_python.hpp", "pyoscode/_pyoscode.hpp"],
    extra_compile_args=['-std=c++17','-Wall']
    )

setup(
    name="pyoscode",
    version="1.3.0",
    description=readme(short=True),
    long_description=readme(),
    url="https://github.com/fruzsinaagocs/oscode",
    project_urls={"Documentation":"https://oscode.readthedocs.io"},
    author="Fruzsina Agocs",
    author_email="fruzsina.agocs@colorado.edu",
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
                'Programming Language :: Python :: 3.x',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Astronomy',
                'Topic :: Scientific/Engineering :: Physics',
                'Topic :: Scientific/Engineering :: Visualization',
                'Topic :: Scientific/Engineering :: Mathematics'
    ],
)

