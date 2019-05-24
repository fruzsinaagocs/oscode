from setuptools import setup, Extension, find_packages
import os
import numpy.distutils.misc_util

def readme(short=False):
    with open("README.rst") as f:
        if short:
            return f.readlines()[1].strip()
        else:
            return f.read()

def get_version(short=False):
    with open("README.rst") as f:
        for line in f:
            if ":Version:" in line:
                ver = line.split(":")[2].strip()
                if short:
                    subver = ver.split(".")
                    return "%s.%s" % tuple(subver[:2])
                else:
                    return ver

pyoscode_module = Extension(
    name="_pyoscode",
    sources=[os.path.join(os.getcwd(),"pyoscode/_pyoscode.cpp")],
    include_dirs=[os.path.join(os.getcwd(),'include')]
    )

setup(
    name="pyoscode",
    version="1.0",
    description=readme(short=True),
    long_description=readme(),
    url="https://github.com/fruzsinaagocs/oscode",
    author="Fruzsina Agocs, Will Handley, Mike Hobson, and Anthony Lasenby",
    author_email="fa325@cam.ac.uk",
    packages=find_packages(),
    install_requires=["numpy", "scipy"],
    extras_require={"plotting:": "matplotlib", "docs":"sphinx"},
    tests_require=["pytest"],
    include_package_data=True,
    license="oscode",
    ext_modules=[pyoscode_module],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
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
                'Topic :: Scientific/Engineering :: Mathematics',
    ],
)

