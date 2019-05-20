from distutils.core import setup, Extension
import os
import numpy.distutils.misc_util

pyoscode_module = Extension(
    name="_pyoscode",
    sources=[os.path.join(os.getcwd(),"pyoscode/_pyoscode.cpp")],
    include_dirs=[os.path.join(os.getcwd(),'include')]
    )

setup(
    name="pyoscode",
    version="1.0",
    description="Python interface to oscode 1.0",
    url="",
    author="Fruzsina Agocs, Will Handley, Mike Hobson, and Anthony Lasenby",
    author_email="fa325@cam.ac.uk",
    packages=[],
    install_requires=[],
    extras_require={},
    ext_modules=[pyoscode_module],
    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
    zip_safe=False
)

