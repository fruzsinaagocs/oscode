========================================================================
oscode: Oscillatory ordinary differential equation solver
========================================================================

.. image:: https://codecov.io/gh/fruzsinaagocs/oscode/branch/joss-paper/graph/badge.svg
        :target: https://codecov.io/gh/fruzsinaagocs/oscode
        :alt: codecov status
.. image:: https://github.com/fruzsinaagocs/oscode/actions/workflows/python-package.yml/badge.svg
        :target: https://github.com/fruzsinaagocs/oscode/actions
        :alt: GH workflow status
.. image:: https://readthedocs.org/projects/oscode/badge/?version=latest
        :target: https://oscode.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status
.. image:: https://badges.gitter.im/oscode-help/community.svg
        :target: https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
        :alt: Chat on gitter
.. image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
        :target: https://opensource.org/licenses/BSD-3-Clause
        :alt: BSD 3-clause license
.. image:: https://img.shields.io/pypi/dm/pyoscode?color=indigo 
        :target: https://pypi.org/project/pyoscode/
        :alt: PyPI downloads
.. image:: https://joss.theoj.org/papers/d4c9396ef9b2b595e2f3881a4f8a7cda/status.svg
        :target: https://joss.theoj.org/papers/d4c9396ef9b2b595e2f3881a4f8a7cda
        :alt: JOSS status
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4322958.svg
        :target: https://doi.org/10.5281/zenodo.4322958
        :alt: Zenodo doi

|
|

.. contents::
   :local:

|

About
-----

Oscode is a C++ tool with a Python interface that solves **osc**\illatory
**o**\rdinary **d**\ifferential **e**\quations efficiently. It is designed to
deal with equations of the form

.. image:: 
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/oscillator.png

where |gamma| (friction term) and |omega| (frequency) can be given as arrays.

.. |gamma| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/gamma.png

.. |omega| image:: https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/omega.png

Oscode makes use of an analytic approximation of x(t) embedded in a
stepping procedure to skip over long regions of oscillations, giving a reduction
in computing time. The approximation is valid when the frequency changes slowly
relative to the timescales of integration, it is therefore worth applying when
this condition holds for at least some part of the integration range. 

For the details of the numerical method used by oscode, see Citation_.


Installation
------------

Dependencies
~~~~~~~~~~~~

Basic requirements for using the C++ interface:

- C++11 or later
- `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`__ (a header-only library included in this source)

The strictly necessary Python dependencies are automatically installed when you use `pip` or the `setup.py`. They are:

- `numpy <https://pypi.org/project/numpy/>`__

The *optional* dependencies are: 

- for tests
    - `scipy <https://pypi.org/project/scipy/>`__ 
    - `pytest <https://docs.pytest.org/en/stable/getting-started.html>`__ 
- for examples/plotting
    - `matplotlib <https://pypi.org/project/matplotlib/>`__
    - `scipy <https://pypi.org/project/scipy/>`__ 
- for generating offline documentation
    - `sphinx <https://pypi.org/project/Sphinx/>`__ 
    - `doxygen <https://www.doxygen.nl/index.html>`__
    - `breathe <https://pypi.org/project/breathe/>`__
    - `exhale <https://pypi.org/project/exhale/>`__


Python
~~~~~~

``pyoscode`` can be installed via pip 

.. code:: bash
   
   pip install pyoscode

or via the setup.py

.. code:: bash

   git clone --recursive https://github.com/fruzsinaagocs/oscode
   cd oscode
   python setup.py install --user

or

.. code:: bash

   git clone --recursive https://github.com/fruzsinaagocs/oscode
   cd oscode
   pip install .

You can then import ``pyoscode`` from anywhere. Omit the ``--user`` option if
you wish to install globally or in a virtual environment. If you have any
difficulties, check out the `FAQs - Installation
<https://github.com/fruzsinaagocs/oscode#installation-1>`__ section below. 

You can check that things are working by running `tests/` (also ran by Travis continuous integration):

.. code:: bash

   pytest tests/

C++
~~~

``oscode`` is a header-only C++ package, it requires no installation.

.. code:: bash

   git clone --recursive https://github.com/fruzsinaagocs/oscode

and then include the relevant header files in your C++ code:

.. code:: c

    #include <oscode/solver.hpp>
    #include <oscode/system.hpp>


Quick start
-----------

Try the following quick examples. They are available in the `examples
<https://github.com/fruzsinaagocs/oscode/tree/master/examples/>`__.

Python
~~~~~~

:Introduction to pyoscode: |intro_binder|
:Cosmology examples: |cosmology_binder|
:Scipy 2020 lecture notebook: |scipy_binder|

.. |intro_binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/master?filepath=examples/introduction_to_pyoscode.ipynb

.. |cosmology_binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/master?filepath=examples/cosmology.ipynb

.. |scipy_binder| image:: https://mybinder.org/badge_logo.svg
 :target: https://mybinder.org/v2/gh/fruzsinaagocs/oscode/master?filepath=examples/pyoscode_scipy.ipynb


.. image::
    https://github.com/fruzsinaagocs/oscode/raw/master/pyoscode/images/spectra.gif
    :width: 800

C++
~~~

:Introduction to oscode: `examples/burst.cpp`
:To plot results from `burst.cpp`: `examples/plot_burst.py`

To compile and run:

.. code:: bash
    
    cd examples/
    g++ -I../include/ -g -Wall -std=c++11 -c -o burst.o burst.cpp
    g++ -I../include/ -g -Wall -std=c++11 -o burst burst.o
    ./burst


Documentation
-------------

Documentation is hosted at `readthedocs <https://oscode.readthedocs.io>`__.

To build your own local copy of the documentation you can run:

.. code:: bash

   cd pyoscode/docs
   make html

Citation
--------

If you use ``oscode`` to solve equations for a publication, please cite:

- `Efficient method for solving highly oscillatory ordinary differential equations with applications to physical systems <https://doi.org/10.1103/PhysRevResearch.2.013030>`__,
- `Dense output for highly oscillatory numerical solutions  <https://arxiv.org/abs/2007.05013>`__

Contributing
------------

Any comments and improvements to this project are welcome. You can contribute
by:

- Opening and `issue <https://www.github.com/fruzsinaagocs/oscode/issues/>`__ to report bugs and propose new features.
- Making a pull request.

Further help
------------

You can get help by submitting an issue or posting a message on `Gitter <https://gitter.im/oscode-help/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge>`__.

FAQs
----

Installation
~~~~~~~~~~~~

1. Eigen import errors:
    .. code:: bash

       pyoscode/_pyoscode.hpp:6:10: fatal error: Eigen/Dense: No such file or directory
        #include <Eigen/Dense>
                  ^~~~~~~~~~~~~

    Try explicitly including the location of your Eigen library via the
    ``CPLUS_INCLUDE_PATH`` environment variable, for example:

    .. code:: bash

       CPLUS_INCLUDE_PATH=/usr/include/eigen3 python setup.py install --user
       # or 
       CPLUS_INCLUDE_PATH=/usr/include/eigen3 pip install pyoscode

    where  ``/usr/include/eigen3`` should be replaced with your system-specific
    eigen location.

Thanks
------

Many thanks to **Will Handley**, **Lukas Hergt**, **Anthony Lasenby**, and **Mike Hobson** for
their support and advice regarding the algorithm behind `oscode`.
There are many packages without which some part of `oscode` (e.g. testing and
examples) wouldn't run as nicely and smoothly, thank you all developers for
making and maintaining these open-source projects. A special thanks goes to the
devs of `exhale <https://pypi.org/project/exhale/>`__ for making the beautiful C++ documentation possible. 


Changelog
---------
- 1.3.0: 
    - Remove deprecated C API functions from the Python wrapper, `_pyoscode.cpp`. 
- 1.2.0:
    - Update the version of Eigen to 3.4.0
- 1.1.2:
    - Dense output bug fix at the C++ interface 
- 1.1.1: 
    - Support for mac and Windows OS at CI.
- 1.1.0:
    - Users can now define w, g as functions in Python (pyoscode) and call the solver via pyoscode.solve_fn(...)
- 1.0.6:
    - Fix issues related to dense output not being correctly generated, e.g. when timepoints at which dense output was asked for are in descending order, etc. 
- 1.0.5:
    - Fixes related to dense output generation
    - Support for w, g to be given as class member functions in C++
    - Switched to GH actions for continuous integration, and fixed code such that unit tests would run again
    - Minor tweaks
- 1.0.4:
    - set minimally required numpy version: numpy>=1.20.0
    - drop Python 2.7 support, instead support 3.8 and 3.9 in addition to 3.7
- 1.0.3: 
    - paper accepted to JOSS
- 1.0.2:
    - Fixed getting correct numpy include directories
- 1.0.1:
    - Added `pyproject.toml` to handle build dependencies (numpy)
- 1.0.0:
    - Dense output
    - Arrays for frequency and damping term need not be evenly spaced
    - Automatic C++ documentation on readthedocs
    - Eigen included in source for pip installability
    - First pip release :)
- 0.1.2:
    - Bug that occurred when beginning and end of integration coincided
      corrected
- 0.1.1:
    - Automatic detection of direction of integration
- 0.1.0:
    - Memory leaks at python interface fixed
    - C++ documentation added 
