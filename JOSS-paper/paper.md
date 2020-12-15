---
title: '(py)oscode: fast solutions of oscillatory ODEs'
tags:
  - Python
  - C++
  - numerical methods
  - ordinary differential equations
authors:
  - name: Fruzsina Julia Agocs
    orcid: 0000-0002-1763-5884
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: Astrophysics Group, Cavendish Laboratory, J. J. Thomson Avenue, Cambridge, CB3 0HE, UK
   index: 1
 - name: Kavli Institute for Cosmology, Madingley Road, Cambridge, CB3 0HA, UK
   index: 2
date: 7 October 2020
bibliography: paper.bib

---

# Summary

Oscillatory differential equations are ubiquitous in physics, chemistry and beyond. They arise in
quantum mechanics, electrical circuitry, suspension systems, molecular dynamics,
and in models of gravitational and electromagnetic waves.
The numerical solution of such systems however can be a computational bottleneck when tackled with conventional methods
available from numerical libraries. 

We present `(py)oscode`, a general-purpose numerical routine for solving a class of highly
oscillatory ordinary differential equations (ODEs) efficiently. The package has
been designed to solve equations which describe a single harmonic oscillator
with a time-dependent frequency and damping term, i.e. are of the form
\begin{equation}\label{eq:eom}
y'' + 2\gamma(x) y' + \omega^2(x) y = 0.
\end{equation}
The frequency $\omega(x)$ and damping $\gamma(x)$ terms do not need
to be explicit functions of $x$ (they can instead be e.g. the result of another
numerical solution of an ODE), as they are supplied as sequences $\omega_j,
\gamma_j$ evaluated at $x_i \leq x_j \leq x_f$, where $(x_i, x_f)$ is the
integration range.

`(py)oscode` is written in C++, but comes with a Python wrapper.
Its Python interface was designed to be similar to those included in `SciPy`'s [@scipy] numerical ODE solution
modules. This is demonstrated in the example below whose output is shown in
\autoref{fig:airy}.

```python
import numpy as np
import scipy.special as sp
import pyoscode

# Set up the Airy equation as an example: y'' + xy = 0
xs = np.linspace(0,40.0,5000)
ws = np.sqrt(xs)
gs = np.zeros_like(xs)
# Initial conditions
xi = 1.0
xf = 40.0
yi = sp.airy(-xi)[0]
dyi = -sp.airy(-xi)[1]
# Get dense output at the following points
t_eval = np.linspace(15,35,600)
# Solve the equation
solution = pyoscode.solve(xs, ws, gs, xi, xf, yi, dyi, t_eval=t_eval)
```

![Numerical solution of the Airy equation, $y'' + xy = 0$, with `pyoscode`. The
increase in step-size of `pyoscode`'s internal steps (orange dots) is due to the
algorithm switching from using the RK method to the WKB approximation in the presence of high-frequency
oscillations. The orange segment shows dense output, the solution at these
points was computed at no additional evaluations of terms in the differential
equation. \label{fig:airy}](../examples/images/airy.png)

# Statement of need 

Even if the terms in \autoref{eq:eom} change slowly, if the frequency of
oscillations in the solution is high enough, standard numerical methods struggle
to solve such equations quickly. Traditional methods have to trace every
oscillation in the solution, taking many steps in $x$ at an enormous
computational cost. The algorithm underlying `(py)oscode`, published in
@oscode and based on @rkwkb-handley, can detect when the solution is oscillatory and switch to a method
based on an analytic approximation (Wentzel--Kramers--Brillouin, WKB) suited for
oscillatory functions, otherwise using a Runge--Kutta (RK) method. Using the WKB
approximation allows the algorithm to skip over several wavelengths of
oscillation in a single step, reducing the number of steps taken drastically. It
adaptively updates its step-size to keep the local numerical error within a
user-specified tolerance. `(py)oscode` is capable of producing a solution estimate
at an arbitrary value of $x$, not just at its internal steps, therefore it can
be used to generate a "continuous" solution, or dense output [@dense-output]. 

# Related research and software

`(py)oscode`'s development was motivated by the need for a significantly more
efficient solver for the evolution of early-universe quantum fluctuations. These
perturbations are thought to have been stretched to macroscopic scales by a
phase of accelerated expansion of the universe (cosmic inflation), to later become the
large-scale structure we see today. To understand the origins of structure it
is therefore essential to model the perturbations and understand the physics
involved in inflation. `(py)oscode` has been used to speed up the numerical evolution of quantum
fluctuations in the early universe, enabling the exploration of models beyond
the standard model of cosmology [@pps-curved]. It served as inspiration for
other numerical methods aiming to extend the range of oscillatory ODEs to solve
[@beyond-rkwkb]. 

The efficient solution of oscillatory ODEs is a long-standing
numerical analysis problem with many existing methods to handle certain
sub-classes of equations. Examples include successful methods by Petzold [@petzold], reviewed in @petzold-review with many references therein, 
Iserles et al. [@condon-deano-iserles; @deano-integrals; @condon-et-al-circuits], and Bremer [@bremer], with code available from @bremer-code.

# Acknowledgements

I thank Lukas Hergt for invaluable discussions during the early development of
`(py)oscode` and his ongoing support. Construction of the algorithm would not have been possible
without the help and guidance of Will Handley, Mike Hobson, and Anthony Lasenby. 
I was supported by the Science and Technology Facilities Council (STFC).

# References
