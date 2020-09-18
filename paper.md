---
title: '(py)oscode: fast solutions of oscillatory ODEs'
tags:
  - Python
  - C/C++
  - numerical methods
  - ordinary differential equations
  - oscillatory 
  - runge-kutta 
authors:
  - name: Fruzsina Agocs^[corresponding author.]
    orcid: 0000-0002-1763-5884
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Will Handley
    orcid: 0000-0002-5866-0445
    affiliation: "1, 2"
  - name: Anthony Lasenby
    orcid: 0000-0002-8208-6332 
    affiliation: "1, 2"
  - name: Mike Hobson
    orcid: 0000-0002-0384-0182 
    affiliation: 1
affiliations:
 - name: Astrophysics Group, Cavendish Laboratory, J. J. Thomson Avenue, Cambridge, CB3 0HE, UK
   index: 1
 - name: Kavli Institute for Cosmology, Madingley Road, Cambridge, CB3 0HA, UK
   index: 2
date: 13 August 2017
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
`(py)oscode` is in C++, but comes with a Python wrapper.

# Statement of need 

Even if the terms in \autoref{eq:eom} change slowly, if the frequency of
oscillations in the solution is high enough, standard numerical methods struggle
to solve such equations quickly. Traditional methods have to trace every
oscillation in the solution, taking many steps in $x$ at an enormous
computational cost. The algorithm underlying `(py)oscode`, published in
[@oscode], can detect when the solution is oscillatory and switch to a method
based on an analytic approximation (Wentzel--Kramers--Brillouin, WKB) suited for
oscillatory functions, otherwise using a Runge--Kutta (RK) method. Using the WKB
approximation allows the algorithm to skip over several wavelengths of
oscillation in a single step, reducing the number of steps taken drastically. It
adaptively updates its step-size to keep the local numerical error within a
user-specified tolerance. `(py)oscode` is capable of producing a solution estimate
at an arbitrary value of $x$, not just at its internal steps, therefore it can
be used to generate a "continuous" solution. 

# Related research and software

`(py)oscode`'s development was motivated by the need for a significantly more
efficient solver for the evolution of early-universe quantum fluctuations. These
perturbations are thought to have been stretched to macroscopic scales by a
phase of accelerated expansion, cosmic inflation, to later become the
large-scale structure we see today. To understand the origins of structure, it
is therfore essential to model the perturbations and understand the physics
involved in inflation. `(py)oscode` has been used to speed up the numerical evolution of quantum
fluctuations in the early universe, enabling the exploration of models beyond
the standard model of cosmology [@curved-pps]. It served as inspiration for
other numerical methods aiming to extend the range of oscillatory ODEs to solve
[@beyond-rkwkb]. 

The efficient solution of oscillatory ODEs is a long-standing
numerical analysis problem with many existing methods to handle certain
sub-classes of equations. Examples include successful methods by Petzold et al. [@petzold], reviewed in [@petzold-review] with many references therein, 
Iserles et al. [@condon-deano-iserles] [@deano-integrals] [@hu-et-al-circuits], and Bremer [@bremer].

# Acknowledgements

We thank Lukas Hergt for invaluable discussions during the early development of
`(py)oscode` and his ongoing support. FA was supported by the Science and
Technology Facilities Council (STFC). WH thanks Gonville & Caius College for
their continuing support via a college research fellowship.

# References
