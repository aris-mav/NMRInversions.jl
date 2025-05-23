---
title: 'NMRInversions.jl, a Julia package for time-domain Nuclear Magnetic Resonance'
date: 04 November 2024 
tags:
  - Magnetic Resonance
  - NMR relaxation
  - NMR diffusion
  - Julia
  - Inverse problems
  - Numerical inversion
authors:
  - name: Aristarchos Mavridis
    affiliation: "1"
    orcid: 0000-0002-6619-2303
  - name: Carmine D'Agostino
    affiliation: "1,2"
    orcid: 0000-0003-3391-8320

affiliations:
  - name: Department of Chemical Engineering, The University of Manchester, Oxford Road, Manchester, UK 
    index: 1
  - name: Dipartimento di Ingegneria Civile, Chimica, Ambientale e dei Materiali (DICAM), Alma Mater Studiorum  Università di Bologna, Via Terracini, 28, 40131 Bologna, Italy
    index : 2
bibliography: paper.bib
---

# Summary

This package provides a library of functions as a user-friendly interface 
for performing time-domain Nuclear Magnetic Resonance (NMR) data processing and visualizations.
It is aimed towards NMR users who are not necessarily familiar with the Julia programming language,
or even programming in general.

Functionality includes importing data from various NMR instrument formats,
performing phase-correction on raw data, fitting multi-exponential curves using numerical inversions,
and interactively visualizing the results in a straightforward manner.

It should be noted that the aforementioned inversions should not be confused with 
harmonic inversions [@mandel1997] used in NMR spectroscopy. These can be handled by 
other packages such as "harminv" [@harminvcite], 
and are not yet implemented in NMRInversions.jl at the time of writing.

# Statement of need

NMR relaxation and diffusion methods are popular in various fields of science 
and engineering, with applications including studying the properties of 
porous rocks, catalyst supports and biological tissues, among many others.
In order to extract useful information from such measurements, numerical 
inversions are often used to uncover the different relaxation times 
or diffusion coefficients contributing to the experimental data.

However, working with the aforementioned inversion methods can be quite daunting towards 
users who are not very familiar with linear algebra, optimization methods, or programming 
environments in general.
This package can cover the needs of users who would like a very simple interface for 
everyday data processing, without assuming any knowledge on the topic, while retaining 
the ability to dive into the details and have control over every step of the process 
if needed, once experience is acquired.

The Julia programming language [@bezanson2017julia] is an excellent option for such a package, since it 
provides user-friendly, high level syntax which is similar to MATLAB, without any 
compromise when it comes to computational performance. 
It's a thriving ecosystem for scientific computing applications, 
with several optimization libraries able to tackle the problems arising in NMR applications. 
Julia's capabilities can enable easy integration with the latest advances from the literature, 
without the need for writing solvers in low-level languages for performance purposes, 
in line with Julia's philosophy of solving the "two language problem".
Additionally, Julia's "multiple dispatch" 
allow us to use the same function names for different 
operations, depending on the function inputs.
This enables our package to be more user-friendly, as it 
reduces the number of functions that need to be memorized by the user.

Instrument manufactures often provide algorithms similar to the ones provided in this package, however, 
they are often part of graphical user interfaces, which don't offer customization opportunities, 
and lack the flexibility of writing scripts to efficiently process large sets of data. 

Furthermore, general efforts should be made in improving the reproducibility of scientific research, 
and having accessible, open-source algorithms is of paramount importance when working towards such goals.


# Theory

Numerical inversion methods have been used extensively in the NMR literature for the last fifty years.
Recent advances in computer performance have made these methods more accessible than ever, 
with the average laptop computer nowadays being easily capable of tackling such computations with ease.

Relaxation and diffusion NMR experiments in fluids yield exponentially-decaying signals, 
according to the BPP theory [@BPP].
Characterizing the rate of these decays can provide valuable information on the physical and 
chemical characteristics of the sample.

Usually, there are multiple relaxation times and/or diffusion coefficients within most samples of interest, 
which can represent different pore environments, functional groups, among others. 
Therefore, multi-exponential forms are commonly used to characterize the data.
Recovering the correct distributions is an inverse problem with many possible solutions.
We must therefore rely on regularization methods, where we are looking to minimize
not only the residuals of the fit, but also a penalty term, such as the norm of the solution.

It is very common to refer to these methods as "Inverse Laplace Transforms" (or "ILT's")
in the NMR literature. However, this terminology can be misleading [@Fordham2017],
and thus the more general term "Numerical Inversions" is preferred.

The methods used in this package are mainly based on [@Mitchell2012] and [@Hansen2010]. 
Mathematical notation in the source code mostly follows the aforementioned sources. 

# Acknowledgements
The authors would like to acknowledge funding from bp-ICAM and the EPSRC, 
and support from the ICAM69 project mentors, Mark Sankey and Kuhan Chellappah. 
Furthermore, A. Mavridis would like to thank the Julia community for numerous 
helpful discussions around the topic.

# References
