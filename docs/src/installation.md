## Installing Julia

You can download the Julia programming language by following [this link](https://Julialang.org/downloads/).
Afterwards, you can either use the Julia console (aka REPL) 
by finding the Julia .exe file in your computer, or you can use VScode,
which provides an environment similar to MATLAB with the use of the 
[Julia extension](https://www.Julia-vscode.org/)
(follow link for installation instructions). 
Of course, if you already have a preferred development workflow 
and you know what you're doing, by all means go for it.

## Installing the package

The package can be installed by running the following command on the Julia console:
```
using Pkg ; Pkg.add("NMRInversions")
```
This usually takes a while, but it needs to be done only once 
(unless you swap environment, more on that 
[here](https://pkgdocs.Julialang.org/v1/environments/)).

Afterwards, you can use the package by running 
```
using NMRInversions
```
in your Julia console, every time you start a new session.


Whenever a new version comes up, you can run:
```
using Pkg ; Pkg.update("NMRInversions")
```
to update the package.

## GLMakie extension
The package provides an extension for interactive visualization,
using the GLMakie package. To gain access to these capabilities,
you also need to install GLMakie:
```
using Pkg ; Pkg.add("GLMakie")
```
And to use it, you need to run
```
using GLMakie
```
in order to access the plotting functions.
