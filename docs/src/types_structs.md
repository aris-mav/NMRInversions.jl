In order to take advantage of Julia's multiple
dispatch, we need to define some structures which
can be used as input to our functions.

!!! info 
    Structures are object which may contain
    multiple fields. The fields can be accessed
    using the dot syntax. For example, if we have
    a structure named `foo`, with two fields, `a`,
    and `b`, we can access the value of `a` using
    `foo.a` and the value of `b` using `foo.b`.


## Pulse Sequences - Data Axes

The following types are used so that the `invert`
function can create the appropriate kernel for the
inversion. Anywhere you see the `seq` keyword, one
of the following types must be used. 

```@docs
IR
SR
CPMG
PFG
Spectrum
FID
FC
```

The types above represent just a special type of
vector. Associating the x-axis of the data with a
pulse sequence seems to be a fairly elegant way to
make the package scale to n-dimensions, thus this
convention is used from version 2.0 onwards.


## Inversion Solvers

These are used to let the invert function know
which solver to use. 

They can be used as input to the `invert` function
as the 'solver' argument. (e.g., `invert(data,
solver=brd())` [Song2002](@cite) or `invert(data,
solver=pdhgm(10,0.1))` [Reci2017](@cite)).

```@docs
brd
pdhgm
cdL1
nnls
optim_nnls
```

The following are also implemented through
external package extensions.

```@docs
jump_nnls
```

## Finding Optimal Alpha

These are methods for finding the optimal
regularization parameter. They can be used as
input to the `invert` function as the 'alpha'
argument (e.g., `invert(data, alpha=gcv())` or
`invert(data, alpha=lcurve(0.001,1,64))` ). If
you'd like to use a particular value of alpha, you
can just use that number instead (`invert(data,
alpha=1)`).

```@docs
gcv
lcurve
```

## Inversion and Expfit Inputs

These can be used as inputs containing all of the
necessary information to run the `invert` and
`expfit` functions (e.g. if `a` is an
`ExperimentData` object, `invert(a)` does the
job).

```@docs
ExperimentData
```

##  Inversion and Expfit Outputs

When you run an inversion function (e.g. ` r =
invert(a)`), the output is an `InversionData`
structure.

```@docs
InversionData
ExpfitData
```

## Kernel Structure

```@docs
Kernel
```
