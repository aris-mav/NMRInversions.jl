
# Types and Structures

In order to take advantage of Julia's multiple dispatch, 
we need to define some structures which can be used as input to our functions.

!!! info
    Structures are object which may contain multiple fields.
    The fields can be accessed using the dot syntax.
    For example, if we have a structure named `foo`, with two fields, `a`, and `b`, 
    we can access the value of `a` using `foo.a` and the value of `b` using `foo.b`.


## Pulse sequences

The following types are used so that the `invert` function 
can create the appropriate kernel for the inversion. 
Anywhere you see the `seq` keyword, one of the following must be used.

```@docs
IR
SR
CPMG
PFG
IRCPMG
```

## Inversion solvers
These are used to let the invert fucntion know which solver to use.
They can be used as input to the `invert` function as the 'solver' argument.
(e.g., `invert(data, solver=brd)` [Song2002](@cite) or 
`invert(data, solver=pdhgm(10,0.1) )`) [reci2017](@cite).

```@docs
brd
pdhgm
optim_nnls
```

## Finding optimal alpha
These are methods for finding the optimal regularization parameter. 
They can be used as input to the `invert` function as the 'alpha' argument
(e.g., `invert(data, alpha=gcv)` or `invert(data, alpha=lcurve(0.001,1,64))` ).
If you'd like to use a particular value of alpha, 
you can just use that number instead (`invert(data, alpha=1)`).
```@docs
gcv
lcurve
```

## Inversion and expfit inputs

These can be used as inputs containing all of the necessary information
to run the inversion and expfit functions 
(e.g. if `a` is an `input1D` object, `invert(a)` does the job).

```@docs
input1D
input2D
```

##  Inversion and expfit outputs

When you run an inversion function (e.g. ` r = invert(a)`),
the output is an `inv_out_1D` or `inv_out_2D` structure.

```@docs
inv_out_1D
inv_out_2D
expfit_struct
```


## Kernel structure

```@docs
svd_kernel_struct
```
