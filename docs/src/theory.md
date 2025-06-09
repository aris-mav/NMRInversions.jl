What follows is a rough explanation of what this package is about.
If you're looking for more rigor, please refer to sources 
[Hansen2010, Mitchell2012, Song2002, Venka2002](@cite)
in the 
[References](references.md) section. If you believe something mentioned below is wrong
or misleading, please [submit an issue](https://github.com/arismavridis/NMRInversions.jl/issues/new).

## Exponentials in NMR
For NMR relaxation and diffusion experiments in liquids, we expect
the results we get to look more or less like an exponential decay 
or recovery, as described by the BPP theory [BPP](@cite).

To obtain useful information from NMR experiments, we can use different 
exponential forms, according to the pulse sequence used.

Examples are : 

| Pulse sequence    | Equation |
|------------------:|:---------| 
| Inversion recovery | `` M = 1 - 2( \exp(-t/T_1) ) `` |
| Saturation recovery| `` M = 1 - \exp(-t/T_1)  `` |
| CPMG               | `` M = \exp(-t/T_2)  `` |
| PGSE               | `` M = \exp(-b * D)  `` |

, where ``M`` is the magnetization (the "signal" you would record in such an 
experiment, normalized so that its maximum value is 1), ``t`` is time, ``T_1`` 
and ``T_2`` are relaxation time constants, ``b`` is the b-factor in the 
Stejskal-Tanner equations, and ``D`` is diffusion coefficient.

This is all good if the sample you want to study is nice and simple 
(e.g., a liquid made up of a single type of molecule), and you can fit your  
data using the equations above to obtain the desired quantity, which would be 
either ``T_1``, ``T_2``, or ``D``. 

## Multi-exponential forms

If the NMR data comes from a system which is to some degree inhomogenous,
e.g. has multiple pore environments, several different liquids, etc., 
we would expect for the data to be a sum of exponentials, one for each of
the underlying components.

For example, if we do a CPMG experiment on a mixture of water and ethanol,
we would expect the equation to look like :

```math
M = \alpha \exp(-t/T_{2, water}) + \beta \exp(-t/T_{2, ethanol}) \ ,
```
where ``\alpha`` and ``\beta`` would represent the fractions of 
water and ethanol spins into the system, respectively.

Of course, there's no reason to stop at two components.
The more general form could be written as :

```math
M = \sum_{k=1}^n f_k \exp(-t/T_{2,k}) \ ,
```

where ``n`` is the total number of different components hidden behind the data,
and ``f_k`` is the fraction of the ``k``th component.

## Exponential fits

Normally, we would use some type of [nonlinear regression]
(https://en.wikipedia.org/wiki/Nonlinear_regression)
to handle cases where we already know that the system we 
study is made up of a limited number of components 
(``n`` would be 1,2 or maybe 3). 
This package covers this case using the `expfit(n, data)` function.


However, in a lot of realistic scenarios, this method fails for a few reasons:

- We don't always know exactly how many components our system is made out of.
- There might be a virtually continuous range of relaxation times (or diffusion 
  coefficients) within the system, e.g. a fluid within a rock sample with 
  many different pore sizes might have a different relaxation time for each pore size.
- When we add too many degrees of freedom by introducing lots of free 
  variables in the fit equation, it's easy for the algorithm to converge
  in nonsensical values which just happen to be close to the starting points 
  and yield a low-enough value for the cost function. The fit could be satisfying,
  but the values might be unrealistic and unreliable.

## Inversions

In linear algebra terms, a discrete version of the equation :
```math
M = \sum_{k=1}^n f_k \exp(-t/T_{2,k})  \ ,
```
can be written as :
```math
Kf = g \ ,
```
where :
- ``g`` is a vector containing the ``M`` values recorded from an experiment,
- ``f`` is a vector containing the ``f_k`` values present in the system,
- ``K`` is a "kernel" matrix, which can be constructed for a predetermined, usually
  logarithmically-spaced range of ``T_2`` values.

What we are interested in is to solve for the ``f`` vector, since 
``g`` is already known and ``K`` can be constructed according to the range of 
``T_2`` values we're looking for (or ``T_1``, or ``D``).

The "naive" way to approach this would be to use the inverse of the ``K`` matrix, and 
solve the problem as 

```math
f = K^{-1}g \ .
``` 

However the resulting ``f`` would be terrible for a few reasons:

- Even small amounts of noise in the data would mess up the results due to the very large 
  [condition number](https://en.wikipedia.org/wiki/Condition_number) of the ``K`` matrix. 
- Since there **will** be some amount of noise in ``g``, we would not even want our fit 
  curve (``Kf``) to be exactly equal to the data (``g``), since that would involve fiting 
  the noise.
- Negative values in ``f`` would be okay with this method, which does not make much 
  sense since we can't have negative amounts of relaxation times 
  (unless you have a reason to allow that [Chand2018](@cite)).

This is a textbook example of an [inverse problem](https://en.wikipedia.org/wiki/Inverse_problem).

The most commmon way to approach it would be [Tikhonov regularization]
(https://en.wikipedia.org/wiki/Ridge_regression#Tikhonov_regularization),
where we solve the following optimization problem : 

```math
\min_{f \geq 0} \ ||Kf-g||_2^2 + \alpha ||f||_2^2 \ ,
```

where ``\alpha`` is a scalar value referred to as the "regularization parameter".

This way, we minimize the residuals (``||Kf-g||``) while preventing the norm of the  
solution (``||f||``) from getting out of control due to the noise in ``g``.

Low values of ``\alpha`` would result in noisy, unregularized solutions with spikey
``T`` or ``D`` distributions, while large ``\alpha``'s would lead to very broad, over-
smoothed distributions.

The `invert` function in this package solves this problem, and automatically finds the 
optimal value for ``\alpha`` using the GCV method. There is also an option to use the 
`L-curve` method to find ``\alpha``, or supply a constant value for it.
