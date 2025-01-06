# Intro
This page contains the documentation for various useful 
functions in the NMRInversions package.

!!! tip
    From the julia command line, you can enter '?', 
    followed by the name of any function, struct, 
    or object you want to learn more about (try it!).
    After typing `using NMRInversions` in the julia console, 
    this feature will work for all the functions mentioned below.

!!! info 
    In julia, function definitions look like this:
    ```
    foo(x, y, z ; a, b)
    ``` 

    For the example above, `foo` is the name of the function, 
    and the contents of the parenteses are the arguments.  
    Within the parenteses, we got two types of arguments:
    - Positional arguments
     
    Everything that appears before the semicolon `;` 
    (x,y and z in this example) is necessary,
    and must be given in a specific order.
    - Keyword arguments.

    Everything that appears after the semicolon `;` 
    (b and a in this example) is optional,
    and can be given in any order, but its name must be specified.

    In the example above, we can call the function by typing `foo(1, 2, 3)`,
    (in which case, x=1, y=2, z=3, and the default value for `a` and `b` will be used). 
    You can also call the function by typing `foo(1, 2, 3, a=3)`, to specify the value for `a`, 
    or by typing `foo(1, 2, 3, b=3, a = 2)`, to specify the value for both `a` and `b`.

    Sometimes, if there are many keyword arguments, we write 
    the function as foo(x, y ; kwargs...). 
    For the necessary arguments before the semicolon, order matters. 
    For the keyword arguments after the semicolon, order does not matter, 
    but the name of each argument must be specified.
    For more information, please refer to [this link](https://docs.julialang.org/en/v1/manual/functions/).

# Importing data functions
This package offers some functions to import NMR experiment data of various formats.
Alternatively, you can of course import your data however you see fit.
If a format you're working with is not yet supported, 
please [submit an issue](https://github.com/arismavridis/NMRInversions.jl/issues/new) 
and we'll work on it.


The most basic use case would be using data saved in a csv format, 
where there are *only* two columns, 
one for your x-axis (time for relaxation and b-factor for diffusion)
and one for your y-axis (signal intensity).


```@docs
import_csv(::Type{<:pulse_sequence1D},::String)
```

If you're using a spinsolve instrument, you can use the `import_spinsolve` function.
This one requires two files as an input. 
The `aqcu.par` is automatically exported by SpinsolveExpert, 
but you might have to export your data file manually in a csv format.

```@docs
import_spinsolve(::String)
```

For geospec instruments, you can export your raw data as a text file.
That text file can be read by the `import_geospec` function.

```@docs
import_geospec(::String)
```


# Inversion functions
The most important function is `invert()`, which is the main function of the package.
It works as follows:

```@docs
invert(::Type{<:pulse_sequence1D}, ::AbstractArray, ::Vector)
invert(::input1D)
```

Due to julia's multiple dispatch, 
it is possible to define a function with the same name
but different arguments, to achieve different results.


Because of that, the inversion function also works for 2D inversions,
if the following arguments are used instead:

```@docs
invert(::Type{<:pulse_sequence2D}, ::AbstractVector, ::AbstractVector, ::AbstractMatrix)
invert(::input2D)
```

# Finding alpha
Here we provide two options for finding the optimal value for alpha, 
namely Generalized Cross Validation (GCV) or L-curve. 
Generally gcv seems slightly more reliable in NMR, but it's far from 
perfect, so it's good to have alternatives and cross-check.
The following methods can be used as inputs for the `alpha` argument in the
`invert` function:

```@docs
gcv()
gcv(::Real; kwargs...)
gcv(::Real, ::Real ; kwargs...)
lcurve(::Real, ::Real, ::Int)
lcurve(::Real; kwargs...)
lcurve(::Real, ::Real ; kwargs...)
```

The `gcv()` method (from Mithcell 2012) usually involves the least amount of 
function calls and it is thus much faster, thus used as the default option.
If you want more precision, the univariate or box methods should be used instead.

# Exponential fit functions

For 1D data, we can use the `expfit` function to perform multiexponential fits.
We can use the function by specifying either the number of exponential components,
or a vector which defines the starting points for the regression.
See below:

```@docs
expfit(::Union{Int, Vector{<:Real}}, ::Type{<:NMRInversions.pulse_sequence1D}, ::Vector, ::Vector)
expfit(::Union{Int, Vector{<:Real}}, ::NMRInversions.input1D)
```  
  
If you have some rough clue about what results you expect, 
it's best to define some starting points close to these.
(especially important if you're using double or tri-exponential fits.)



# Plotting functions

This package offers plotting capabilities, using its GLMakie extension.
Simply use `using GLMakie` before or after `using NMRInversions`
to load the package, and these functions become available.

We basically take the `plot()` function offered by GLMakie and extend it to types from the package.

For inversions results, if your `inv_out_1D` or `inv_out_2D` struct is saved as
a variable named `r`, using `plot(r)` will run a GUI to interractively visualize
your data and select regions to characterize.

Using a matrix or vector of results (e.g. plot([r1, r2, r3]) will plot all of
them on the same figure, but without the selection options.

For example, for 1D inversions, we have:

```@docs
plot(res::inv_out_1D)
plot(res_mat::VecOrMat{inv_out_1D};kwargs...)
```

And for 2D inversions:

```@docs
plot(::NMRInversions.inv_out_2D)
plot(results::AbstractVecOrMat{inv_out_2D}; kwargs...)
```

For expfits, we have:

```@docs
plot(::NMRInversions.expfit_struct)
plot!(::Union{Makie.Figure,Makie.GridPosition}, ::NMRInversions.expfit_struct )
```

# Miscellaneous functions

Once you have selected some peaks in your inversion results through the GUI,
you might want to extract the weighted averages of these selected peaks,
to get the underlying relaxation times or diffusion coefficients.
The following functions do the job:

```@docs
weighted_averages(r::inv_out_1D)
weighted_averages(r::inv_out_2D)
```
