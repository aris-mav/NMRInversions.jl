## Intro

This page will give you a basic idea about how to use the NMRInversions package.
For more details, it's best to refer to the [functions](functions.md) page.

!!! hint
    All of the commands mentioned below should be typed in the julia console, 
    or saved in a text file with the .jl extension, to be used 
    from a terminal with `julia file.jl`, or through an IDE such as VSCode.


# Performing an inversion

Suppose we're working with data coming from a Spinsolve instrument
Then we can do the following:

(you can find some example data in the files of this package, 
look at the github page for the example_data folder).

```julia
using NMRInversions

data = import_spinsolve()
```

!!! info
    Since we called the `import_spinsolve` function without an argument, 
    it'll open a file dialog for us to select the files we want to import.
    (note that `import_spinsolve` requires two files, the `aqcu.par` file,
    plus the file containing the experiment data.)

Now we have the data imported, the inversion can be performed using a single line of code!

```julia
results = invert(data)
```

!!! info
    The `results` variable above is an `inv_out_1D` or `inv_out_2D` structure, 
    which contains all the relevant information produced by the `inversion` function.
    To access that information, we can look at the fields of the structure using the dot notation.
    The field names can be shown by using the REPL help mode (typing ? at the julia> prompt), 
    and typing the variable's name (in this case ?results). 
    Alternatively, running `@doc results` will also give you the same answers.
The results can easily be visualised through the GLMakie extension of the package.

```julia
using GLMakie
plot(results)
```

!!! info
    The `plot` function of GLMakie is modified by this package 
    to work with results from the invert function as arguments.
    It's really easy to use, but if you want more control 
    on how your plots look, it's best to create them from scratch 
    using all the tools available in GLMakie.

Or, if the plot is all you need, one line of code is enough:

```julia
using NMRInversions, GLMakie
plot(invert(import_spinsolve()))
```

Note that the workflow above can work for both 1D and 2D inversions!


# Using the expfit function

In a similar way, we can perform various exponential 
fits to the imported data using the `expfit` function.

```julia
using NMRInversions, GLMakie

data = import_spinsolve()

a = expfit(1, data)  # mono-exponential fit
b = expfit(2, data)  # bi-exponential fit

plot(a,b)  # Visualize both on the same plot

a.eqn  # Print the equation of the mono-exponential fit
b.eqn  # Print the equation of the bi-exponential fit
```
