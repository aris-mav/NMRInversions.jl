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
    (note that `import_spinsolve` requires two files, the `aqcu.par` file
    containing aqcuistion parameters, plus the `.dat` or `.csv` file 
    containing the experiment data.)

Now we have the data imported, the inversion can be performed using a single line of code!

```julia
results = invert(data)
```

!!! info
    The `results` variable above is an `inv_out_1D` or `inv_out_2D` structure, 
    which contains all the relevant information produced by the `inversion` function.
    To access that information, we can look at the fields of the structure using the dot notation.
    The field names contained in the structure can be shown by using the REPL help mode 
    (typing ? at the julia> prompt), and typing the variable's name (in this case, `?results`). 
    Alternatively, running `@doc results` will also give you the same answers.

The results can easily be visualised through the GLMakie extension of the package.

```julia
using GLMakie
plot(results)
```
This will open a GUI with tools to interactively extract some information from the inversion results,
by selecting regions and labelling them accordingly.

!!! info
    The `plot` function of GLMakie is modified by this package 
    to work with results from the invert function as arguments.
    It's really easy to use, but if you want more control 
    on how your plots look, it's best to create them from scratch 
    using all the tools available in GLMakie.

The process above can also be achieved by a single line of code:
```julia
using NMRInversions, GLMakie
plot(invert(import_spinsolve()))
```

Note that the workflow above can work for both 1D and 2D inversions!

This is how some full examples would look like:

```julia
using NMRInversions, GLMakie

path = ".../NMRInversions.jl/example_data/csv_files/graphene_CPMG.csv"
data = import_csv(IR, path)
results = invert(data)
plot(results)
```
The resulting plot will look like:

![Resulting plot](./assets/1D_gui.png)

Notice that benath the ``T_2`` distribution there's a slider.
You can move the ends of it to select a region within the limits
defined by the red veritical lines.
Then you can use the following options:
- `Label current selection` will highlight the selected region 
  and add some text in the plot with the weighed average ``T_2``
  of that region.
- `Filter-out current selection` will remove the selected region 
  from the distribution, and it will update the fit and the residuals 
  accordingly on the plot.
- `Reset selections` gets you back where you started, removes any 
  selections and brings back filtered-out regions.
- `Save and exit` will bring up a window so that you can save your 
   plot as a .png (without the buttons and the slider).

Let's look at a 2D example as well:


```julia
using NMRInversions, GLMakie

paths = [".../NMRInversions.jl/example_data/spinsolve_IRCPMG/T1IRT2.dat",
         ".../NMRInversions.jl/example_data/spinsolve_IRCPMG/aqcu.par"]
data = import_spinsolve(paths)
results = invert(data)
plot(results)
```
![Resulting plot](./assets/2D_gui.png)

Similarly, now we can select regions by left-clicking at points within 
the plot which enclose a region of interest. Then, we can :

- `Label current selection` will highlight the selected polygon 
  with a dashed line and add some text in the plot with the weighed 
  average ``T_1/T_2`` of that region, as well as the volume fraction of it.
- `Filter-out unselected` will remove everything outside the polygon.
- `Cancel current selection` will discard the current polygon in 
   case you do not like what you selected.
- `Reset everything` will get you back to the original state, 
   removing selections and restoring filtered-out regions.
- `Save and exit` will bring up a window so that you can save your 
   plot as a .png (without the buttons).

There are also some options to change the appearance of the plot, in 
terms of colormap, contour levels and toggling between filled and non-filled
contours. You can also add a title to the plot, which will be saved.

If you have multiple results, you can pass them as a matrix or vector 
into the `plot()` function as:

```julia
plot([results  results ; resutls results])
```
![Resulting plot](./assets/multiple_plots.png)

Of course, it would be more interesting if it wasn't four
copies of the same results, but you get the point.

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
