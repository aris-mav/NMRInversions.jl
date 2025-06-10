Julia provides the Serialization.jl package to save and load data in a binary format.
This provides us with a straightforward way to save our results in a file, and
load them later if we need to check what's inside them.

!!! info
    Serialization.jl and DelimitedFiles.jl are both part of the 
    standard library, and thus it can be used without installation.

## Saving results 

Assuming that we have a variable `results` which contains our data,
(e.g., `results` is a `inv_out_2D` structure),
we can use the following example to save them in a file.

```julia
using Serialization

# Save the results in a file named "results.dat"
# using the serialize function
serialize("results.dat", results)

# Load the results
# using the deserialize function
results = deserialize("results.dat")

```

Note that we use a ".dat" file extension above, but you can 
use whatever extension you like (or none), as long as it 
makes sense to you.

Writing into plain-text files in delimited formats is also possible,
using the `DelimitedFiles` library and the `writedlm` function.

We have to choose which data to save. In this example below we 
use the relaxation time (or diffusion coefficient) distribution.

```julia
using DelimitedFiles

data_matrix = [results.X results.f]

# The third parameter in the function below
# determines the delimiter, in this case a comma
writedlm("data.csv", data_matrix, ',')
```

This might be useful if you want to export the data into any
non-julia software.
For more elaborate formats, you could have a look at packages 
such as [TOML.jl](https://docs.julialang.org/en/v1/stdlib/TOML/)
or [JSON.jl](https://github.com/JuliaIO/JSON.jl).

## Saving plots

You can use the 'save' function provided by `GLMakie` to save plots.

```julia
using NMRInversions, GLMakie

data = import_csv(IR)
results = invert(data)
p = plot(results)

save("results.png", p)

```
This will save the plot `p` in the file path `results.png`.

If you'd rather use a file dialog to specify the path interactively, you can do something like:

```julia
save(NMRInversions.save_file("", filterlist = "png") , p)

```
This will open a file dialog, where you can interractively select where to save
the image.

Note that for 2D data, the plot GUI already has a "Save and exit" button which 
can be used to save the plot in the same way.
