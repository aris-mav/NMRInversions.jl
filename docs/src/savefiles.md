Julia provides the Serialization.jl package to
save and load data in a binary format. This
provides us with a straightforward way to save our
results in a file, and load them later if we need
to check what's inside them.

!!! info 
    Serialization.jl and DelimitedFiles.jl are
    both part of the standard library, and thus it
    can be used without installation.

## Saving results 

Assuming that we have a variable `results` which
contains our data, (e.g., `results` is a
`InversionData` structure), we can use the
following example to save them in a file.

```julia
using NMRInversions

data = import_csv(IR)
results = invert(data)

using Serialization

# Save the results in a file named "results.dat"
# using the serialize function
serialize("results.dat", results)

# Load the results later
# using the deserialize function
results = deserialize("results.dat")
```

Note that we use a ".dat" file extension above,
but you can use whatever extension you like (or
none), as long as it makes sense to you.

Writing into plain-text, readable files in
delimited formats is also possible, using the
`DelimitedFiles` library and the `writedlm`
function.

We have to choose which data to save. In this
example below we use the relaxation time (or
diffusion coefficient) distribution.

```julia
using DelimitedFiles

data_matrix = [results.axes[1] results.data]

# The third parameter in the function below
# determines the delimiter, in this case a comma
writedlm("data.csv", data_matrix, ',')
```

This might be useful if you want to export the
data into any non-julia software. For more
options, you could have a look at packages such as
[TOML.jl](https://docs.julialang.org/en/v1/stdlib/TOML/)
or [JSON.jl](https://github.com/JuliaIO/JSON.jl).
