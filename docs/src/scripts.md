# Writing scripts

The greatest benefit of processing data using
command line interfaces is the ability to write
scripts, automating repetitive processes. This is
especially useful if you're dealing with large
amounts of files.

Below is an example script which imports some data
from a folder full of `.tnt` files, processes
them, and plots them in a presentable manner.

Keep in mind that this one is intentionally
complicated, in order to showcase many of the
things you could do. Chances are, you'd be writing
something shorter for most tasks.

```julia
using NMRInversions, GLMakie

# read all files in current directory
# keep only the ones which contain ".tnt" in their name
files = filter(contains(".tnt"), readdir())

# use broadcasting "." to apply the import function to an array of input strings
data = import_tnt.(CPMG, files)

# use an array comprehension to cut the first 6 data points 
# from each element of the data array
data = [x[6:end] for x in data]

# use broadcasting to invert the whole array of data, 
# returning an array of results
results = invert.(data, lims=(-4, 1, 256))

# remove the .tnt extension from the filename and use it as a title
titles = chopsuffix.(files, ".tnt")

# fill the "title" field in the result structures
for (r, t) in zip(results, titles)

    r.title = t
end

# create figure and axes
f = Figure()
ax = Axis(
    f[1, 1],
    xlabel=L"T_2 \, \textrm{(s)}", 
    xlabelsize = 20,
    xscale=log10,
    yticklabelsvisible = false,
    limits = (
        minimum(results[1].X), maximum(results[1].X), -0.1, length(results)
    ),
    yticks= 0:length(results),
)

# find the largest value in the distributions (used below to set them apart)
max_value = maximum(val for r in results for val in r.f)

# plot the results, one stacked on top of the other
for (i, r) in enumerate(results)
    lines!(
        ax, 
        r.X, 
        r.f ./ max_value .+ (i-1), 
        colormap=:tab10, colorrange=(1, 10), color=i,
    )
    text!(
        ax, 
        3e-4, (i-0.8),
        text = r.title,
        colormap=:tab10, colorrange=(1, 10), color=i,
    )
end

#save the figure as a png file, 
# with a resolution of 3 pixels per unit
save("my_results.png", f, px_per_unit=3)
```
