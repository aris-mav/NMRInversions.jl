These are some IR (inversion recovery) data points corresponding to Cheshire sandstone, for testing the invert function.

You can use them as

```julia

using NMRInversions, GLMakie

data = import_csv(IR, "path/to/data/file.csv")
results = invert(data)
plot(results)
```
