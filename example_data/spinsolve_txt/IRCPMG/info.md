These are some IR-CPMG (T1-T2 , inversion recovery - cpmg) 
data points corresponding to Berea sandstone, 
for testing the invert function.

You can use them as

```julia

using NMRInversions, GLMakie

data = import_spinsolve(["path/to/data/T1IRT2.data","path/to/data/acqu.par"] )
results = invert(data)
plot(results)
```

