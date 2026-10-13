# Using NMRInversions.jl from Python

To get started, install
[juliacall](https://pypi.org/project/juliacall/)
using your preferred method (`pip`, `uv`, etc.) in
a Python [virtual
environment](https://docs.python.org/3/tutorial/venv.html).
Afterwards, execute the following from `Python`:

```python
import juliapkg as pkg
pkg.add("NMRInversions","55c20db2-0166-4687-95c3-62a9c7afb29b")
pkg.add("GLMakie","e9467ef8-e4e7-5192-8a1a-b1aee30e663a") #optional, for plots
pkg.resolve()
```

The code above installs Julia (if it is not
already available in your system), as well as the
`NMRInversions` package. This can take several
minutes the first time.

Afterwards, you can use the package as follows:

```python
from juliacall import Main as jl
jl.seval("using NMRInversions")
jl.seval("using GLMakie") # optional
```

Now you can call any Julia function or object from
python, as long as you prefix it with `jl.`.

Here we [import](functions.md#Importing-data) some
[CPMG](types_structs.md#Pulse-Sequences-Data-Axes)
data from a `.csv` file,
[invert](functions.md#Inversion), and plot the
results using `GLMakie`:

```python
...

data = jl.import_csv(jl.CPMG)
results = jl.invert(data)

fig = jl.plot(results)
screen = jl.display(fig)
jl.wait(screen)  # keeps the window alive until you close it
```

Of course, you can use whatever package you like
to plot the data:

```python
...

import matplotlib.pyplot as plt

plt.semilogx(results.dimensions[0],
results.distribution)
plt.xlabel(r"$T_2$ (s)")
plt.show()
```

If you have your NMR data already imported in
Python (e.g., in the form of `numpy` arrays, `t`
for time and `s` for signal) you can still invert
them by creating a `data` structure as follows:

```python
...

data = jl.ExperimentData(jl.CPMG(t), s) 
results = jl.invert(data)
```

Similarly, for 2D cases, suppose you have a T1-T2
experiment with CPMG as the direct dimension and
inversion recovery as the indirect one. You can
do:

```python
...

data = jl.ExperimentData(
    (jl.CPMG(t_direct), jl.IR(t_indirect)),
    data_matrix,
)
results = jl.invert(data)
```

Note that Julia structures are not automatically
converted to Python equivalents, but doing so is
straightforward (if needed):

```python
...

import numpy as np

d = np.array(results.distribution)
```

The package(s) installed here via Julia can be
updated to the latest version by:

```python
import juliapkg as pkg
pkg.update()  
```
