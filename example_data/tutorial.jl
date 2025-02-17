using NMRInversions, GLMakie


# 1D example

data = import_csv(IR) # select the csv IR file in the dialog

a = expfit(1, data)  # mono-exponential fit
b = expfit(2, data)  # bi-exponential fit

plot(a,b)  # Visualize both on the same plot

a.eqn  # Print the equation of the mono-exponential fit
b.eqn  # Print the equation of the bi-exponential fit

results = invert(data) # do an inversion
plot(results) # visualize and interact with results

results = invert(data, alpha = lcurve(0.001, 10)) #check to see whether lcurve produces different results from gcv
plot(results) # visualize and interact with results



# 2D example

data = import_spinsolve() # select the spinsolve IRCPMG file in the dialog
results = invert(data)
plot(results)
