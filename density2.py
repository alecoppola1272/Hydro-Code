import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker #show the numbers in x axis

density = np.loadtxt("density.dat")
r =density[:,0]

plt.plot(r, density[:,5], label = "Non-Isothermal profile + BCG", color = "black")
plt.plot(r, density[:,6], label = "Rebusco profile", color = "purple")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("r (kpc)")
plt.ylabel(r"$\rho$ (g cm$^{-3}$)") # the r before stands for "raw string", needed to correctly interpret the backslashes
plt.title("Rebusco and non-isothermal profiles")

plt.gca().xaxis.set_major_formatter(ticker.LogFormatter()) #show the kpc scale not in scientific notation

plt.grid()
plt.legend()
plt.show()