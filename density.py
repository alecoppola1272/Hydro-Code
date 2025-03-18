import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker #show the numbers in x axis

density = np.loadtxt("density.dat")
r =density[:,0]

plt.plot(r, density[:,1], label = "Isothermal profile (NO BCG)", color = "red")
plt.plot(r, density[:,2], label = "NFW analytical profile", color = "green")
plt.plot(r, density[:,3], label = "Isothermal profile + BCG", color = "blue")
plt.plot(r, density[:,5], label = "Non-Isothermal profile + BCG", color = "black")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("r (kpc)")
plt.ylabel(r"$\rho$ (g cm$^{-3}$)") # the r before stands for "raw string", needed to correctly interpret the backslashes
plt.title("Density profiles")

plt.gca().xaxis.set_major_formatter(ticker.LogFormatter()) #show the kpc scale not in scientific notation

plt.grid()
plt.legend()
plt.show()