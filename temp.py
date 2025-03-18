import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker #show the numbers in x axis

temp = np.loadtxt("temperature.dat")
r =temp[:,0]

plt.plot(r, temp[:,1], label = "Non-Isothermal profile", color = "black")
plt.plot(r, temp[:,2], label = "Isothermal profile", color = "red")
plt.plot(r, temp[:,3], label = "Rebusco profile", color = "blue")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("r (kpc)")
plt.ylabel("T (K)")
plt.title("Temperature profiles")

plt.gca().xaxis.set_major_formatter(ticker.LogFormatter()) #show the kpc scale not in scientific notation

plt.legend()
plt.grid()
plt.show()