import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker #show the numbers in x axis

mass = np.loadtxt("masse.dat")
r = mass[:,0]

plt.plot(r, mass[:,1], label = "NFW profile", color = "red")
plt.plot(r, mass[:,2], label = "NFW analytical profile", color = "green")
plt.plot(r, mass[:,3], label = "BCG Hernquist profile", color = "blue")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("r (kpc)")
plt.ylabel("M ($M_{\odot}$)")
plt.title("Mass profiles")

plt.gca().xaxis.set_major_formatter(ticker.LogFormatter()) #show the kpc scale not in scientific notation

plt.grid()
plt.legend()
plt.show()
