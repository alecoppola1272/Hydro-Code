import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker #show the numbers in x axis

fe= np.loadtxt("zfe_initial.dat")
r =fe[:,0]

plt.plot(r, fe[:,1], label = "Fe abundance - $Z_{Fe,out} $ ", color = "red")
plt.plot(r, fe[:,2], label = "Fe abudance", color = "blue")

plt.xscale("log")
#plt.yscale("log")

plt.xlabel("r (kpc)")
plt.ylabel(r"$Z_{Fe}$ ($Z_{Fe,\odot}$)") # the r before stands for "raw string", needed to correctly interpret the backslashes
plt.title("Iron profiles")

plt.gca().xaxis.set_major_formatter(ticker.LogFormatter()) #show the kpc scale not in scientific notation

plt.grid()
plt.legend()
plt.show()