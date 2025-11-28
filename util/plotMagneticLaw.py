import numpy as np
import matplotlib.pyplot as plt


u = np.arange(0, 1E-3, 1E-5)
a1 = 1E3
a3 = 1E5
f = a1 * u + a3 * u**3

plt.plot(u, f)
plt.show()