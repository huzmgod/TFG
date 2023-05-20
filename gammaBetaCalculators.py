from scipy.stats import gamma, beta
import numpy as np
import matplotlib.pyplot as plt

# Carga los valores de rango y sill
ranges = np.loadtxt("range_sill.dat", usecols=[0])
sills = np.loadtxt("range_sill.dat", usecols=[1])

# Ajusta los valores de rango a una distribución gamma
shape, loc, scale = gamma.fit(ranges)
ranges_gamma = gamma.pdf(ranges, shape, loc, scale)

# Ajusta los valores de sill a una distribución beta
a, b, loc, scale = beta.fit(sills)
sills_beta = beta.pdf(sills, a, b, loc, scale)

# Gráfica de los resultados
plt.hist(ranges, bins=100, density=True, alpha=0.5, label="ranges")
plt.plot(ranges, ranges_gamma, label="gamma")

plt.hist(sills, bins=100, density=True, alpha=0.5, label="sills")
plt.plot(sills, sills_beta, label="beta")

plt.legend()
plt.show()
