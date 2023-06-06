import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, gamma, norm, kstest, chisquare
import scipy.stats as stats

"""
Carga de datos y plot de histograma
"""

df = pd.read_csv('F:/EscritorioPC/BayesianKriging/TFG/newResearch/sills.dat', sep='\s+', names=['sill'])

# Plots de rangos no normalizados
""" plt.hist(df['sill'], bins=30, color='lightblue')
plt.xlabel('sill')
plt.ylabel('Frequency')
plt.show() """

#------------------------------------------------------------

"""
Ajuste sobre una distribución normal
"""

# Normalización
df['sill_norm'] = (df['sill'] - df['sill'].mean()) / df['sill'].std()

# Ajusta una distribución normal
mu, std = norm.fit(df['sill'])
D_norm, p_norm = kstest(df['sill'], 'norm', args=(mu, std))
chisq_norm, p_chi_norm = chisquare(df['sill'])

# Ajusta una distribución normal
mu, std = norm.fit(df['sill_norm'])

# Crea el histograma
plt.hist(df['sill_norm'], bins=30, density=True, color='lightblue', alpha=0.6)

# Crea un rango de valores x para la función de densidad de probabilidad
x = np.linspace(df['sill_norm'].min(), df['sill_norm'].max(), 100)

# Calcula y dibuja la función de densidad de probabilidad
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

# Etiqueta los ejes
plt.xlabel('sill (normalized)')
plt.ylabel('Density')

plt.show()

print('Media (mu):', mu)
print('Desviación estándar (std):', std)

"""
Residuos
"""
# Calcula los valores ajustados
df['fitted'] = norm.pdf(df['sill_norm'], mu, std)

# Calcula los residuos
df['residuals'] = df['sill_norm'] - df['fitted']

# Crea un histograma de los residuos
plt.hist(df['residuals'], bins=30, density=True, color='lightblue', alpha=0.6)
plt.xlabel('Residuos')
plt.ylabel('Densidad')
plt.show()

# Realiza una prueba de normalidad en los residuos
k2, p = stats.normaltest(df['residuals'])
print("Normality test p-value:", p)

# QQ Graph
stats.probplot(df['residuals'], plot=plt)
plt.show()

#------------------------------------------------------------

""" # Ajusta una distribución beta
alpha, beta_, loc, scale = beta.fit(df['sill'])
D_beta, p_beta = kstest(df['sill'], 'beta', args=(alpha, beta_, loc, scale))
chisq_beta, p_chi_beta = chisquare(df['sill'])

# Ajusta una distribución gamma
alpha, loc, scale = gamma.fit(df['sill'])
D_gamma, p_gamma = kstest(df['sill'], 'gamma', args=(alpha, loc, scale))
chisq_gamma, p_chi_gamma = chisquare(df['sill'])
 """

# Tests de bondad
""" print('Kolmogorov-Smirnov test:')
print('Beta: D =', D_beta, ', p =', p_beta)
print('Gamma: D =', D_gamma, ', p =', p_gamma)
print('Normal: D =', D_norm, ', p =', p_norm)

print('Chi-square test:')
print('Beta: χ2 =', chisq_beta, ', p =', p_chi_beta)
print('Gamma: χ2 =', chisq_gamma, ', p =', p_chi_gamma)
print('Normal: χ2 =', chisq_norm, ', p =', p_chi_norm) """


