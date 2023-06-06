import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta, gamma, norm, kstest, chisquare
import scipy.stats as stats

"""
Carga de datos y plot de histograma
"""

df = pd.read_csv('F:/EscritorioPC/BayesianKriging/TFG/newResearch/ranges.dat', sep='\s+', names=['range'])

# Plots de rangos no normalizados
plt.hist(df['range'], bins=30, color='lightblue')
plt.xlabel('Range')
plt.ylabel('Frequency')
plt.show()

#------------------------------------------------------------

"""
Ajuste sobre una distribución normal
"""

""" # Normalización
df['range_norm'] = (df['range'] - df['range'].mean()) / df['range'].std()

# Ajusta una distribución normal
mu, std = norm.fit(df['range_norm'])

# Crea el histograma
plt.hist(df['range_norm'], bins=30, density=True, color='lightblue', alpha=0.6)

# Crea un rango de valores x para la función de densidad de probabilidad
x = np.linspace(df['range_norm'].min(), df['range_norm'].max(), 100) """
# Ajusta una distribución normal a los datos originales
mu, std = norm.fit(df['range'])

# Crea el histograma normalizado (density=True hace que la suma del área sea 1)
plt.hist(df['range'], bins=30, density=True, color='lightblue', alpha=0.6)

# Crea un rango de valores x para la función de densidad de probabilidad
x = np.linspace(df['range'].min(), df['range'].max(), 100)

p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)

plt.xlabel('Range (normalized)')
plt.ylabel('Density')

plt.show()

print('Media (mu):', mu)
print('Desviación estándar (std):', std)

"""
Residuos
"""
# Calcula los valores ajustados
df['fitted'] = norm.pdf(df['range'], mu, std)

# Calcula los residuos
df['residuals'] = df['range'] - df['fitted']

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
alpha, beta_, loc, scale = beta.fit(df['range'])
D_beta, p_beta = kstest(df['range'], 'beta', args=(alpha, beta_, loc, scale))
chisq_beta, p_chi_beta = chisquare(df['range'])

# Ajusta una distribución gamma
alpha, loc, scale = gamma.fit(df['range'])
D_gamma, p_gamma = kstest(df['range'], 'gamma', args=(alpha, loc, scale))
chisq_gamma, p_chi_gamma = chisquare(df['range'])
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


