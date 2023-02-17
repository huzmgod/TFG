import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

# Leer valores de rango y sill de range_sill.dat
df = pd.read_csv('range_sill.dat', sep='\t', header=None, names=['rango', 'sill'])


# Crear un histograma para los valores de rango y sill
plt.hist(df['rango'], bins=20, edgecolor='black', alpha=0.7, label='Rango')
plt.hist(df['sill'], bins=20, edgecolor='black', alpha=0.7, label='Sill')
plt.legend()
plt.show()

# Ajustar una distribución normal a los valores de rango y sill
rango_fit = stats.norm.fit(df['rango'])
sill_fit = stats.norm.fit(df['sill'])

# Crear un histograma para los valores de rango y sill ajustados a una distribución normal
plt.hist(df['rango'], bins=20, edgecolor='black', alpha=0.7, label='Rango', density=True)
plt.hist(df['sill'], bins=20, edgecolor='black', alpha=0.7, label='Sill', density=True)
x = np.linspace(df['rango'].min(), df['rango'].max(), 100)
y = stats.norm.pdf(x, loc=rango_fit[0], scale=rango_fit[1])
plt.plot(x, y, label='Distribución normal')
plt.legend()
plt.show()

# Ajustar una chi cuadrado a los valores de rango y sill
rango_fit = stats.chi2.fit(df['rango'])
sill_fit = stats.chi2.fit(df['sill'])

# Crear un histograma para los valores de rango y sill ajustados a una chi cuadrado
plt.hist(df['rango'], bins=20, edgecolor='black', alpha=0.7, label='Rango', density=True)
plt.hist(df['sill'], bins=20, edgecolor='black', alpha=0.7, label='Sill', density=True)
x = np.linspace(df['rango'].min(), df['rango'].max(), 100)
y = stats.chi2.pdf(x, df=rango_fit[0], scale=rango_fit[1])
plt.plot(x, y, label='Chi cuadrado')
