
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

# Leer valores de rango y sill de range_sill.dat
df = pd.read_csv('range_sill.dat', sep=' ', header=None, names=['rango', 'sill'])
#df["rango"] = pd.to_numeric(df["rango"], downcast="float")
#df["sill"] = pd.to_numeric(df["sill"], downcast="float")
# Normalizar los valores de rango y sill
r = (df['rango'] - df['rango']) / df['rango'].std()
s = (df['sill'] - df['sill'].mean()) / df['sill'].std()

# Crear un histograma para los valores de rango y sill
plt.hist(df['rango'], bins=20, edgecolor='black', alpha=0.7, label='Rango')

plt.legend()
plt.show()


plt.hist(df['sill'], bins=20, edgecolor='black', alpha=0.7, label='Sill')
plt.show()