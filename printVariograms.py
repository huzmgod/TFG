import numpy as np
import matplotlib.pyplot as plt

def teorico_esferico(h, sill, rango):
    return sill * (1.5 * (h / rango) - 0.5 * (h / rango) ** 3)

def variograma_teorico(dists, sill, rango):
    return teorico_esferico(dists, sill, rango)

def ajuste_variograma(df, sill, rango):
    df["teorico"] = variograma_teorico(df["dist"], sill, rango)
    return df

# cargar valores de rango y sill
rango_sill = np.loadtxt("range_sill.dat")

n = 1000
k = 100
for j in range(n // k):
    
    plt.figure()
    for i in range(j * k, (j + 1) * k):
        # Generar una lista de distancias para el eje x
        dists = np.arange(0, 20000, 50)
        rango = rango_sill[i, 0]
        sill = rango_sill[i, 1]
        teorico = variograma_teorico(dists, sill, rango)
        plt.plot(dists, teorico, color="blue", alpha=0.1)

    plt.xlabel("Distancia (m)")
    plt.ylabel("Variograma teórico")
    plt.title("Ajuste de variograma teórico esférico")
    plt.show()


