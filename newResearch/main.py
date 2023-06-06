import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import curve_fit
absPath = 'F:/EscritorioPC/BayesianKriging/TFG/'

df = pd.read_csv(absPath + 'datafiles/vel_recortada_blanked.dat', sep='\s+', names=['x', 'y', 'speed'])
df = df.dropna()

with open(f'{absPath}ranges.dat', 'w') as rangesFile:
    with open(f'{absPath}sills.dat', 'w') as sillsFile:
        for i in range(10000):
            df_sample = df.sample(n=1000)
            #print(df_sample)

            '''
            Se emplea la siguiente fórmula para calcular el variograma experimental:
            γ(h) = (1 / 2N(h)) ∑_{i=1}^{N(h)} [Z(x_i + h) - Z(x_i)]^2

            Donde:

            γ(h) es el variograma (semivariograma) a la distancia h
            Z(x_i + h) y Z(x_i) son los valores de velocidad en las ubicaciones x_i + h y x_i, respectivamente
            N(h) es el número total de pares de puntos que están separados por la distancia h
            La suma es sobre todos los pares de puntos que están a una distancia h el uno del otro
            '''
            def variograma_experimental(h, bin_width, x, y, z):
                dists = pdist(df_sample[['x', 'y']].values)

                diffs = pdist(df_sample[['speed']].values, 'sqeuclidean')

                bins = np.arange(0, h, bin_width)

                hist_with_weights, _ = np.histogram(dists, bins=bins, weights=diffs)
                hist_without_weights, _ = np.histogram(dists, bins=bins)

                with np.errstate(divide='ignore', invalid='ignore'):
                    vario = hist_with_weights / hist_without_weights
                    vario[hist_without_weights == 0] = np.nan

                bins = bins[:-1][~np.isnan(vario)]
                vario = vario[~np.isnan(vario)]

                return bins, vario / 2.0

            def variograma_esferico(h, nugget, sill, range_):
                # Para h < rango, usa el modelo esférico
                # Para h >= rango, la semivarianza es constante (= sill)
                nugget = 50
                return np.where(h < range_, nugget + sill * (1.5 * h / range_ - 0.5 * (h / range_) ** 3), nugget + sill)

            # variogram params
            h = 7500
            bin_width = 300

            # Calcula el variograma experimental
            distances, semivariances = variograma_experimental(h, bin_width, df_sample['x'], df_sample['y'], df_sample['speed'])

            # Ajusta el modelo esférico al variograma experimental
            p0 = [0.2, 0.8, 1000]
            popt, pcov = curve_fit(variograma_esferico, distances, semivariances, p0, bounds=(0, np.inf))

            nugget, sill, range_ = popt
            print(f"Nugget: {nugget}")
            print(f"Sill: {sill}")
            print(f"Range: {range_}")

            rangesFile.write(f"{range_}\n")
            sillsFile.write(f"{sill}\n")

            # Dibuja el variograma experimental y el modelo esférico
            """ plt.figure(figsize=(10, 7))
            plt.plot(distances, semivariances, marker='o', label='Experimental')
            plt.plot(distances, variograma_esferico(distances, *popt), label='Modelo esférico')
            plt.title('Variograma experimental y Modelo esférico')
            plt.xlabel('Distancia (h)')
            plt.ylabel('Semivarianza')
            plt.legend()
            plt.grid(True)
            plt.show() """
