# TFG

En el fichero posiciones_tyczki.dat tenemos las posiciones de las estacasde medida de velocidad. En coordenadas UTM WGS84 33N.
Columnas:
X Y Z ID

En el fichero vel_darek.dat tenemos las velocidades horizontales diarias (en m/día) de las estacas en el periodo desde el 1 de mayo de 2005 hasta el 30 de abril de 2011. 
Columnas:
1 2 3 4 5 6 7 8 9 10 11 t t t t t Número de día Día Mes Año

En el fichero vel_recortada.dat tenemos el campo de velocidad horizontal obtenido como media de las velocidades medidas en el periodo de observación, de enero de 2013 a agosto de 2014. En coordenadas UTM WGS84 33N
Columnas:
X Y Vel(m/día) Vel(m/año)

En el fichero contorno_hansbreen.dat tiene los puntos del contorno de Hansbreen. En coordenadas UTM WGS84 33N
Columnas:
ID X Y Las demás no son interesantes para ti.

En el fichero cuadricula_superficie_Hansbreen.dat se tiene la matriz rectangular Z(y,x) de alturas del glaciar, con filas correspondientes a valores de 'y', en orden 
descendente, y columnas correspondientes a valores de 'x', en orden ascendente. Es decir, la matriz dibuja la forma del glaciar de norte a sur y de este a oeste, tal como se vería en un mapa.

En el fichero vel_interpolada_DEM_ordenada.dat se establecen los valores de las velocidades interpoladas del fichero vel_recortada.dat sobre el DEM en el que se está trabajando.
COLUMNAS: 
X   Y   PUNTOS_DE_INTERPOLACIÓN   DISTANCIAS_A_PUNTOS   VELOCIDADES_PUNTOS_INTERPOLACIÓN    VELOCIDAD_ESTIMADA

En los ficheros gradient_x_dem.dat y gradient_y_dem.dat se encuentran los valores del gradiente de la matriz Z del DEM (se carga desde el fichero cuadricula_superficie_Hansbreen.dat). Filas según valor de 'y' en orden descendente y columnas según valor de 'x' en orden ascendente. 

En el fichero gradient_in_stick.dat se calcula media ponderada por el inverso de la distancia de las componentes de gradiente, calculadas para el DEM, de los 3 puntos del DEM que más cerca están de la posición de la estaca. 
COLUMNAS:
X   Y   PUNTOS_DE_INTERPOLACIÓN   DISTANCIAS_A_PUNTOS   COMPONENTE_X_GRADIENTE    COMPONENTE_Y_GRADIENTE

En el fichero speed_components se calcula, para cada coordenada del DEM, las componentes de la velocidad (anuales) en función de los valores del gradiente calculado. 
Las dos columnas intermedias son auxiliares, no son de interés.
COLUMNAS:
X   Y   X_MATRIZ    Y_MATRIZ    VELOCIDAD_X     VELOCIDAD_Y    

En el fichero vel_darek_mensuales.dat se calcula la suma por cada 30 días de las velocidades diarias de las estacas.
COLUMNAS: 
1 2 3 4 5 6 7 8 9 10 11 t t t t t

En el fichero updated_stick_positions.dat se obtiene la posición de las estacas actualizada por meses desde mayo de 2005 hasta abril de 2011.
COLUMNAS:
1 2 3 4 5 6 7 8 9 10 11 t t t t t
