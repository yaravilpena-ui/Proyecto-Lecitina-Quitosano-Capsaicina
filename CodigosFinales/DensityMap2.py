import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def center_of_mass_pbc(posiciones, longitud_caja):
    "Calcula el centro de masa en coordenadas circulares"
    theta = (posiciones / longitud_caja) * 2 * np.pi
    xi_mean = np.mean(np.cos(theta))
    zeta_mean = np.mean(np.sin(theta))
    theta_mean = np.arctan2(zeta_mean, xi_mean)
    com = (theta_mean / (2 * np.pi)) * longitud_caja
    return com % longitud_caja

print("Cargando topología y trayectoria...")
u = mda.Universe("trayectoria.lammpstrj", "trayectoria.lammpstrj", topology_format="LAMMPSDUMP", format="LAMMPSDUMP", in_memory=False)

lipidos = u.select_atoms("type 2 3 4")

# Longitud caja
L = 36.840314 

# El gráfico mas grande que la partícula
limite_grafico = 25 
rango_caja = [[-limite_grafico, limite_grafico], [-limite_grafico, limite_grafico]]

bins = 250
grosor_corte = 3 

densidad_total = np.zeros((bins, bins))
num_frames = len(u.trajectory)

print(f"Calculando mapa 2D en {num_frames} frames...")

for ts in u.trajectory:
    pos = lipidos.positions
    
    # Calcula centros masa 
    com_x = center_of_mass_pbc(pos[:, 0], L)
    com_y = center_of_mass_pbc(pos[:, 1], L)
    com_z = center_of_mass_pbc(pos[:, 2], L)
    
    # Hacer corte
    z_centrado = pos[:, 2] - com_z
    z_centrado = z_centrado - L * np.round(z_centrado / L)
    mascara_ecuador = (z_centrado > -grosor_corte) & (z_centrado < grosor_corte)
    
    # Extraemos X y Y 
    x_ecuador = pos[mascara_ecuador, 0]
    y_ecuador = pos[mascara_ecuador, 1]
    
    # Centramos origen (0,0)
    x_centrado = x_ecuador - com_x
    x_centrado = x_centrado - L * np.round(x_centrado / L)
    y_centrado = y_ecuador - com_y
    y_centrado = y_centrado - L * np.round(y_centrado / L)
    
    # Calcular el histograma
    densidad_frame, _, _ = np.histogram2d(x_centrado, y_centrado, bins=bins, range=rango_caja)
    densidad_total += densidad_frame

# Promedio
densidad_promedio = densidad_total / num_frames
densidad_suavizada = gaussian_filter(densidad_promedio.T, sigma=3.0) 

print("Generando imagen final...")

# Grafica
plt.figure(figsize=(8, 8)) 
plt.imshow(densidad_suavizada, origin='lower', 
           extent=[-limite_grafico, limite_grafico, -limite_grafico, limite_grafico], 
           cmap='inferno', 
           interpolation='bicubic')
barra = plt.colorbar()
barra.set_label('Densidad (partículas / unidad de área)', rotation=90, labelpad=15)
plt.xlabel('x (Å)')
plt.ylabel('y (Å)')
plt.title('Mapa de densidad 2D (coordenadas desenvolviendo)')
plt.savefig("mapa_densidad_2D_Desenrollado.png", dpi=300, bbox_inches='tight')
print("Gráfica guardada como mapa_densidad_2D.png")
