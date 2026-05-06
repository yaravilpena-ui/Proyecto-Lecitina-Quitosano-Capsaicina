import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

print("Cargando topología y trayectoria...")
u = mda.Universe("trayectoria.lammpstrj", "trayectoria.lammpstrj", topology_format="LAMMPSDUMP", format="LAMMPSDUMP", in_memory=False)

lipidos = u.select_atoms("type 2 3 4")

# Longitud de la caja
L = 36.840314
limite = L / 2
rango_caja = [-limite, limite]
bins = 200

perfil_x = np.zeros(bins)
perfil_y = np.zeros(bins)
perfil_z = np.zeros(bins)

num_frames = len(u.trajectory)
print(f"Desenrollando y calculando perfiles en {num_frames} frames...")

def center_of_mass_pbc(posiciones, longitud_caja):
    "Calcula el centro de masa"
    # Convertimos las posiciones lineales en ángulos (de 0 a 360 grados / 2*Pi)
    theta = (posiciones / longitud_caja) * 2 * np.pi
    
    # Calculamos el promedio circular usando senos y cosenos
    xi_mean = np.mean(np.cos(theta))
    zeta_mean = np.mean(np.sin(theta))
    theta_mean = np.arctan2(zeta_mean, xi_mean)
    
    #  Regresamos el ángulo a coordenadas lineales
    com = (theta_mean / (2 * np.pi)) * longitud_caja
    return com % longitud_caja

for ts in u.trajectory:
    pos = lipidos.positions
    
    # Calculamos el centro para X, Y y Z
    com_x = center_of_mass_pbc(pos[:, 0], L)
    com_y = center_of_mass_pbc(pos[:, 1], L)
    com_z = center_of_mass_pbc(pos[:, 2], L)
    
    # Centramos las posiciones restando centro de masa
    x_centrado = pos[:, 0] - com_x
    y_centrado = pos[:, 1] - com_y
    z_centrado = pos[:, 2] - com_z
    
    # Todo en la caja
    x_centrado = x_centrado - L * np.round(x_centrado / L)
    y_centrado = y_centrado - L * np.round(y_centrado / L)
    z_centrado = z_centrado - L * np.round(z_centrado / L)
    
    # Calculamos los histogramas con las coordenadas curadas
    hx, bordes = np.histogram(x_centrado, bins=bins, range=rango_caja)
    hy, _ = np.histogram(y_centrado, bins=bins, range=rango_caja)
    hz, _ = np.histogram(z_centrado, bins=bins, range=rango_caja)
    
    perfil_x += hx
    perfil_y += hy
    perfil_z += hz

# Promedios
perfil_x /= num_frames
perfil_y /= num_frames
perfil_z /= num_frames

# Normalizar
perfil_x = perfil_x / np.max(perfil_x)
perfil_y = perfil_y / np.max(perfil_y)
perfil_z = perfil_z / np.max(perfil_z)
eje_posiciones = ((bordes[:-1] + bordes[1:]) / 2) + limite

# Gráfica
fig, axs = plt.subplots(1, 3, figsize=(15, 4)) 
limite_superior_y = 1.5
limite_superior_x = L 

# Eje X
axs[0].plot(eje_posiciones, perfil_x, color='red', linewidth=1.5)
axs[0].set_xlabel('x*')
axs[0].set_ylabel('ρ*(x)')
axs[0].set_xlim(0, limite_superior_x)
axs[0].set_ylim(0, limite_superior_y)

# Eje Y
axs[1].plot(eje_posiciones, perfil_y, color='lime', linewidth=1.5)
axs[1].set_xlabel('y*')
axs[1].set_ylabel('ρ*(y)')
axs[1].set_xlim(0, limite_superior_x)
axs[1].set_ylim(0, limite_superior_y)

# Eje Z
axs[2].plot(eje_posiciones, perfil_z, color='blue', linewidth=1.5)
axs[2].set_xlabel('z*')
axs[2].set_ylabel('ρ*(z)')
axs[2].set_xlim(0, limite_superior_x)
axs[2].set_ylim(0, limite_superior_y)

plt.tight_layout()
plt.savefig("perfiles_densidad_1D_Curados.png", dpi=300)
print("Gráfico guardado como perfiles_densidad_1D_Curados.png")
