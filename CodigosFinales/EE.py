import MDAnalysis as mda
import numpy as np

#Centro de masa
def center_of_mass_pbc(posiciones, longitud_caja):
    "Calcula el centro geométrico"
    theta = (posiciones / longitud_caja) * 2 * np.pi
    xi_mean = np.mean(np.cos(theta))
    zeta_mean = np.mean(np.sin(theta))
    theta_mean = np.arctan2(zeta_mean, xi_mean)
    com = (theta_mean / (2 * np.pi)) * longitud_caja
    return com % longitud_caja

print("1. Cargando trayectoria...")
u = mda.Universe("trayectoria.lammpstrj", "trayectoria.lammpstrj", topology_format="LAMMPSDUMP", format="LAMMPSDUMP", in_memory=False)

# Configuración del sistema
L = 36.840314          # caja
radio_corte = 13.0     #Radio aproximado de la nanoparticula por el perfil de densidad(26 / 2)
frames_a_analizar = 50 # Promediar ultimos frames

#Crea los grupos para poder analizarlos
lipidos = u.select_atoms("type 2 3 4")
capsaicina = u.select_atoms("type 7 8 9")
total_capsaicina = len(capsaicina)

porcentajes_ee = []

print(f"2. Calculando EE% en los últimos {frames_a_analizar} frames...")

# Iterar solo con los ultimos frames
for ts in u.trajectory[-frames_a_analizar:]:
    pos_lipidos = lipidos.positions
    pos_capsaicina = capsaicina.positions
    
    # Encontrar el centro (X, Y, Z) del liposoma
    cx = center_of_mass_pbc(pos_lipidos[:, 0], L)
    cy = center_of_mass_pbc(pos_lipidos[:, 1], L)
    cz = center_of_mass_pbc(pos_lipidos[:, 2], L)
    
    # Ver que tan lejos esta cada molecula de capsaicina de ese centro, para ver si esta adentro del liposoma o no, con el radio que se aproximo antes
    dx = pos_capsaicina[:, 0] - cx
    dx = dx - L * np.round(dx / L)
    
    dy = pos_capsaicina[:, 1] - cy
    dy = dy - L * np.round(dy / L)
    
    dz = pos_capsaicina[:, 2] - cz
    dz = dz - L * np.round(dz / L)
    
    distancias = np.sqrt(dx**2 + dy**2 + dz**2)
    
    # Contar cuántas distancias son menores o iguales al radio del liposoma
    capsaicina_encapsulada = np.sum(distancias <= radio_corte)
    
    # Calcular el porcentaje en el frame
    ee_frame = (capsaicina_encapsulada / total_capsaicina) * 100
    porcentajes_ee.append(ee_frame)

#Resultados
ee_promedio = np.mean(porcentajes_ee)
ee_desviacion = np.std(porcentajes_ee)

print("Resultados eficiencia encapsulación (EE%)")
print(f"Radio límite considerado : {radio_corte} rc")
print(f"Total de 'beads' fármaco: {total_capsaicina}")
print(f"EE% Promedio             : {ee_promedio:.2f} %")
print(f"Desviación estándar (±)  : {ee_desviacion:.2f} %")
