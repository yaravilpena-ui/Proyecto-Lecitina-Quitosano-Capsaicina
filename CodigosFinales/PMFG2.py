import MDAnalysis as mda
from MDAnalysis.analysis import rdf
import numpy as np
import matplotlib.pyplot as plt
import warnings

print("1. Cargando topología y trayectoria...")
u = mda.Universe("trayectoriaG2.lammpstrj", "trayectoriaG2.lammpstrj", 
                 topology_format="LAMMPSDUMP", format="LAMMPSDUMP", in_memory=False)

print("2. Definiendo los grupos de moléculas...")
# Crear grupos de atomos para mapear
lecitina = u.select_atoms("type 2 3 4")
quitosano = u.select_atoms("type 5 6")
capsaicina = u.select_atoms("type 7 8 9")

# Rango de distancia para el análisis
distancia_max = 12.0
bines = 150

print("3. Calculando Funciones de Distribución Radial g(r)...")
# (A) Lecitina - Quitosano (Lecithin-CS)
rdf_A = rdf.InterRDF(lecitina, quitosano, nbins=bines, range=(0.0, distancia_max))

# (B) Lecitina - Capsaicina (Lecithin-Capsaicin)
rdf_B = rdf.InterRDF(lecitina, capsaicina, nbins=bines, range=(0.0, distancia_max))

# (C) Quitosano - Capsaicina (CS-Capsaicin)
rdf_C = rdf.InterRDF(quitosano, capsaicina, nbins=bines, range=(0.0, distancia_max))

# Analiza solo los últimos frames, para que sea más rápido
total_frames = len(u.trajectory)
inicio_equilibrio = max(0, total_frames - 100) #por si hay menos de 100

print(f"Calculando g(r) solo para los últimos frames (desde {inicio_equilibrio} hasta {total_frames})...")

rdf_A.run(start=inicio_equilibrio)
rdf_B.run(start=inicio_equilibrio)
rdf_C.run(start=inicio_equilibrio)
# ------------------------------------------------------------------------

def calcular_pmf(g_r):
    "Convierte el g(r) en PMF aplicando el logaritmo natural"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore") # Ignorar advertencias por log(0)
        pmf = -np.log(g_r)
    
    "Cuando g(r) es 0 (las moléculas no se tocan), el PMF es infinito"
    pmf[np.isinf(pmf)] = 5.0
    pmf[np.isnan(pmf)] = 5.0
    return pmf

print("4. Convirtiendo g(r) a Potencial de Fuerza Media (PMF)...")
pmf_A = calcular_pmf(rdf_A.results.rdf)
pmf_B = calcular_pmf(rdf_B.results.rdf)
pmf_C = calcular_pmf(rdf_C.results.rdf)
distancias = rdf_A.results.bins

print("5. Generando gráficas de resultados...")
plt.figure(figsize=(10, 6))

# Trazar las 3 curvas
plt.plot(distancias, pmf_A, label='(A) Lecitina - Quitosano', color='blue', linewidth=2)
plt.plot(distancias, pmf_B, label='(B) Lecitina - Capsaicina', color='red', linewidth=2)
plt.plot(distancias, pmf_C, label='(C) Quitosano - Capsaicina', color='green', linewidth=2)

# Configuración visual de la gráfica
plt.axhline(0, color='black', linestyle='--', linewidth=1) # Línea de energía cero
plt.ylim(-3.0, 5.0) # Ajusta estos límites según tus resultados
plt.xlim(0, distancia_max)

plt.title('Potencial de Fuerza Media (PMF)', fontsize=14, fontweight='bold')
plt.xlabel(r'Distancia radial $r$ (unidades $r_c$)', fontsize=12)
plt.ylabel(r'PMF ($k_B T$)', fontsize=12)
plt.legend(loc='upper right', fontsize=10)
plt.grid(True, linestyle=':', alpha=0.7)

# Guardar la imagen
plt.savefig("Grafica_PMF_Resultados.png", dpi=300, bbox_inches='tight')
print("La gráfica ha sido guardada como Grafica_PMF_Resultados.png")
