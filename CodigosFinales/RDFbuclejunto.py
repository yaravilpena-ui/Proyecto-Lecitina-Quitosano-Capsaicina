import numpy as np
import matplotlib.pyplot as plt
import freud
import gsd.hoomd

# Configuración
bins = 600
rmax = 18 #por el tamaño de la caja calculado con la densidad
num_archivos = 11
archivos = [f'RDFTipo1.{i}' for i in range(num_archivos)]  # nombres de archivos exportados de OVITO

# Lista para guardar los resultados
rdfs = []
bin_centers = None

# Procesar cada archivo
for nombre in archivos:
    try:
        traj = gsd.hoomd.open(nombre, 'r')
    except FileNotFoundError:
        print(f"Archivo {nombre} no encontrado, omitiendo...")
        continue

    rdf = freud.density.RDF(bins=bins, r_max=rmax)
    for frame in traj:
        rdf.compute(system=frame, reset=False)
    rdfs.append(rdf.rdf)
    if bin_centers is None:
        bin_centers = rdf.bin_centers
    traj.close()

# Graficar todas las RDF en la misma figura
plt.figure(figsize=(10, 6))
for i, rdf_vals in enumerate(rdfs):
    plt.plot(bin_centers, rdf_vals, label=f'Archivo {archivos[i]}')

plt.title("Radial Distribution Function - Comparación")
plt.xlabel("$r$")
plt.ylabel("$g(r)$")
plt.legend()
plt.grid(True)
plt.savefig("RDFTipo1_comparativa.png", dpi=300)
plt.show()
