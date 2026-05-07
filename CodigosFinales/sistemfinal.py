import math
import random

# --- 1. CONFIGURACIÓN DEL SISTEMA ---
archivo_entrada = "final.liposome.data"
archivo_salida = "Sistema_Completo.data"

# Margen de seguridad: La caja de LAMMPS es de 18.42. 
# Generamos átomos hasta 17.0 para evitar que los enlaces crucen la frontera.
limite_caja_seguro = 17.0       

# Aumentamos el radio de exclusión ligeramente para evitar solapamientos 
# con los grupos fosfato de la lecitina.
radio_exclusion = 16.0   

num_capsaicina = 250
num_quitosano = 50

# Tipos: 1,2,3=Lecitina | 4=Agua | 5,6=Quitosano | 7,8,9=Capsaicina
secuencia_quitosano = [5, 6, 6, 6, 5, 6, 6, 6, 5, 6] * 5  

# --- 2. LECTURA EN MEMORIA POR BLOQUES ---
with open(archivo_entrada, "r") as f:
    lineas = f.readlines()

secciones = {}
seccion_actual = "Header"
secciones[seccion_actual] = []

for linea in lineas:
    if linea.startswith("Masses"): seccion_actual = "Masses"
    elif linea.startswith("Pair Coeffs"): seccion_actual = "Pair Coeffs"
    elif linea.startswith("Bond Coeffs"): seccion_actual = "Bond Coeffs"
    elif linea.startswith("Angle Coeffs"): seccion_actual = "Angle Coeffs"
    elif linea.startswith("Atoms"): seccion_actual = "Atoms"
    elif linea.startswith("Velocities"): seccion_actual = "Velocities"
    elif linea.startswith("Bonds"): seccion_actual = "Bonds"
    elif linea.startswith("Angles"): seccion_actual = "Angles"
    elif linea.startswith("Dihedrals"): seccion_actual = "Dihedrals"
    
    if seccion_actual not in secciones:
        secciones[seccion_actual] = []
    secciones[seccion_actual].append(linea)

# Extraer conteos originales del encabezado
old_atoms, old_bonds, old_angles = 0, 0, 0
for linea in secciones["Header"]:
    if " atoms" in linea: old_atoms = int(linea.split()[0])
    elif " bonds" in linea: old_bonds = int(linea.split()[0])
    elif " angles" in linea: old_angles = int(linea.split()[0])

atom_id = old_atoms + 1
mol_id = 200000  
bond_id = old_bonds + 1
angle_id = old_angles + 1

# --- 3. FUNCIONES ESPACIALES ESTRICTAS ---
def es_valida(x, y, z):
    """Verifica que un punto específico no rompa las fronteras ni penetre el liposoma."""
    if abs(x) > limite_caja_seguro or abs(y) > limite_caja_seguro or abs(z) > limite_caja_seguro:
        return False
    if math.sqrt(x**2 + y**2 + z**2) <= radio_exclusion:
        return False
    return True

def coord_segura():
    """Genera la coordenada para el átomo cabeza (semilla)."""
    while True:
        x = random.uniform(-limite_caja_seguro, limite_caja_seguro)
        y = random.uniform(-limite_caja_seguro, limite_caja_seguro)
        z = random.uniform(-limite_caja_seguro, limite_caja_seguro)
        if es_valida(x, y, z):
            return x, y, z

def vector_enlace(longitud=0.7):
    theta = random.uniform(0, 2*math.pi)
    phi = math.acos(random.uniform(-1, 1))
    return longitud * math.sin(phi) * math.cos(theta), longitud * math.sin(phi) * math.sin(theta), longitud * math.cos(phi)

# --- 4. CONSTRUCCIÓN DE NUEVAS MOLÉCULAS ---
str_atomos, str_vels, str_bonds, str_angles = [], [], [], []

print("Generando moléculas de Capsaicina (Topología 4 beads: C1-C2-C3-C3)...")
for _ in range(num_capsaicina):
    valido = False
    while not valido:
        # Generación de la coordenada base (C1)
        c_x, c_y, c_z = coord_segura()
        
        # 3 Vectores para los 3 enlaces (C1-C2, C2-C3, C3-C3)
        d1x, d1y, d1z = vector_enlace()
        d2x, d2y, d2z = vector_enlace()
        d3x, d3y, d3z = vector_enlace()
        
        # Proyección de los beads
        p2x, p2y, p2z = c_x+d1x, c_y+d1y, c_z+d1z          # C2
        p3x, p3y, p3z = p2x+d2x, p2y+d2y, p2z+d2z          # C3 (Primero)
        p4x, p4y, p4z = p3x+d3x, p3y+d3y, p3z+d3z          # C3 (Segundo)
        
        # Validación estricta para que toda la molécula esté en el espacio permitido
        if es_valida(p2x, p2y, p2z) and es_valida(p3x, p3y, p3z) and es_valida(p4x, p4y, p4z):
            valido = True
            
    # Asignación de 4 IDs consecutivos
    a1, a2, a3, a4 = atom_id, atom_id+1, atom_id+2, atom_id+3
    atom_id += 4
    
    # Escritura: C1(7), C2(8), C3(9), C3(9)
    str_atomos.append(f"{a1} {mol_id} 7 0.0 {c_x:.4f} {c_y:.4f} {c_z:.4f} 0 0 0\n")
    str_atomos.append(f"{a2} {mol_id} 8 0.0 {p2x:.4f} {p2y:.4f} {p2z:.4f} 0 0 0\n")
    str_atomos.append(f"{a3} {mol_id} 9 0.0 {p3x:.4f} {p3y:.4f} {p3z:.4f} 0 0 0\n")
    str_atomos.append(f"{a4} {mol_id} 9 0.0 {p4x:.4f} {p4y:.4f} {p4z:.4f} 0 0 0\n")
    
    # Velocidades iniciales
    str_vels.extend([
        f"{a1} 0.0 0.0 0.0\n", 
        f"{a2} 0.0 0.0 0.0\n", 
        f"{a3} 0.0 0.0 0.0\n", 
        f"{a4} 0.0 0.0 0.0\n"
    ])
    
    # Conectividad: 3 Enlaces
    str_bonds.append(f"{bond_id} 1 {a1} {a2}\n"); bond_id += 1 # C1 - C2
    str_bonds.append(f"{bond_id} 1 {a2} {a3}\n"); bond_id += 1 # C2 - C3
    str_bonds.append(f"{bond_id} 1 {a3} {a4}\n"); bond_id += 1 # C3 - C3
    
    # Conectividad: 2 Ángulos
    str_angles.append(f"{angle_id} 3 {a1} {a2} {a3}\n"); angle_id += 1 # C1 - C2 - C3
    str_angles.append(f"{angle_id} 3 {a2} {a3} {a4}\n"); angle_id += 1 # C2 - C3 - C3
    
    mol_id += 1

print("Generando cadenas de Quitosano...")
for _ in range(num_quitosano):
    cadena_valida = False
    while not cadena_valida:
        coords = []
        x, y, z = coord_segura()
        coords.append((x, y, z))
        exito = True
        
        for _ in range(1, len(secuencia_quitosano)):
            intentos_paso = 0
            paso_exito = False
            while intentos_paso < 50:
                dx, dy, dz = vector_enlace()
                nx, ny, nz = x + dx, y + dy, z + dz
                if es_valida(nx, ny, nz):
                    x, y, z = nx, ny, nz
                    coords.append((x, y, z))
                    paso_exito = True
                    break
                intentos_paso += 1
            
            if not paso_exito:
                exito = False
                break 
        
        if exito:
            cadena_valida = True
            
    ids_cadena = []
    for i, (cx, cy, cz) in enumerate(coords):
        tipo = secuencia_quitosano[i]
        str_atomos.append(f"{atom_id} {mol_id} {tipo} 0.0 {cx:.4f} {cy:.4f} {cz:.4f} 0 0 0\n")
        str_vels.append(f"{atom_id} 0.0 0.0 0.0\n")
        ids_cadena.append(atom_id)
        atom_id += 1
        
    for i in range(len(ids_cadena) - 1):
        str_bonds.append(f"{bond_id} 1 {ids_cadena[i]} {ids_cadena[i+1]}\n")
        bond_id += 1
    for i in range(len(ids_cadena) - 2):
        str_angles.append(f"{angle_id} 2 {ids_cadena[i]} {ids_cadena[i+1]} {ids_cadena[i+2]}\n")
        angle_id += 1
    mol_id += 1

# --- 5. PURGA DE AGUA Y CONECTIVIDAD HUÉRFANA ---
print("Purgando agua y ajustando topología...")
agua_a_eliminar = len(str_atomos)
agua_eliminada = 0
atomos_eliminados = set()

# Purgar Atoms
nuevos_atoms = []
for linea in secciones.get("Atoms", []):
    partes = linea.split()
    if len(partes) >= 7 and partes[0].isdigit() and partes[2] == "4" and agua_eliminada < agua_a_eliminar:
        agua_eliminada += 1
        atomos_eliminados.add(partes[0])
    else:
        nuevos_atoms.append(linea)

idx = len(nuevos_atoms)
while idx > 0 and nuevos_atoms[idx-1].strip() == "": idx -= 1
secciones["Atoms"] = nuevos_atoms[:idx] + str_atomos + nuevos_atoms[idx:]

# Purgar Velocidades
if "Velocities" in secciones:
    nuevos_vels = []
    for linea in secciones["Velocities"]:
        partes = linea.split()
        if len(partes) >= 3 and partes[0].isdigit() and partes[0] in atomos_eliminados:
            continue
        nuevos_vels.append(linea)
    
    idx = len(nuevos_vels)
    while idx > 0 and nuevos_vels[idx-1].strip() == "": idx -= 1
    secciones["Velocities"] = nuevos_vels[:idx] + str_vels + nuevos_vels[idx:]

bonds_eliminados = 0
angles_eliminados = 0

# Purgar Bonds huérfanos
if "Bonds" in secciones:
    nuevos_bonds = []
    for linea in secciones["Bonds"]:
        partes = linea.split()
        if len(partes) >= 4 and partes[0].isdigit():
            if partes[2] in atomos_eliminados or partes[3] in atomos_eliminados:
                bonds_eliminados += 1
                continue
        nuevos_bonds.append(linea)
        
    idx = len(nuevos_bonds)
    while idx > 0 and nuevos_bonds[idx-1].strip() == "": idx -= 1
    secciones["Bonds"] = nuevos_bonds[:idx] + str_bonds + nuevos_bonds[idx:]

# Purgar Angles huérfanos
if "Angles" in secciones:
    nuevos_angles = []
    for linea in secciones["Angles"]:
        partes = linea.split()
        if len(partes) >= 5 and partes[0].isdigit():
            if partes[2] in atomos_eliminados or partes[3] in atomos_eliminados or partes[4] in atomos_eliminados:
                angles_eliminados += 1
                continue
        nuevos_angles.append(linea)
        
    idx = len(nuevos_angles)
    while idx > 0 and nuevos_angles[idx-1].strip() == "": idx -= 1
    secciones["Angles"] = nuevos_angles[:idx] + str_angles + nuevos_angles[idx:]

# --- 6. REESCRITURA TOTAL DE COEFICIENTES ---
if "Masses" in secciones:
    secciones["Masses"] = ["Masses\n\n"] + [f"{j} 1.0\n" for j in range(1, 10)] + ["\n"]
if "Pair Coeffs" in secciones:
    secciones["Pair Coeffs"] = ["Pair Coeffs\n\n"] + [f"{j} 78.33 4.5\n" for j in range(1, 10)] + ["\n"]
if "Angle Coeffs" in secciones:
    secciones["Angle Coeffs"] = ["Angle Coeffs\n\n", "1 25 170\n", "2 25 170\n", "3 25 170\n", "\n"]

# --- 7. ACTUALIZACIÓN DINÁMICA DEL HEADER ---
for i, linea in enumerate(secciones["Header"]):
    if " atoms" in linea:
        secciones["Header"][i] = f"{old_atoms - agua_eliminada + len(str_atomos)} atoms\n"
    elif " bonds" in linea:
        secciones["Header"][i] = f"{old_bonds - bonds_eliminados + len(str_bonds)} bonds\n"
    elif " angles" in linea:
        secciones["Header"][i] = f"{old_angles - angles_eliminados + len(str_angles)} angles\n"
    elif " atom types" in linea:
        secciones["Header"][i] = "9 atom types\n"
    elif " angle types" in linea:
        secciones["Header"][i] = "3 angle types\n"

# --- 8. ESCRITURA FINAL ---
orden_secciones = ["Header", "Masses", "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", 
                   "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals"]

with open(archivo_salida, "w") as f_out:
    for nombre in orden_secciones:
        if nombre in secciones:
            f_out.write("".join(secciones[nombre]))

print(f"Éxito total: Se insertaron {num_capsaicina*4 + num_quitosano*50} nuevos átomos (Capsaicina ajustada a 4 beads).")
print(f"Se purgaron {agua_eliminada} átomos de agua, {bonds_eliminados} enlaces y {angles_eliminados} ángulos.")
