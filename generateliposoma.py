import math
import random   # Necesario para generar posiciones aleatorias del agua

# ------------------------------------------------------------
# 1. PARÁMETROS DEL SISTEMA (basados en el artículo)
# ------------------------------------------------------------

# Número total de moléculas de lecitina en el liposoma.
# 4929 corresponde a una concentración de 0.60 M (ver artículo, Tabla 1 y Figura S2)
total_moleculas = 4929

# Mitad de las moléculas irán en la capa externa y la otra mitad en la capa interna.
# La capa externa tiene una molécula más (para números impares) y así se reparte mejor.
mol_ext = total_moleculas // 2 + 1   # 2465 moléculas externas
mol_int = total_moleculas // 2       # 2464 moléculas internas

# ------------------------------------------------------------
# Cálculos para la cabecera del archivo de datos de LAMMPS
# ------------------------------------------------------------
# Cada molécula de lecitina tiene 3 átomos (L1, L2, L3)
btot = total_moleculas * 3               # número total de átomos de lecitina
num_agua = 192000 - total_moleculas * 3  # calculado para densidad global = 3
total_atoms = btot + num_agua            # átomos totales (lecitina + agua)

bonds = total_moleculas * 2              # cada molécula tiene 2 enlaces (L1-L2 y L2-L3)
angles = total_moleculas * 1             # cada molécula tiene 1 ángulo (L1-L2-L3)

atypes = 4        # tipos de átomo: 1=L1, 2=L2, 3=L3, 4=agua (W)
btypes = 1        # un solo tipo de enlace (armónico)
antypes = 1       # un solo tipo de ángulo (para la lecitina)

# ------------------------------------------------------------
# Radios de las capas (en unidades DPD) - valores ajustados
# ------------------------------------------------------------
# La bicapa tiene dos capas: externa (cabezas hacia fuera) e interna (cabezas hacia dentro).
# Las cabezas (L1) están en los radios más externo (10.4) e interno (7.6).
# Los cuellos (L2) en radios intermedios.
# Las colas (L3) se enfrentan en un radio común de 9.0, formando la zona hidrofóbica central.

# Capa externa (cabeza hacia afuera)
r_L1_ext = 10.4   # cabeza externa
r_L2_ext = 9.7    # cuello externo
r_L3_ext = 9.0    # cola externa

# Capa interna (cabeza hacia adentro)
r_L3_int = 9.0    # cola interna (mismo radio que la cola externa)
r_L2_int = 8.3    # cuello interno
r_L1_int = 7.6    # cabeza interna

# ------------------------------------------------------------
# Algoritmo de la espiral de Fibonacci para distribuir puntos uniformemente en una esfera
# ------------------------------------------------------------
angulo_dorado = math.pi * (3 - math.sqrt(5))   # ≈ 2.39996 rad

# ------------------------------------------------------------
# 2. APERTURA DEL ARCHIVO DE SALIDA Y ESCRITURA DE LA CABECERA
# ------------------------------------------------------------
with open("data.liposome", "w") as f:
    # Primera línea: título (obligatorio)
    f.write("LAMMPS data file for Liposome (lecithin + water)\n\n")
    
    # Número de átomos, enlaces y ángulos
    f.write(f"{total_atoms} atoms\n")
    f.write(f"{bonds} bonds\n")
    f.write(f"{angles} angles\n\n")
    
    # Número de tipos
    f.write(f"{atypes} atom types\n")
    f.write(f"{btypes} bond types\n")
    f.write(f"{antypes} angle types\n\n")
    
    # Caja de simulación: cúbica de lado 40, centrada en el origen
    # Esto da un volumen de 64000 u^3, y con densidad 3 tenemos 192000 partículas.
    f.write("-20.0 20.0 xlo xhi\n")
    f.write("-20.0 20.0 ylo yhi\n")
    f.write("-20.0 20.0 zlo zhi\n\n")
    
    # Masas: todos los beads tienen masa 1.0 (unidades DPD)
    f.write("Masses\n\n")
    f.write("1 1.0\n")   # L1
    f.write("2 1.0\n")   # L2
    f.write("3 1.0\n")   # L3
    f.write("4 1.0\n\n") # Agua
    
    # ------------------------------------------------------------
    # 3. SECCIÓN ATOMS
    # ------------------------------------------------------------
    f.write("Atoms\n\n")
    atom_id = 1        # identificador único para cada átomo
    mol_id = 1         # identificador de molécula (para la lecitina)
    
    # ------------------------------------------------------------
    # 3.1 CAPA EXTERNA (2465 moléculas)
    # ------------------------------------------------------------
    for i in range(mol_ext):
        # Coordenada y en el rango [1, -1] (Fibonacci sphere)
        y = 1 - (i / float(mol_ext - 1)) * 2
        theta = angulo_dorado * i
        radio_plano = math.sqrt(1 - y * y)
        
        # Vector dirección unitario (x, y, z) sobre la esfera de radio 1
        dir_x = math.cos(theta) * radio_plano
        dir_y = y
        dir_z = math.sin(theta) * radio_plano
        
        # Átomo L1 (cabeza, tipo 1) en el radio externo
        f.write(f"{atom_id} {mol_id} 1 0.0 {dir_x * r_L1_ext:.4f} {dir_y * r_L1_ext:.4f} {dir_z * r_L1_ext:.4f}\n")
        atom_id += 1
        
        # Átomo L2 (cuello, tipo 2) en radio intermedio externo
        f.write(f"{atom_id} {mol_id} 2 0.0 {dir_x * r_L2_ext:.4f} {dir_y * r_L2_ext:.4f} {dir_z * r_L2_ext:.4f}\n")
        atom_id += 1
        
        # Átomo L3 (cola, tipo 3) en radio común 9.0
        f.write(f"{atom_id} {mol_id} 3 0.0 {dir_x * r_L3_ext:.4f} {dir_y * r_L3_ext:.4f} {dir_z * r_L3_ext:.4f}\n")
        atom_id += 1
        
        mol_id += 1   # siguiente molécula
    
    # ------------------------------------------------------------
    # 3.2 CAPA INTERNA (2464 moléculas)
    # ------------------------------------------------------------
    for i in range(mol_int):
        y = 1 - (i / float(mol_int - 1)) * 2
        theta = angulo_dorado * i
        radio_plano = math.sqrt(1 - y * y)
        
        dir_x = math.cos(theta) * radio_plano
        dir_y = y
        dir_z = math.sin(theta) * radio_plano
        
        # Átomo L1 (cabeza, tipo 1) en el radio interno (ahora apunta hacia el interior)
        f.write(f"{atom_id} {mol_id} 1 0.0 {dir_x * r_L1_int:.4f} {dir_y * r_L1_int:.4f} {dir_z * r_L1_int:.4f}\n")
        atom_id += 1
        
        # Átomo L2 (cuello, tipo 2) en radio intermedio interno
        f.write(f"{atom_id} {mol_id} 2 0.0 {dir_x * r_L2_int:.4f} {dir_y * r_L2_int:.4f} {dir_z * r_L2_int:.4f}\n")
        atom_id += 1
        
        # Átomo L3 (cola, tipo 3) en el mismo radio común 9.0 (se enfrenta a la cola externa)
        f.write(f"{atom_id} {mol_id} 3 0.0 {dir_x * r_L3_int:.4f} {dir_y * r_L3_int:.4f} {dir_z * r_L3_int:.4f}\n")
        atom_id += 1
        
        mol_id += 1
    
    # ------------------------------------------------------------
    # 3.3 AGUA (partículas tipo 4) distribuidas aleatoriamente en toda la caja
    # ------------------------------------------------------------
    # El agua debe llenar tanto el exterior como el interior del liposoma.
    # Las posiciones aleatorias garantizan que, durante la simulación, el agua penetre
    # en la cavidad interna y rodee la bicapa.
    for _ in range(num_agua):
        x = random.uniform(-20.0, 20.0)
        y = random.uniform(-20.0, 20.0)
        z = random.uniform(-20.0, 20.0)
        # mol_id = 0 porque el agua son partículas independientes (no forman moléculas con enlaces)
        f.write(f"{atom_id} 0 4 0.0 {x:.4f} {y:.4f} {z:.4f}\n")
        atom_id += 1
    
    # ------------------------------------------------------------
    # 4. SECCIÓN BONDS (enlaces entre L1-L2 y L2-L3 para cada molécula de lecitina)
    # ------------------------------------------------------------
    f.write("\nBonds\n\n")
    bond_id = 1
    for m in range(total_moleculas):
        # Calculamos los ID de los tres átomos de la molécula m
        # Los átomos se numeraron consecutivamente: para la molécula 0, átomos 1,2,3; para la 1, 4,5,6; etc.
        id_L1 = (m * 3) + 1
        id_L2 = id_L1 + 1
        id_L3 = id_L2 + 1
        
        # Enlace L1-L2 (tipo 1)
        f.write(f"{bond_id} 1 {id_L1} {id_L2}\n")
        bond_id += 1
        # Enlace L2-L3 (tipo 1)
        f.write(f"{bond_id} 1 {id_L2} {id_L3}\n")
        bond_id += 1
    
    # ------------------------------------------------------------
    # 5. SECCIÓN ANGLES (ángulo L1-L2-L3 para cada molécula de lecitina)
    # ------------------------------------------------------------
    f.write("\nAngles\n\n")
    angle_id = 1
    for m in range(total_moleculas):
        id_L1 = (m * 3) + 1
        id_L2 = id_L1 + 1
        id_L3 = id_L2 + 1
        
        # Ángulo con vértice en L2 (tipo 1)
        f.write(f"{angle_id} 1 {id_L1} {id_L2} {id_L3}\n")
        angle_id += 1

# ------------------------------------------------------------
# Fin del script
# ------------------------------------------------------------
print("¡Archivo 'data.liposome' generado con éxito!")
print(f"Contiene {total_moleculas} moléculas de lecitina y {num_agua} partículas de agua.")
print(f"Total de átomos: {total_atoms} (densidad = {total_atoms / 64000:.2f})")
