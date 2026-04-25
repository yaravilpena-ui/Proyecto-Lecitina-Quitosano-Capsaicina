# -*- coding: utf-8 -*-
"""
Generador de configuración aleatoria para DPD (formato full).
Modelo: L1–L2–L3–L3–L3  (5 cuentas por lípido)
Tipos de cuenta: 1=W, 2=L1, 3=L2, 4=L3
Ángulos: UN SOLO TIPO (1) con k=50.0, θ0=170°
"""

import numpy as np
import math

# ========== PARÁMETROS ==========
N_lip = 4929
rho_global = 3.0
N_total = 150000
N_lip_beads = 5 * N_lip
N_water = N_total - N_lip_beads
L = (N_total / rho_global) ** (1.0/3.0)

print(f"Caja: {L:.3f}³, densidad: {rho_global}")
print(f"Cuentas: {N_total} total = {N_lip_beads} lípido + {N_water} agua")

r0 = 0.7
ks = 100.0
theta_deg = 170.0
bend_angle = math.radians(180.0 - theta_deg)   # 10°

atoms = []      # (id, mol, tipo, q, x, y, z)
bonds = []      # (tipo_enlace, id1, id2)
angles = []     # (tipo_angulo, id1, id2, id3)   <- ahora todos tipo 1

atom_id = 1
mol_id = 1

for i in range(N_lip):
    mol = mol_id + i

    # Posición aleatoria de L1
    l1 = np.random.uniform(-L/2, L/2, 3)

    # Dirección L1 -> L2
    dir_prev = np.random.randn(3)
    dir_prev /= np.linalg.norm(dir_prev)
    l2 = l1 + r0 * dir_prev

    # --- L3_1 ---
    refer = -dir_prev
    perp = np.random.randn(3)
    perp -= perp.dot(refer) * refer
    norm = np.linalg.norm(perp)
    if norm < 1e-12:
        perp = np.cross(refer, np.array([1,0,0]))
        norm = np.linalg.norm(perp)
    perp /= norm
    cos_b = math.cos(bend_angle)
    sin_b = math.sin(bend_angle)
    dir_new = cos_b * refer + sin_b * perp
    dir_new /= np.linalg.norm(dir_new)
    l3_1 = l2 + r0 * dir_new

    # --- L3_2 ---
    refer = -dir_new
    perp = np.random.randn(3)
    perp -= perp.dot(refer) * refer
    norm = np.linalg.norm(perp)
    if norm < 1e-12:
        perp = np.cross(refer, np.array([0,1,0]))
        norm = np.linalg.norm(perp)
    perp /= norm
    dir_next = cos_b * refer + sin_b * perp
    dir_next /= np.linalg.norm(dir_next)
    l3_2 = l3_1 + r0 * dir_next

    # --- L3_3 ---
    refer = -dir_next
    perp = np.random.randn(3)
    perp -= perp.dot(refer) * refer
    norm = np.linalg.norm(perp)
    if norm < 1e-12:
        perp = np.cross(refer, np.array([0,0,1]))
        norm = np.linalg.norm(perp)
    perp /= norm
    dir_last = cos_b * refer + sin_b * perp
    dir_last /= np.linalg.norm(dir_last)
    l3_3 = l3_2 + r0 * dir_last

    # IDs
    idL1 = atom_id; atom_id += 1
    idL2 = atom_id; atom_id += 1
    idL3a = atom_id; atom_id += 1
    idL3b = atom_id; atom_id += 1
    idL3c = atom_id; atom_id += 1

    q = 0.0
    atoms.append((idL1, mol, 2, q, l1[0], l1[1], l1[2]))
    atoms.append((idL2, mol, 3, q, l2[0], l2[1], l2[2]))
    atoms.append((idL3a, mol, 4, q, l3_1[0], l3_1[1], l3_1[2]))
    atoms.append((idL3b, mol, 4, q, l3_2[0], l3_2[1], l3_2[2]))
    atoms.append((idL3c, mol, 4, q, l3_3[0], l3_3[1], l3_3[2]))

    # Enlaces (tipo 1)
    bonds.append((1, idL1, idL2))
    bonds.append((1, idL2, idL3a))
    bonds.append((1, idL3a, idL3b))
    bonds.append((1, idL3b, idL3c))

    # Ángulos: TODOS tipo 1 (porque k y θ0 son idénticos)
    angles.append((1, idL1, idL2, idL3a))
    angles.append((1, idL2, idL3a, idL3b))
    angles.append((1, idL3a, idL3b, idL3c))

# Agua
for _ in range(N_water):
    pos = np.random.uniform(-L/2, L/2, 3)
    atoms.append((atom_id, mol_id, 1, 0.0, pos[0], pos[1], pos[2]))
    atom_id += 1
    mol_id += 1

# ----- Archivo data -----
with open("data.random_lecithin_5beads_170_full", "w") as f:
    f.write("Config aleatoria lecitina 5 beads – 1 tipo de ángulo\n\n")
    f.write(f"{atom_id-1} atoms\n")
    f.write(f"{len(bonds)} bonds\n")
    f.write(f"{len(angles)} angles\n")
    f.write("0 dihedrals\n")
    f.write("0 impropers\n\n")
    f.write("4 atom types\n")
    f.write("1 bond types\n")
    f.write("1 angle types\n\n")            # ← solo 1 tipo de ángulo
    f.write(f"{-L/2:.6f} {L/2:.6f} xlo xhi\n")
    f.write(f"{-L/2:.6f} {L/2:.6f} ylo yhi\n")
    f.write(f"{-L/2:.6f} {L/2:.6f} zlo zhi\n\n")
    f.write("Masses\n\n")
    f.write("1 1.0\n2 1.0\n3 1.0\n4 1.0\n\n")

    f.write("Atoms # full\n\n")
    for a in atoms:
        f.write(f"{a[0]} {a[1]} {a[2]} {a[3]:.6f} {a[4]:.6f} {a[5]:.6f} {a[6]:.6f}\n")

    f.write("\nBonds\n\n")
    for idx, b in enumerate(bonds, start=1):
        f.write(f"{idx} {b[0]} {b[1]} {b[2]}\n")

    f.write("\nAngles\n\n")
    for idx, a in enumerate(angles, start=1):
        f.write(f"{idx} {a[0]} {a[1]} {a[2]} {a[3]}\n")

print("Archivo 'data.random_lecithin_5beads_170_full' creado.")
