#Importar librerias
import hoomd
import gsd.hoomd
import numpy as np

#Creara las posiciones random del agua
n_agua=191950 #el numero que se calsulo de acuerdo a la densidad
l_caja=40
posiciones_agua = np.random.uniform(low=-l_caja/2, high=l_caja/2, size=(n_agua, 3)) #crea la matriz con 3 coordenadas random para cada molécula de agua

frame = gsd.hoomd.Frame()

# Armar los beads y las posiciones del polímero
frame.particles.N = 192000  #el numero de particulas totales del sistema para culplir con la densidad
frame.particles.types = ["G","A","W"] #Se definen los nombres de los beads que habrá en el sistema
bloques_base=[1,0,0,0,1,0,0,0,1,0] #en la ista anterior se guarda como G=0, A=1, W=2. Entonces qui le decimos que el polimero tiene un orden de 1,0... Este es el bloque base, como el monómero
lista_poli=bloques_base*5 #Se crea la lista del polímero, que son 5 veces lo de arriba
lista_agua=[2]*n_agua #creamos la lista para el agua
frame.particles.typeid = lista_poli+lista_agua #Se asigna un nombre a cada bead segun su posición en la lista
frame.configuration.box = [l_caja, l_caja, l_caja, 0, 0, 0] #Define la caja
posiciones_poli = [[i*0.7-17.15,0,0] for i in range (50)] #Le estamos diciendo como acomodar los beads del polímero
posiciones = np.concatenate((posiciones_poli,posiciones_agua)) #unimos las listas de las posiciones
frame.particles.position = posiciones
# Connect particles with bonds
frame.bonds.N = 49 #porque solo queremos conectar los del quitosano
frame.bonds.types = ["Ch-Ch"]
frame.bonds.typeid = [0] * 49
frame.bonds.group = [[i,i+ 1] for i in range (49)]

with gsd.hoomd.open(name="Ch_H2O.gsd", mode="x") as f:
    f.append(frame)

# Apply the harmonic potential on the bonds.
harmonic = hoomd.md.bond.Harmonic()
harmonic.params["Ch-Ch"] = dict(k=100, r0=0.7)

#Interaction parameters
dpd = hoomd.md.pair.DPD(nlist=hoomd.md.nlist.Cell(buffer=0.4), kT=1.0, default_r_cut=1.0)

# Interacciones con el agua (W)
dpd.params[('G', 'W')] = dict(A=80.98, gamma=4.5)
dpd.params[('A', 'W')] = dict(A=83.41, gamma=4.5)

# Interacciones iguales
dpd.params[('W', 'W')] = dict(A=78.33, gamma=4.5)
dpd.params[('G', 'G')] = dict(A=78.33, gamma=4.5)
dpd.params[('A', 'A')] = dict(A=78.33, gamma=4.5)

# Interacciones A-G
dpd.params[('A', 'G')] = dict(A=78.68, gamma=4.5)

# Perform the MD simulation.
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_gsd(filename="Ch_H2O.gsd")
#langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
# En lugar de Langevin, usamos esto para DPD
integrator = hoomd.md.Integrator(dt=0.03, forces=[dpd, harmonic])
integrator.methods.append(hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All()))
gsd_writer = hoomd.write.GSD(
    filename="Ch_H2Otrajectory.gsd", trigger=hoomd.trigger.Periodic(1000), mode="xb"
)
sim.operations.integrator = integrator
sim.operations.writers.append(gsd_writer)
sim.run(5e6)

