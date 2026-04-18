import hoomd
import gsd.hoomd
import numpy as np

frame = gsd.hoomd.Frame()

# Place a polymer in the box.
frame.particles.N = 200
frame.particles.types = ["W"]
frame.particles.typeid = [0] * 200
frame.configuration.box = [20, 20, 20, 0, 0, 0]
posiciones= np.random.uniform(low=-10,high=10,size=(200,3))
frame.particles.position = posiciones

# Connect particles with bonds
#frame.bonds.N = 49
#frame.bonds.types = ["Ch-Ch"]
#frame.bonds.typeid = [0] * 49
#frame.bonds.group = [[i,i+ 1] for i in range (49)]

with gsd.hoomd.open(name="H2O.gsd", mode="x") as f:
    f.append(frame)

# Apply the harmonic potential on the bonds.
#harmonic = hoomd.md.bond.Harmonic()
#harmonic.params["Ch-Ch"] = dict(k=100, r0=0.7)

# Perform the MD simulation.
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_gsd(filename="H2O.gsd")
langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
integrator = hoomd.md.Integrator(dt=0.005, methods=[langevin])
gsd_writer = hoomd.write.GSD(
    filename="H2O_trajectory.gsd", trigger=hoomd.trigger.Periodic(1000), mode="xb"
)
sim.operations.integrator = integrator
sim.operations.writers.append(gsd_writer)
sim.run(10e3)


