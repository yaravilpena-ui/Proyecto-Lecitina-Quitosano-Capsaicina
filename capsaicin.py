import hoomd
import gsd.hoomd

frame = gsd.hoomd.Frame()

# Place a polymer in the box.
frame.particles.N = 4
frame.particles.types = ["C1","C2","C3"]
bloques_base=[0,1,2,2]
frame.particles.typeid = bloques_base
frame.configuration.box = [70, 70, 70, 0, 0, 0]
frame.particles.position = [[i*0.7,0,0] for i in range (4)]

# Connect particles with bonds
frame.bonds.N = 3
frame.bonds.types = ["Cp-Cp"]
frame.bonds.typeid = [0] * 3
frame.bonds.group = [[i,i+ 1] for i in range (3)]

with gsd.hoomd.open(name="Cp.gsd", mode="x") as f:
    f.append(frame)

# Apply the harmonic potential on the bonds.
harmonic = hoomd.md.bond.Harmonic()
harmonic.params["Cp-Cp"] = dict(k=100, r0=0.7)

# Perform the MD simulation.
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_gsd(filename="Cp.gsd")
langevin = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=1.0)
integrator = hoomd.md.Integrator(dt=0.005, methods=[langevin], forces=[harmonic])
gsd_writer = hoomd.write.GSD(
    filename="Cp_trajectory.gsd", trigger=hoomd.trigger.Periodic(1000), mode="xb"
)
sim.operations.integrator = integrator
sim.operations.writers.append(gsd_writer)
sim.run(10e3)
                     
