import os
import gmsh
import numpy as np
from objectgmsh import Model, Shape, MeshControlExponential, cut
from pyelmer import elmerkw as elmer, execute, post


if not os.path.exists("Lorentz_force/simdata"):
    os.makedirs("Lorentz_force/simdata")

occ = gmsh.model.occ

####################
# mesh setup
####################

model = Model()
cylinder = Shape(model, 2, "cylinder", [occ.add_rectangle(0, -0.1, 0, 0.05, 0.2)])
coil_0deg = Shape(model, 2, "coil_0deg", [occ.add_rectangle(0.2, -0.3, 0, 0.01, 0.01)])
coil_60deg = Shape(
    model, 2, "coil_60deg", [occ.add_rectangle(0.2, -0.182, 0, 0.01, 0.01)]
)
coil_120deg = Shape(
    model, 2, "coil_120deg", [occ.add_rectangle(0.2, -0.064, 0, 0.01, 0.01)]
)
coil_180deg = Shape(
    model, 2, "coil_180deg", [occ.add_rectangle(0.2, 0.054, 0, 0.01, 0.01)]
)
coil_240deg = Shape(
    model, 2, "coil_240deg", [occ.add_rectangle(0.2, 0.172, 0, 0.01, 0.01)]
)
coil_300deg = Shape(
    model, 2, "coil_300deg", [occ.add_rectangle(0.2, 0.29, 0, 0.01, 0.01)]
)

surrounding = Shape(model, 2, "surrounding", [occ.add_rectangle(0, -0.7, 0, 0.7, 1.4)])
surrounding.geo_ids = cut(
    surrounding.dimtags,
    cylinder.dimtags
    + coil_0deg.dimtags
    + coil_60deg.dimtags
    + coil_120deg.dimtags
    + coil_180deg.dimtags
    + coil_240deg.dimtags
    + coil_300deg.dimtags,
    remove_tool=False,
)
model.synchronize()
outside = Shape(
    model,
    1,
    "outside",
    [surrounding.bottom_boundary, surrounding.right_boundary, surrounding.top_boundary],
)
cylinder_surface = Shape(model, 1, "cylinder_surface", [cylinder.right_boundary])

model.synchronize()
model.make_physical()

model.deactivate_characteristic_length()
for shape in [
    coil_0deg,
    coil_60deg,
    coil_120deg,
    coil_180deg,
    coil_240deg,
    coil_300deg,
]:
    MeshControlExponential(model, shape, 0.0025, exp=2.2)

model.generate_mesh()
model.write_msh("Lorentz_force/simdata/case.msh")
# model.show()

####################
# simulation setup
####################
omega = 2 * np.pi * 500  # 1/s
current_dty = 1000 / 0.01**2  # A/m^2

sim = elmer.load_simulation("axi-symmetric_steady", "Lorentz_force/config_elmer.yml")
solver_mgdyn = elmer.load_solver(
    "MagnetoDynamics2DHarmonic", sim, "Lorentz_force/config_elmer.yml"
)
solver_mgdyn.data.update({"Angular frequency": omega})
solver_mgdyncalc = elmer.load_solver(
    "MagnetoDynamicsCalcFields", sim, "Lorentz_force/config_elmer.yml"
)
solver_mgdyncalc.data.update({"Angular frequency": omega})
solver_output = elmer.load_solver(
    "ResultOutputSolver", sim, "Lorentz_force/config_elmer.yml"
)

equation = elmer.Equation(sim, "eqn", [solver_mgdyn, solver_mgdyncalc])

insulator = elmer.Material(
    sim,
    "insulator",
    {
        "Electric conductivity": 0,
        "Relative Permeability": 1,
    },
)
conductor = elmer.Material(
    sim,
    "conductor",
    {
        "Electric conductivity": 1e6,
        "Relative Permeability": 1,
    },
)

for coil, phase in zip(
    [coil_0deg, coil_60deg, coil_120deg, coil_180deg, coil_240deg, coil_300deg],
    [0, 60, 120, 180, 240, 300],
):
    body = elmer.Body(sim, coil.name, [coil.ph_id])
    force = elmer.BodyForce(sim, coil.name + "_current")
    current_dty_re = np.cos(2 * np.pi * phase / 360) * current_dty
    current_dty_im = np.sin(2 * np.pi * phase / 360) * current_dty
    force.data.update(
        {
            "Current density": f"real {current_dty_re}",
            "Current density Im": f"real {current_dty_im}",
        }
    )
    body.body_force = force
    body.material = insulator
    body.equation = equation

cylinder = elmer.Body(sim, "cylinder", [cylinder.ph_id])
cylinder.material = conductor
cylinder.equation = equation

surrounding = elmer.Body(sim, "surrounding", [surrounding.ph_id])
surrounding.material = insulator
surrounding.equation = equation

outside = elmer.Boundary(sim, "outside", [outside.ph_id])
outside.zero_potential = True
sim.write_sif("Lorentz_force/simdata")


####################
# execute simulation
####################
execute.run_elmer_grid("Lorentz_force/simdata", "case.msh")
execute.run_elmer_solver("Lorentz_force/simdata")
err, warn, stats = post.scan_logfile("Lorentz_force/simdata")
print("Errors:", err)
print("Warnings:", warn)
print("Statistics:", stats)
