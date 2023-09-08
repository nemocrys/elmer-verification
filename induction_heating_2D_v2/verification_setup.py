import os
import gmsh
import numpy as np
import yaml
from objectgmsh import Model, Shape, MeshControlExponential
from pyelmer import elmerkw as elmer
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile


config_file = "induction_heating_2D_v2/config.yml"
use_statmag = True

with open(config_file) as f:
    config = yaml.safe_load(f)

wdir = f"induction_heating_2D_v2/simdata_air-size={config['air_size_factor']*(config['coil_d']+config['coil_r_i'])}_mesh-size={config['mesh_size_factor']}"
if not os.path.exists(wdir):
    os.makedirs(wdir)

########################################
# geometry modeling
########################################
model = Model()
occ = gmsh.model.occ

# parameters
air_size = config["air_size_factor"] * (config["coil_r_i"] + config["coil_d"])

# geometry modeling
cylinder = Shape(model, 2, "cylinder", [occ.addRectangle(0, -config["height"]/2, 0, config["graphite_r"], config["height"])])
coil = Shape(model, 2, "coil", [occ.addRectangle(config["coil_r_i"], -config["height"]/2, 0, config["coil_d"], config["height"])])
air = Shape(model, 2, "air", [occ.addRectangle(0, -config["height"]/2, 0, air_size, config["height"])])
air.geo_ids = [x[1] for x in occ.cut(air.dimtags, cylinder.dimtags + coil.dimtags, removeTool=False)[0]]
occ.synchronize()

# boundaries
cylinder_surf = Shape(model, 1, "cylinder_surf", [cylinder.right_boundary])
outside_surf = Shape(model, 1, "outside_surf", [air.right_boundary])

model.make_physical()

# mesh
mesh_size_base = config["height"] * 0.2
mesh_size_min = config["height"] * 0.01
cylinder.mesh_size = mesh_size_base
air.mesh_size = mesh_size_base
coil.mesh_size = mesh_size_base
MeshControlExponential(model, coil, mesh_size_min)

model.deactivate_characteristic_length()
model.set_const_mesh_sizes()
model.generate_mesh(size_factor=config["mesh_size_factor"])
model.write_msh(f"{wdir}/case.msh")

########################################
# simulation setup
########################################

sim = elmer.load_simulation("axi-symmetric_steady", config_file)

omega = 2*np.pi*config["frequency"]
if use_statmag:
    solver_statmag = elmer.load_solver("StatMagSolver", sim, config_file)
    solver_statmag.data.update({"Angular Frequency": omega})
else:
    solver_mgdyn = elmer.load_solver("StatMMagnetoDynamics2DHarmonicagSolver", sim, config_file)
    solver_mgdyn.data.update({"Angular Frequency": omega})
    solver_mgdyncalc = elmer.load_solver("MagnetoDynamicsCalcFields", sim, config_file)
    solver_mgdyncalc.data.update({"Angular Frequency": omega})
solver_heat = elmer.load_solver("HeatSolver", sim, config_file)
solver_res = elmer.load_solver("ResultOutputSolver", sim, config_file)
solver_save_line = elmer.load_solver("SaveLine", sim, config_file)

eqn_main = elmer.Equation(sim, "eqn_main", [])
if use_statmag:
    eqn_main.solvers = [solver_statmag, solver_heat]
else:
    eqn_main.solvers = [solver_mgdyn, solver_mgdyncalc, solver_heat]

graphite = elmer.load_material("graphite", sim, config_file)
insulator = elmer.load_material("insulator", sim, config_file)

current_dty = elmer.BodyForce(sim, "current_dty", {"Current Density": config["current"] / (config["height"] * config["coil_d"])})
joule_heat = elmer.BodyForce(sim, "joule_heat")
joule_heat.joule_heat = True

cylinder = elmer.Body(sim, "cylinder", [cylinder.ph_id])
cylinder.material = graphite
cylinder.equation = eqn_main
cylinder.body_force = joule_heat

air = elmer.Body(sim, "air", [air.ph_id])
air.material = air
air.equation = eqn_main

coil = elmer.Body(sim, "coil", [coil.ph_id])
coil.material = insulator
coil.equation = eqn_main
coil.body_force = current_dty

cylinder_surf = elmer.Boundary(sim, "cylinder_surf", [cylinder_surf.ph_id])
cylinder_surf.radiation_idealized = True
cylinder_surf.T_ext = config["T_amb"]
cylinder_surf.heat_transfer_coefficient = config["htc_graphite"]

outside_surf = elmer.Boundary(sim, "outside_surf", [outside_surf.ph_id])
outside_surf.zero_potential = True

sim.write_sif(wdir)

run_elmer_grid(wdir, "case.msh")
run_elmer_solver(wdir)

err, warn, stats = scan_logfile(wdir)
print("Errors:", err)
print("Warnings", warn)
print("Statistics:", stats)
