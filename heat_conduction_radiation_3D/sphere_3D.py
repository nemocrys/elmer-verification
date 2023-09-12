import os
import gmsh
from pyelmer import elmerkw as elmer
from pyelmer.execute import run_elmer_grid, run_elmer_solver
from pyelmer.post import scan_logfile
from objectgmsh import Model, Shape, rotate
import gmsh
occ = gmsh.model.occ


##########
heater_r_in = 0.4
heater_r_out = 0.5
insulation_r_in = 0.9
insulation_r_out = 1
mesh_size = 0.01
size_factor = 1
##########

simdir = f"./simdata_sf-{size_factor}"

model = Model()

# modeling 2D
heater_1d = occ.addCircle(0, 0, 0, heater_r_out)
heater = occ.addSurfaceFilling(occ.addCurveLoop([heater_1d]))
heater_hole_1d = occ.addCircle(0, 0, 0, heater_r_in)
heater_hole = occ.addSurfaceFilling(occ.addCurveLoop([heater_hole_1d]))
occ.cut([(2, heater)], [(2, heater_hole)])

insulation_1d = occ.addCircle(0, 0, 0, insulation_r_out)
insulation = occ.addSurfaceFilling(occ.addCurveLoop([insulation_1d]))
insulation_hole_1 = occ.addCircle(0, 0, 0, insulation_r_in)
insulation_hole = occ.addSurfaceFilling(
    occ.addCurveLoop([insulation_hole_1])
)
occ.cut([(2, insulation)], [(2, insulation_hole)])
box = occ.addRectangle(0, -1, 0, -1, 2)
occ.cut([(2, heater), (2, insulation)], [(2, box)])

# modeling 3D
heater = Shape(model, 3, "heater", [rotate(heater)])
insulation = Shape(model, 3, "insulation", [rotate(insulation)])

model.synchronize()

heater_outside = Shape(model, 2, "heater_outside", heater.boundaries[:2])
heater_inside = Shape(model, 2, "heater_inside", heater.boundaries[2:])

insulation_outside = Shape(model, 2, "insulation_outside", insulation.boundaries[:2])
insulation_inside = Shape(model, 2, "insulation_inside", insulation.boundaries[2:])

model.make_physical()

# visualize & export mesh
heater.mesh_size = mesh_size
insulation.mesh_size = mesh_size
model.deactivate_characteristic_length()
model.set_const_mesh_sizes()
model.generate_mesh(3, size_factor=size_factor, optimize="Netgen")
if not os.path.exists(simdir):
    os.mkdir(simdir)
model.write_msh(f"{simdir}/case.msh")



sim = elmer.Simulation()
sim.settings = {
    "Max Output Level": 4,
    "Coordinate System": "Cartesian",
    "Simulation Type": "Steady state",
    "Steady State Max Iterations": 10,
    "Output File": "case.result",
}

solver_heat = elmer.load_solver("HeatSolver", sim, "./solvers.yml")
solver_output = elmer.load_solver("ResultOutputSolver", sim, "./solvers.yml")
solver_save_line = elmer.load_solver("SaveLine", sim, "./solvers.yml")
equation = elmer.Equation(sim, "main", [solver_heat])

initial_temp = elmer.InitialCondition(sim, "T_init", {"Temperature": 300})
heating = elmer.BodyForce(sim, "heating")
heating.data.update({"Heat Source": 30000, "Integral Heat Source": 30000})

mat_heater = elmer.Material(sim, "mat_heater")
mat_heater.data = {
    "Heat Conductivity": 20.0,
    "Density": 1.0,
    "Heat Capacity": 1.0,
    "Emissivity": 0.8,
}
mat_insulation = elmer.Material(sim, "mat_insulation")
mat_insulation.data = {
    "Heat Conductivity": 0.5,
    "Density": 1.0,
    "Heat Capacity": 1.0,
    "Emissivity": 0.5,
}

bdy_heater = elmer.Body(sim, "heater", [heater.ph_id])
bdy_heater.material = mat_heater
bdy_heater.equation = equation
bdy_heater.initial_condition = initial_temp
bdy_heater.body_force = heating

bdy_insulation = elmer.Body(sim, "insulation", [insulation.ph_id])
bdy_insulation.material = mat_insulation
bdy_insulation.equation = equation
bdy_insulation.initial_condition = initial_temp

bnd_heater_in = elmer.Boundary(sim, "heater_in", [heater_inside.ph_id])
bnd_heater_in.radiation = True
bnd_heater_out = elmer.Boundary(sim, "heater_out", [heater_outside.ph_id])
bnd_heater_out.radiation = True
bnd_insulation_in = elmer.Boundary(sim, "insulation_in", [insulation_inside.ph_id])
bnd_insulation_in.radiation = True
bnd_insulation_out = elmer.Boundary(sim, "insulation_out", [insulation_outside.ph_id])
bnd_insulation_out.radiation_idealized = True
bnd_insulation_out.T_ext = 300

sim.write_sif(simdir)

run_elmer_grid(simdir, "case.msh")
run_elmer_solver(simdir)

warn, err, stats = scan_logfile(simdir)
print("Warnings:", warn)
print("Errors:", err)
print("Statistics:", stats)
