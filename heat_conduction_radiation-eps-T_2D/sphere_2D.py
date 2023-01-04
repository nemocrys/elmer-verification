import os
import gmsh
from pyelmer import elmer


def mesh(heater_r_in, heater_r_out, insulation_r_in, insulation_r_out, mesh_size):
    # gmsh initialization
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("2D-sphere")
    factory = gmsh.model.occ

    # modeling
    heater_1d = factory.addCircle(0, 0, 0, heater_r_out)
    heater = factory.addSurfaceFilling(factory.addCurveLoop([heater_1d]))
    heater_hole_1d = factory.addCircle(0, 0, 0, heater_r_in)
    heater_hole = factory.addSurfaceFilling(factory.addCurveLoop([heater_hole_1d]))
    factory.cut([(2, heater)], [(2, heater_hole)])

    insulation_1d = factory.addCircle(0, 0, 0, insulation_r_out)
    insulation = factory.addSurfaceFilling(factory.addCurveLoop([insulation_1d]))
    insulation_hole_1 = factory.addCircle(0, 0, 0, insulation_r_in)
    insulation_hole = factory.addSurfaceFilling(
        factory.addCurveLoop([insulation_hole_1])
    )
    factory.cut([(2, insulation)], [(2, insulation_hole)])

    box = factory.addRectangle(0, -1, 0, -1, 2)
    factory.cut([(2, heater), (2, insulation)], [(2, box)])

    factory.synchronize()

    ph_heater = gmsh.model.addPhysicalGroup(2, [heater])
    gmsh.model.setPhysicalName(2, ph_heater, "heater")
    ph_insulation = gmsh.model.addPhysicalGroup(2, [insulation])
    gmsh.model.setPhysicalName(2, ph_insulation, "insulation")

    # arc tags from visual output
    heater_inside = [3, 4]
    ph_heater_in = gmsh.model.addPhysicalGroup(1, heater_inside)
    gmsh.model.setPhysicalName(1, ph_heater_in, "heater_inside")
    heater_outside = [1, 6]
    ph_heater_out = gmsh.model.addPhysicalGroup(1, heater_outside)
    gmsh.model.setPhysicalName(1, ph_heater_out, "heater_outside")
    insulation_inside = [9, 10]
    ph_insulation_in = gmsh.model.addPhysicalGroup(1, insulation_inside)
    gmsh.model.setPhysicalName(1, ph_insulation_in, "insulation_inside")
    insulation_outside = [7, 12]
    ph_insulation_out = gmsh.model.addPhysicalGroup(1, insulation_outside)
    gmsh.model.setPhysicalName(1, ph_insulation_out, "insulation_outside")

    # visualize & export mesh
    factory.synchronize()
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    gmsh.model.mesh.generate(2)
    # gmsh.fltk.run()
    if not os.path.exists("./simdata"):
        os.mkdir("./simdata")
    gmsh.write("./simdata/2d_sphere.msh")
    gmsh.finalize()

    return (
        ph_heater,
        ph_insulation,
        ph_heater_in,
        ph_heater_out,
        ph_insulation_in,
        ph_insulation_out,
    )


def elmer_setup(
    ph_heater,
    ph_insulation,
    ph_heater_in,
    ph_heater_out,
    ph_insulation_in,
    ph_insulation_out,
):
    sim = elmer.Simulation()
    sim.settings = {
        "Max Output Level": 4,
        "Coordinate System": "Axi Symmetric",
        "Simulation Type": "Steady state",
        "Steady State Max Iterations": 10,
        "Output File": "case.result",
    }

    solver_heat = elmer.load_solver("HeatSolver", sim, "./solvers.yml")
    solver_output = elmer.load_solver("ResultOutputSolver", sim, "./solvers.yml")
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

    bdy_heater = elmer.Body(sim, "heater", [ph_heater])
    bdy_heater.material = mat_heater
    bdy_heater.equation = equation
    bdy_heater.initial_condition = initial_temp
    bdy_heater.body_force = heating

    bdy_insulation = elmer.Body(sim, "insulation", [ph_insulation])
    bdy_insulation.material = mat_insulation
    bdy_insulation.equation = equation
    bdy_insulation.initial_condition = initial_temp

    bnd_heater_in = elmer.Boundary(sim, "heater_in", [ph_heater_in])
    bnd_heater_in.data = {"Radiation": "Diffuse Gray"}
    bnd_heater_out = elmer.Boundary(sim, "heater_out", [ph_heater_out])
    bnd_heater_out.data = {"Radiation": "Diffuse Gray"}
    bnd_insulation_in = elmer.Boundary(sim, "insulation_in", [ph_insulation_in])
    bnd_insulation_in.data = {"Radiation": "Diffuse Gray"}
    bnd_insulation_out = elmer.Boundary(sim, "insulation_out", [ph_insulation_out])
    bnd_insulation_out.data = {
        "Radiation": "Idealized",
        "External Temperature": 300,
    }
    sim.write_startinfo("./simdata")
    sim.write_sif("./simdata")


if __name__ == "__main__":
    heater_r_in = 0.4
    heater_r_out = 0.5
    insulation_r_in = 0.9
    insulation_r_out = 1
    mesh_size = 0.01

    (
        ph_heater,
        ph_insulation,
        ph_heater_in,
        ph_heater_out,
        ph_insulation_in,
        ph_insulation_out,
    ) = mesh(heater_r_in, heater_r_out, insulation_r_in, insulation_r_out, mesh_size)
    elmer_setup(
        ph_heater,
        ph_insulation,
        ph_heater_in,
        ph_heater_out,
        ph_insulation_in,
        ph_insulation_out,
    )

    from pyelmer.execute import run_elmer_grid, run_elmer_solver

    run_elmer_grid("./simdata", "2d_sphere.msh")
    run_elmer_solver("./simdata")
    from pyelmer.post import scan_logfile

    warn, err, stats = scan_logfile("./simdata")
    print("Warnings:", warn)
    print("Errors:", err)
    print("Statistics:", stats)
