import os
import shutil
import subprocess
import gmsh
import pyelmer.elmerkw as elmer
from pyelmer.gmsh import *


def geometry(l, n_l, r_e, r_i, d, alpha_mesh=1, with_cylinder=True):

    # initialize
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add('Induction Verification 2D')
    factory = gmsh.model.occ

    # geometry modeling
    l_ends = (n_l - 1) / 2 * l
    if with_cylinder:
        cylinder_bottom_end = factory.addRectangle(0, -l_ends, 0, r_e, l_ends)
        cylinder = factory.addRectangle(0, 0, 0, r_e, l)
        factory.fragment([(2, cylinder_bottom_end)], [(2, cylinder)])
        cylinder_top_end = factory.addRectangle(0, l, 0, r_e, l_ends)
        factory.fragment([(2, cylinder)], [(2, cylinder_top_end)])
    coil = factory.addRectangle(r_i, -l_ends, 0 , d, n_l*l)
    air = factory.addRectangle(0, -5*l - l_ends , 0, 10*l, 2*l_ends + 11* l)
    if with_cylinder:
        bodies = [(2, cylinder_bottom_end), (2, cylinder), (2, cylinder_top_end), (2, coil)]
    else:
        bodies = [(2, coil)]
    factory.cut([(2, air)], bodies, removeTool=False)
    factory.synchronize()

    # gmsh.fltk.run()

    # bodies
    if with_cylinder:
        ph_cylinder = add_physical_group(2, [cylinder], 'cylinder')
        ph_cylinder_ends = add_physical_group(2, [cylinder_bottom_end, cylinder_top_end], 'cylinder_ends')
    ph_coil = add_physical_group(2, [coil], 'coil')
    ph_air = add_physical_group(2, [air], 'air')

    # boundaries
    if with_cylinder:
        surfs = []
        surfs.append(get_cylinder_boundary(2, cylinder, r_e, l, 0))  # cylinder side
        ph_cylinder_surf = add_physical_group(1, surfs, 'cylinder_surf')
    surfs = []
    surfs.append(get_ring_boundary(2, air, 10*l,-5*l - l_ends))  # air bottom
    surfs.append(get_ring_boundary(2, air, 10*l, 6*l + l_ends))  # air top
    surfs.append(get_cylinder_boundary(2, air, 10*l, 2*l_ends + 11* l, -5*l - l_ends))
    ph_outside_surfs = add_physical_group(1, surfs, 'outside_surfs')


    # mesh
    lc_base = l * 0.5
    lc_min = l * 0.01
    mesh = gmsh.model.mesh
    gmsh.option.setNumber('Mesh.CharacteristicLengthFromPoints', 0)
    gmsh.option.setNumber('Mesh.CharacteristicLengthFromCurvature', 0)
    gmsh.option.setNumber('Mesh.CharacteristicLengthExtendFromBoundary', 0)
    mesh.setSize(gmsh.model.getEntities(0), lc_base)

    # generate fields for mesh size control
    fields = []
    fields.append(restricted_const_field(air, lc_base))
    if with_cylinder:
        fields.append(restricted_const_field(cylinder, lc_base))
        fields.append(restricted_const_field(cylinder_bottom_end, lc_base))
        fields.append(restricted_const_field(cylinder_top_end, lc_base))
    fields.append(exp_field(coil, lc_min, exp=1.8, fact=1))

    min_field = mesh.field.add('Min')
    mesh.field.setNumbers(min_field, 'FieldsList', fields)
    mesh.field.setAsBackgroundMesh(min_field)
    gmsh.option.setNumber('Mesh.CharacteristicLengthFactor', alpha_mesh)
    mesh.generate(2)


    # visualize & export
    gmsh.fltk.run()
    gmsh.write('./simdata/induction_verification.msh2')
    gmsh.finalize()

    if with_cylinder:
        return ph_cylinder, ph_cylinder_ends, ph_coil, ph_air, ph_cylinder_surf, ph_outside_surfs
    else:
        return 0, 0, ph_coil, ph_air, 0, ph_outside_surfs

def geometry_quad(l, n_l, r_e, r_i, d, alpha_mesh=1):

    # initialize
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add('Induction Verification 2D')
    factory = gmsh.model.occ

    # geometry modeling
    l_ends = (n_l - 1) / 2 * l
    cylinder_bottom_end = factory.addRectangle(0, -l_ends, 0, r_e, l_ends)
    cylinder = factory.addRectangle(0, 0, 0, r_e, l)
    factory.fragment([(2, cylinder_bottom_end)], [(2, cylinder)])
    cylinder_top_end = factory.addRectangle(0, l, 0, r_e, l_ends)
    factory.fragment([(2, cylinder)], [(2, cylinder_top_end)])
    coil = factory.addRectangle(r_i, -l_ends, 0 , d, n_l*l)

    air1 = factory.addRectangle(r_e, -l_ends, 0, r_i-r_e, n_l*l)
    factory.fragment([(2, air1)],[(2, cylinder), (2, cylinder_bottom_end), (2, cylinder_top_end), (2, coil)])
    air2 = factory.addRectangle(0, l_ends + l, 0, 10*l, 5*l)
    factory.fragment([(2, air2)],[(2, cylinder_top_end), (2, coil), (2, air1)])
    air3 = factory.addRectangle(0, -l_ends - 5*l, 0, 10*l, 5*l)
    factory.fragment([(2, air3)],[(2, cylinder_bottom_end), (2, coil), (2, air1)])
    air4 = factory.addRectangle(r_i + d, -l_ends, 0, 10*l - r_i - d, n_l*l)
    factory.fragment([(2, air4)],[(2, coil), (2, air2), (2, air3)])

    factory.synchronize()

    # bodies
    ph_cylinder = add_physical_group(2, [cylinder], 'cylinder')
    ph_cylinder_ends = add_physical_group(2, [cylinder_bottom_end, cylinder_top_end], 'cylinder_ends')
    ph_coil = add_physical_group(2, [coil], 'coil')
    ph_air = add_physical_group(2, [air1, air2, air3, air4], 'air')

    # boundaries
    surfs = []
    surfs.append(get_cylinder_boundary(2, cylinder, r_e, l, 0))  # cylinder side
    ph_cylinder_surf = add_physical_group(1, surfs, 'cylinder_surf')
    surfs = []
    surfs.append(get_ring_boundary(2, air2, 10*l, 6*l + l_ends))  # air top
    surfs.append(get_ring_boundary(2, air3, 10*l,-5*l - l_ends))  # air bottom
    surfs.append(get_cylinder_boundary(2, air2, 10*l, 5*l, l + l_ends))
    surfs.append(get_cylinder_boundary(2, air3, 10*l, 5*l, -l_ends - 5*l))
    surfs.append(get_cylinder_boundary(2, air4, 10*l, 2*l_ends + l, - l_ends))
    ph_outside_surfs = add_physical_group(1, surfs, 'outside_surfs')

    # mesh
    mesh = gmsh.model.mesh
    lc_base = l * 0.5
    lc_min = l * 0.01
    # gmsh.option.setNumber('Mesh.CharacteristicLengthFromPoints', 0)
    # gmsh.option.setNumber('Mesh.CharacteristicLengthFromCurvature', 0)
    # gmsh.option.setNumber('Mesh.CharacteristicLengthExtendFromBoundary', 0)
    mesh.setSize(gmsh.model.getEntities(0), lc_base)

    # generate fields for mesh size control
    fields = []
    fields.append(restricted_const_field(air1, lc_base))
    fields.append(restricted_const_field(air2, lc_base))
    fields.append(restricted_const_field(air3, lc_base))
    fields.append(restricted_const_field(air4, lc_base))
    fields.append(restricted_const_field(cylinder, lc_base))
    fields.append(restricted_const_field(cylinder_bottom_end, lc_base))
    fields.append(restricted_const_field(cylinder_top_end, lc_base))

    fields.append(threshold_field(coil, lc_min, lc_base, min_dist=lc_min, max_dist=lc_base*3))
    fields.append(threshold_field(air1, lc_min, lc_base, min_dist=lc_min, max_dist=lc_base*3))
    fields.append(threshold_field(cylinder, lc_min, lc_base, min_dist=lc_min, max_dist=lc_base*3))
    fields.append(threshold_field(cylinder_bottom_end, lc_min, lc_base, min_dist=lc_min, max_dist=lc_base*3))
    fields.append(threshold_field(cylinder_top_end, lc_min, lc_base, min_dist=lc_min, max_dist=lc_base*3))

    min_field = mesh.field.add('Min')
    mesh.field.setNumbers(min_field, 'FieldsList', fields)
    mesh.field.setAsBackgroundMesh(min_field)
    gmsh.option.setNumber('Mesh.CharacteristicLengthFactor', alpha_mesh)
    #gmsh.option.setNumber("Mesh.Algorithm", 8)
    mesh.generate(2)
    mesh.recombine()


    # visualize & export
    gmsh.fltk.run()
    gmsh.write('./simdata/induction_verification.msh2')
    gmsh.finalize()

    return ph_cylinder, ph_cylinder_ends, ph_coil, ph_air, ph_cylinder_surf, ph_outside_surfs


def sif(ph_cylinder, ph_cylinder_ends, ph_coil, ph_air, ph_cylinder_surf, ph_outside_surfs, omega, I, rho, d, l_tot, with_cylinder=True, mgdyn=False):
    sim = elmer.Simulation()
    sim.settings = {
        "Max Output Level": 4,
        "Coordinate System": "Axi Symmetric",
        "Simulation Type": "Steady state",
        "Steady State Max Iterations": 10,
        "Output File": "case.result"
    }
    sim.settings.update({'Angular Frequency': omega})  # for mgdyn solver
    # solvers
    if mgdyn:
        solver_mgdyn = elmer.Solver(sim, 'mgdyn')
        solver_mgdyn.data = {
            'Equation': 'MgDyn2DHarmonic',
            'Procedure': '"MagnetoDynamics2D" "MagnetoDynamics2DHarmonic"',
            'Variable': 'Potential[Potential Re:1 Potential Im:1]',
            'Variable Dofs': 2
        }
        solver_calcfields = elmer.Solver(sim, 'calcfields')
        solver_calcfields.data = {
            'Equation': 'CalcFields',
            'Procedure': '"MagnetoDynamics" "MagnetoDynamicsCalcFields"',
            'Potential Variable': 'Potential',
            'Calculate Joule Heating': True,
            'Calculate Magnetic Field Strength': True,
            'Calculate Electric Field': True,
            'Calculate Current Density': True,
            'Calculate Nodal Fields': 'Logical False',
            'Calculate Elemental Fields': 'Logical True'
        }
    else:
        solver_statmag = elmer.Solver(sim, 'statmag')
        solver_statmag.data = {
            'Equation': 'Potential Solver',
            'Procedure': '"StatMagSolve" "StatMagSolver"',
            'Variable': 'Potential',
            'Variable DOFs': 2,
            'Angular Frequency': omega,
            'Calculate Joule Heating': 'Logical True',
            'Calculate Magnetic Flux': 'Logical True'
        }

    if with_cylinder:
        solver_heat = elmer.Solver(sim, 'heat')
        solver_heat.data = {
            'Equation': 'Heat Equation',
            'Procedure': '"HeatSolve" "HeatSolver"',
            'Variable': '"Temperature"',
            'Variable Dofs': 1,
            'Calculate Loads': True
        }
        solver_save_scalars = elmer.Solver(sim, 'save_scalars')
        solver_save_scalars.data = {
            'Exec Solver': 'After timestep',
            'Equation': 'SaveScalars',
            'Procedure': '"SaveData" "SaveScalars"',
            'Filename': '"scalars.dat"',
            'Operator 1': 'boundary sum',
            'Variable 1': 'Temperature Loads',
        }
    solver_output = elmer.Solver(sim, 'output')
    solver_output.data = {
        'Exec Solver': 'After timestep',
        'Equation': '"ResultOutput"',
        'Procedure': '"ResultOutputSolve" "ResultOutputSolver"',
        'VTU Format': True,
        'Save Geometry Ids': 'Logical True',
        'Scalar Field 1': '"Temperature"',
        'Scalar Field 2': '"Heat Conductivity"',
        'Scalar Field 3': '"Temperature Loads"',
        'Scalar Field 4': '"Potential"',
        'Scalar Field 5': '"Joule Heating"',
        'Scalar Field 6': '"Magnetic Flux Density"',
        'Scalar Field 7': '"Phase Surface"',
        'Scalar Field 8': '"Joule Heating E"',
        'Scalar Field 9': '"Magnetic Flux Density E"',
    }

    solver_save_line = elmer.Solver(sim, 'save_line')
    solver_save_line.data={
        'Equation': 'SaveLine',
        'Procedure': 'File "SaveData" "SaveLine"',
        'Filename': 'line.dat',
        'Polyline Coordinates(2,2)': '0.0 0.025 0.04 0.025',
        'Polyline Divisions(1)': '200',
        'Exact Coordinates': True,
        'Optimize Node Ordering': 'Logical False',
    }

    # equations
    if mgdyn:
        eqn_mgdyn = elmer.Equation(sim, 'mgdyn', [solver_mgdyn, solver_calcfields])
    else:
        eqn_statmag = elmer.Equation(sim, 'mgdyn', [solver_statmag])
    if with_cylinder:
        if mgdyn:
            eqn_mgdyn_heat = elmer.Equation(sim, 'mgdyn_heat', [solver_mgdyn, solver_calcfields, solver_heat])
        else:
            eqn_statmag_heat = elmer.Equation(sim, 'mgdyn_heat', [solver_statmag, solver_heat])

    # materials
    if with_cylinder:
        material_cylinder = elmer.Material(sim, 'material_cylinder')
        material_cylinder.data = {
            'Density': 1,
            'Electric Conductivity': 1 / rho,
            'Heat Capacity': 1,
            'Heat Conductivity': 100,
            'Relative Permeability': 1,
            'Relative Permittivity': 1,
        }
    material_coil = elmer.Material(sim, 'material_coil')
    material_coil.data = {
        'Density': 1,
        'Electric Conductivity': 0,
        'Relative Permeability': 1,
        'Relative Permittivity': 1
    }
    material_air = elmer.Material(sim, 'material_air')
    material_air.data = {
        'Density': 1,
        'Electric Conductivity': 0,
        'Relative Permeability': 1,
        'Relative Permittivity': 1
    }

    # forces
    current = elmer.BodyForce(sim, 'current')
    current.current_density = I / (d * l_tot)
    if with_cylinder:
        joule_heat = elmer.BodyForce(sim, 'joule_heat')
        joule_heat.joule_heat = True

    # initial condition
    if with_cylinder:
        t0 = elmer.InitialCondition(sim, 'T0')
        t0.data = {'Temperature': 293.15}

    # bodies
    if with_cylinder:
        cylinder = elmer.Body(sim, 'cylinder', [ph_cylinder])
        cylinder.material = material_cylinder
        if mgdyn:
            cylinder.equation = eqn_mgdyn_heat
        else:
            cylinder.equation = eqn_statmag_heat
        cylinder.initial_condition = t0
        cylinder.body_force = joule_heat

        cylinder_ends = elmer.Body(sim, 'cylinder_ends', [ph_cylinder_ends])
        cylinder_ends.material = material_cylinder
        if mgdyn:
            cylinder_ends.equation = eqn_mgdyn
        else:
            cylinder_ends.equation = eqn_statmag

    coil = elmer.Body(sim, 'coil', [ph_coil])
    coil.material = material_coil
    if mgdyn:
        coil.equation = eqn_mgdyn
    else:
        coil.equation = eqn_statmag
    coil.body_force = current

    air = elmer.Body(sim, 'air', [ph_air])
    air.material = material_air
    if mgdyn:
        air.equation = eqn_mgdyn
    else:
        air.equation = eqn_statmag

    # boundaries
    if with_cylinder:
        cylinder_surf = elmer.Boundary(sim, 'cylinder_surf', [ph_cylinder_surf])
        cylinder_surf.fixed_temperature = 293.15
        cylinder_surf.save_scalars = True
    outside_surf = elmer.Boundary(sim, 'outside_surf', [ph_outside_surfs])
    outside_surf.zero_potential = True

    # write sif
    sim.write_sif('./simdata/')


def run_elmer():
    # Elmer Grid
    with open('./simdata/elmergrid.log', 'w') as f:
        args = ['C:/Program Files/Elmer 8.4-Release/bin/ElmerGrid.exe', '14', '2', './simdata/induction_verification.msh2']
        subprocess.run(args, stdout=f, stderr=f)
    files = os.listdir('./simdata/induction_verification/')
    for f in files:
        if os.path.exists('./simdata/' + f):
            os.remove('./simdata/' + f)
        shutil.move('./simdata/induction_verification/' + f, './simdata/')
    shutil.rmtree('./simdata/induction_verification/')

    # Elmer Solver
    with open('./simdata/elmersolver.log', 'w') as f:
        args = ['C:/Program Files/Elmer 8.4-Release/bin/ElmerSolver.exe', './case.sif']
        subprocess.run(args, cwd='./simdata', stdout=f, stderr=f)


if __name__ == "__main__":
    l = 0.05
    r_e = 0.035
    r_i = 0.04
    d = 0.01
    n_l = 10
    l_tot = n_l * l
    omega = 84200
    I = 89.225
    rho = 1.7141032172810484e-05
    if not os.path.exists('./simdata'):
        os.mkdir('./simdata')
    ph_cylinder, ph_cylinder_ends, ph_coil, ph_air, ph_cylinder_surf, ph_outside_surfs = geometry(l, n_l, r_e, r_i, d)
    sif(ph_cylinder, ph_cylinder_ends, ph_coil, ph_air, ph_cylinder_surf, ph_outside_surfs, omega, I, rho, d, l_tot, with_cylinder=True, mgdyn=True)
