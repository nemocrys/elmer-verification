## Verification of impedance boundary condition (work in progress)

The simulation setup is based on https://github.com/ElmerCSC/elmerfem/tree/devel/fem/tests/mgdyn_harmonic_wire_impedanceBC_circuit

The analytical solution is based on *E.J. Davies, Conduction and Induction Heating, IEE Power Engineering Series 11, 1990*.

### Simulations

- current-driven-wire_impedanceBC: fixed current density at wire end, impedance BC
- current-driven-wire_no_impedanceBC: fixed current density at wire end, no impedance BC (resolved skin layer)
- wire_impedanceBC: fixed potential at wire end, no impedance BC (resolved skin layer), no circuits model
- wire_impedanceBC_circuits: fixed potential at wire end, no impedance BC (resolved skin layer), circuits model
- wire_no_impedanceBC: fixed potential at wire end, no impedance BC (resolved skin layer)
- wire_no_impedanceBC_circuits: fixed potential at wire end, no impedance BC (resolved skin layer), circuits model
