from sympy import init_printing, symbols, Eq, solve, pprint, pi


init_printing() 

a1, a2 = symbols("a1, a2")  # surface area
eps_e1, eps_e2 = symbols("epsilon_e1, epsilon_e2")  # emission for outgoing radiation
eps_r1, eps_r2 = symbols("epsilon_r1, epsilon_r2")  # emission for reflection / absorption
eps_1, eps_2 = symbols("epsilon_1, epsilon_2")  # emission for old case
sigma = symbols("sigma")  # Stefan-Boltzmann constant
t1, t2 = symbols("T_1, T_2")  # temperatures
r1, r2 = symbols("r1, r2")  # radiosities
p = symbols("p")  # power

# eq1 = Eq(a1 / (1-eps_1)* (sigma*t1**4 - r1), p)  # old case
eq1 = Eq(a1 / (1-eps_r1)* (sigma*eps_e1*t1**4 - eps_r1 *r1), p)  # new case

r1 = solve([eq1], [r1])[r1]

eq2 = Eq(a1*(r1 - r2), p)
r2 = solve([eq2], [r2])[r2]

# eq3 = Eq(a1 / (1-eps_1)* (sigma*t1**4 - r1), -a2 / (1-eps_2)* (sigma*t2**4 - r2))  # old case
eq3 = Eq(a1 / (1-eps_r1)* (sigma*eps_e1*t1**4 - eps_r1 *r1), -a2 / (1-eps_r2)* (sigma*eps_e2*t2**4 - eps_r2 *r2))  # new case

sol = solve([eq3], [t1])
pprint(sol[-1])
print(sol[-1])
t1_val = sol[-1][0].subs([(eps_e1, 0.8), (eps_r1, 0.8), (eps_e2, 0.5), (eps_r2, 0.5), (sigma, 5.670374419e-8), (p, 30000), (a1, 4*pi* 0.5**2), (a2, 4*pi* 0.9**2), (t2, 1081.7104112004408)])
print(t1_val.evalf())
