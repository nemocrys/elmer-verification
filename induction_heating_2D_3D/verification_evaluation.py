from pyelmer.post import dat_to_dataframe
import matplotlib.pyplot as plt
import yaml
import numpy as np
from scipy.special import jv

plt.rcParams.update({'font.size': 14})

#############################
# input parameters
config_file = "induction_heating_2D_v2/config.yml"
mesh_study = {
    "1 mm": "induction_heating_2D_v2\simdata_air-size=0.8_mesh-size=2",
    "0.5 mm": "induction_heating_2D_v2\simdata_air-size=0.8_mesh-size=1",
    "0.25 mm": "induction_heating_2D_v2\simdata_air-size=0.8_mesh-size=0.5"
}
radius_study = {
    0.4: "induction_heating_2D_v2\simdata_air-size=0.4_mesh-size=1",
    0.8: "induction_heating_2D_v2\simdata_air-size=0.8_mesh-size=1",
    1.6: "induction_heating_2D_v2\simdata_air-size=1.6_mesh-size=1"
}
temperature_analytical = 719.681573588678  # from analytical_solution_temperature.ipynb

with open(config_file) as f:
    config = yaml.safe_load(f)

r_e = config["graphite_r"]  # = crucible radius
l = config["height"]  # = crucible height
I = config["current"]  # A
omega = 2 * np.pi* config["frequency"]  # similar to experiment
sigma = config["graphite"]["Electric Conductivity"]  # S/m


fig, ax = plt.subplots(figsize=(6, 4))

#############################
# compute and plot analytical solution
rho = 1/sigma
H = I/l
delta = (2 * rho / (omega * 4e-7 * np.pi))**0.5
print("delta =", delta)
m = 2**0.5 * r_e / delta
print('m =', m)
print(f"H = {H} A/m")

# Analytical solution according to Lupi2017
def J(xi, H, delta, m):
    J = (-1j)**0.5 * H * 2**0.5 / delta * jv(1, (-1j)**0.5 * m * xi) / jv(0, (-1j)**0.5 * m)
    return J

xi = np.linspace(0, 1, 1000)
J_xi = J(xi, H, delta, m)
w = rho * np.abs(J_xi)**2 / 2

ax.plot(xi * r_e *1000, w/1e6, label="Analytical")

#############################
# compute and plot numerical solution
for name, simdir in mesh_study.items():
    df = dat_to_dataframe(f"{simdir}/results/save_line.dat")
    ax.plot(df["coordinate 1"]*1e3, df["joule heating"]/1e6, label=name)

ax.legend()
ax.set_xlabel( 'Cylinder radius in mm')
ax.set_ylabel('Joule heat in $\\frac{\\mathrm{MW}}{\mathrm{m}^3}$')
ax.grid(linestyle=":")
fig.tight_layout()
fig.savefig("induction_heating_2D_v2/mesh_influence.svg")
plt.close(fig)

### zoom ###
r_zoom = 59.9
plt.rcParams.update({'font.size': 18})
fig, ax = plt.subplots(figsize=(4, 3))
ax.grid(linestyle=":")
radius = np.array(xi * r_e * 1000)
power = np.array(w/1e6)
filter = radius >= r_zoom
line, = ax.plot(radius[filter], power[filter])
line.set_label('analytical')

for name, simdir in mesh_study.items():
    df = dat_to_dataframe(f"{simdir}/results/save_line.dat")
    df_filtered = df.loc[df["coordinate 1"] >= r_zoom *1e-3]
    ax.plot(df_filtered["coordinate 1"]*1e3, df_filtered["joule heating"]/1e6, label=name)
fig.tight_layout()
fig.savefig("induction_heating_2D_v2/mesh_influence_zoom.svg")
plt.close(fig)

#############################
# compare error at r_max for boundary radius
plt.rcParams.update({'font.size': 14})

fig, ax = plt.subplots(figsize=(6, 4))

heating_analytical = max(power)*1e6  # in W
error_heating = []
error_temperature = []
radii = []

for radius, simdir in radius_study.items():
    df = dat_to_dataframe(f"{simdir}/results/save_line.dat")
    df.sort_values("coordinate 1", inplace=True)
    df.reset_index(inplace=True)
    heating = df.loc[len(df)-1, "joule heating"]
    temperature = df.loc[len(df)-1, "temperature"]
    error_heating.append((heating_analytical - heating)/heating_analytical)
    error_temperature.append((temperature_analytical - temperature)/temperature_analytical)
    radii.append(radius)
ax.plot(radii, np.array(error_heating)*100, "x-", label="Maximum joule heat")
ax.plot(radii, np.array(error_temperature)*100, "x-", label="Surface temperature")
ax.grid(linestyle=":")
ax.set_xlabel("Boundary radius in m")
ax.set_ylabel("Releative error in %")
ax.legend()
fig.tight_layout()
fig.savefig("induction_heating_2D_v2/radius_influence.pdf")
print("Error heating", np.array(error_heating)*100)
print("Error temperature", np.array(error_temperature)*100)
