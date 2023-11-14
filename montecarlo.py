import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

plt.style.use('informes')
np.random.seed(42)


# DEFINICION DE VARIABLES

l = 39.5 * 1e-2 # m
sigma_l = 0.1/4 * 1e-2 # m

L = 43 * 1e-2  # m
sigma_L = 0.1/4 * 1e-2  # m

k_n = 1.875


M = 82.9345 * 1e-3  # kg
sigma_M = 0.0001 * 1e-3  # kg

d = 5.1 * 1e-3  # m
sigma_d = 0.1/4 * 1e-3 # m

f = 15.2  # Hz
sigma_f = 0.01667146597566566  # Hz

alpha = 0  # rad/s
sigma_alpha = 0  # rad/s


# GENERAR DATOS

N = 1000000

l = norm.rvs(loc=l, scale=sigma_l, size=N)
L = norm.rvs(loc=L, scale=sigma_L, size=N)
M = norm.rvs(loc=M, scale=sigma_M, size=N)
d = norm.rvs(loc=d, scale=sigma_d, size=N)
f = norm.rvs(loc=f, scale=sigma_f, size=N)

E_i = 64*M*l**4*(alpha**2 + 4*np.pi**2*f**2)/(np.pi*L*d**4*k_n**4) * 1e-9
E = np.mean(E_i)
sigma_E = np.sqrt(np.sum((E - E_i)**2) / (N**2 - N))

print(E, sigma_E)


# GRAFICAR

fig, ax = plt.subplots()

ax.hist(E_i, density=True)

plt.show()