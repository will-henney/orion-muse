from __future__ import print_function
import sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import astropy.constants as C
import astropy.units as u
# import k_B, Ryd, h, c
figname = sys.argv[0].replace('.py', '.pdf')
D = '../Storey-OIII/oiii/'
kT_Joules = 1e4*C.k_B 
Ryd_Joules = C.Ryd*C.h*C.c

# ID, term, stat weight for ground level
lower = [
    ['1', '3P0', 1],
    ['2', '3P1', 3],
    ['3', '3P2', 5],
]

# ID, term for excited levels
upper = [
    ['4', '1D'],
    ['5', '1S'],
]

fig, ax = plt.subplots(1, 1)
for id_upper, term_upper in upper:
    label = '3P - ' + term_upper
    E_stack = []
    sigma_stack = []
    E_common = []
    omegas = []
    for id_lower, term_lower, omega in lower:
        s = '{}_{}'.format(id_lower, id_upper)
        E_ryd, Omega = np.loadtxt(D + 'OMEGA_{}_OIII.dat'.format(s), unpack=True)
        E_Joules = E_ryd*Ryd_Joules
        sigma = 2*np.pi*(C.h/(2*np.pi))**2 * (Omega/omega) / (C.m_e*E_Joules)
        sigma = sigma.to(u.cm*u.cm)
        E_over_kT = E_Joules/kT_Joules
        omegas.append(omega)
        E_stack.append(E_over_kT.value)
        sigma_stack.append(sigma.value)
        E_common.extend(list(E_over_kT.value))
    E_common = np.array(sorted(list(set(E_common))))
    sigma_common = np.zeros_like(E_common)
    omega_sum = 0.0
    for E, sigma, omega in zip(E_stack, sigma_stack, omegas):
        sigma_common += omega*np.interp(E_common, E, sigma)
        omega_sum += omega
    sigma_common /= omega_sum
    ax.plot(E_common, sigma_common, lw=0.6, label=label)
ax.set_xlabel('Electron energy / $k T$')
ax.set_ylabel('Cross section, cm$^2$')
ax.set_xlim(0.0, 20)
ax.set_ylim(3e-18, 3e-15)
ax.set_yscale('log')
ax.legend()
fig.savefig(figname)
print(figname)
