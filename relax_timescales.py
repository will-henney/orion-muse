from __future__ import print_function
import sys
from scipy.special import erf
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

def erfd(x):
    '''First derivative of error function'''
    return 2*np.exp(-x**2)/np.sqrt(np.pi)


def tau_F(xi):
    '''Velocity-space friction timescale'''
    return xi**2 / (erf(xi) - xi*erfd(xi))


def tau_s(xi):
    '''Slowing down timescale'''
    return xi*tau_F(xi)


def tau_d(xi):
    '''Deflection timescale'''
    return 2*xi**5 / ((2*xi**2 - 1)*erf(xi) + xi*erfd(xi))


if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1)
    xi = np.logspace(-1.5, 1.5, 500)
    ax.plot(xi, tau_F(xi), label=r'$\tau_F$')
    ax.plot(xi, tau_s(xi), label=r'$\tau_s$')
    ax.plot(xi, tau_d(xi), label=r'$\tau_d$')
    ax.set_xlabel(r'$\xi = u\, /  \langle u \rangle$')
    ax.set_ylabel(r'$\tau \, / \, \tau_0 $')
    ax.set_xlim(xi[0], xi[-1])
    ax.set_ylim(1e-3, 1e3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc='lower right', fontsize='large')
    plotfile = sys.argv[0].replace('.py', '.pdf')
    fig.set_size_inches(3.5, 3.5)
    fig.tight_layout()
    fig.savefig(plotfile)
    print(plotfile)
