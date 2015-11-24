from scipy.special import gamma
import numpy as np
from numpy import exp, sqrt

def A_kappa(kappa):
    return gamma(kappa+1)/gamma(kappa-0.5)/(kappa-1.5)**1.5


def f_M(E):
    return sqrt(E) * exp(-E)


def f_CH(E, a, b):
    return sqrt(E) * (exp(-E) + (b/a**1.5)*exp(-E/a))/(1 + b)

def f_kappa(E, kappa):
    return A_kappa(kappa) * sqrt(E) / (1 + E/(kappa - 1.5))**(kappa + 1)
