import numpy as np
from matplotlib import pyplot as plt
import pyneb
h1 = pyneb.RecAtom('H', 1)
tems = [2000, 4000, 8000, 12000, 16000, 20000]
den = np.logspace(2.0, 6.0)
linepairs = [
    [(2, 4), (2, 3)],
    [(2, 3), (3, 9)],
    [(3, 10), (3, 9)],
    [(3, 11), (3, 9)],
    [(3, 12), (3, 9)],
    [(3, 13), (3, 9)],
    [(3, 14), (3, 9)],
    [(3, 15), (3, 9)],
    [(3, 16), (3, 9)],
    [(3, 17), (3, 9)],
]
for levelsA, levelsB in linepairs:
    wavA = h1.getWave(*levelsA)
    wavB = h1.getWave(*levelsB)
    fig, ax = plt.subplots(1, 1)
    for T in tems:
        # Level order is backward for emissivity
        emA = h1.getEmissivity(T, den, *levelsA[::-1])
        emB = h1.getEmissivity(T, den, *levelsB[::-1])
        ax.plot(den, emA/emB, label='T = {:.0f} K'.format(T))
    ax.set_xscale('log')
    ax.set_xlim(den.min(), den.max())
    ax.set_ylim(0.0, None)
    ax.set_xlabel('Density, pcc')
    ax.set_ylabel('Line Ratio H({}-{})/H({}-{})'.format(*(levelsA[::-1]+levelsB[::-1])))
    ax.set_title('H I {:.2f} / {:.2f}'.format(wavA, wavB))
    ax.legend(loc='lower right', fontsize='small')
    fig.set_size_inches(7, 7)
    fig.savefig('intrinsic-ratio-H_I_{}_{}.pdf'.format(int(wavA+0.5), int(wavB+0.5)))
