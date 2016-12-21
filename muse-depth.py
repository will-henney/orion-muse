from __future__ import print_function
import sys
import numpy as np
from astropy.io import fits
import pyneb 

# Convert muse line surface brightness to erg/s/cm2/sr
MUSE_CONSTANT = 1e-20 * 0.85 / (0.2 / 206265)**2

ESTEBAN_ABUNDANCES = {
    'H': 12.00,  'He': 10.99, 'C': 8.42, 'N': 7.73, 'O': 8.65,
    'Ne': 8.05, 'S': 7.22, 'Cl': 5.46, 'Ar': 6.62, 'Fe': 6.0
}

ESTEBAN_ABUNDANCES_T2_EQUALS_0 = {
    'H': 12.00,  'He': 10.991, 'C': 8.42, 'N': 7.65, 'O': 8.51,
    'Ne': 7.78, 'S': 7.06, 'Cl': 5.33, 'Ar': 6.5, 'Fe': 5.86
}

LINE_DICT = {
    'N_II-6583': {
        'element': 'N',
        'atom': pyneb.Atom(atom='N2'),
        'wave': '6584A',
        'Te': 'Te',
        'Ne': 'Ne',
    },
    'S_II-6731': {
        'element': 'S',
        'atom': pyneb.Atom(atom='S2'),
        'wave': '6731A',
        'Te': 'Te',
        'Ne': 'Ne',
    },
    'S_III-9069': {
        'element': 'S',
        'atom': pyneb.Atom(atom='S3'),
        'wave': '9069A',
        'Te': 'Te-iii',
        'Ne': 'Ne-iii',
    },
    'Cl_III-5538': {
        'element': 'Cl',
        'atom': pyneb.Atom(atom='Cl3'),
        'wave': '5538A',
        'Te': 'Te-iii',
        'Ne': 'Ne-iii',
    },
    'Cl_II-8579': {
        'element': 'Cl',
        'atom': pyneb.Atom(atom='Cl2'),
        'wave': '8579A',
        'Te': 'Te',
        'Ne': 'Ne',
    },
    'Cl_IV-8046': {
        'element': 'Cl',
        'atom': pyneb.Atom(atom='Cl4'),
        'wave': '8579A',
        'Te': 'Te-iii',
        'Ne': 'Ne-iii',
    },
    'O_III-5007': {
        'element': 'O',
        'atom': pyneb.Atom(atom='O3'),
        'wave': '5007A',
        'Te': 'Te-iii',
        'Ne': 'Ne-iii',
    },
    'O_II-7330': {
        'element': 'O',
        'atom': pyneb.Atom(atom='O2'),
        'wave': '7330A',
        'Te': 'Te',
        'Ne': 'Ne',
    },
    'O_I-6300': {
        'element': 'O',
        'atom': pyneb.Atom(atom='O1'),
        'wave': '6300A',
        'Te': 'Te',
        'Ne': 'Ne',
    },
    'H_I-6563':{
        'element': 'H',
        'atom': pyneb.RecAtom('H', 1),
        'label': '3_2',
        'Te': 'Te-iii',
        'Ne': 'Ne-iii',
    },

}


D_M42_PC = 440.0
AU_cm = 1.49597870691e13
iHDU = 'SCALED'
abundances = ESTEBAN_ABUNDANCES_T2_EQUALS_0

def find_emissivity_map(lineid, suffix):
    d = LINE_DICT[lineid]
    Te = fits.open('LineMaps/muse-derived-' + d['Te'] + suffix)[iHDU].data
    Ne = fits.open('LineMaps/muse-derived-' + d['Ne'] + suffix)[iHDU].data
    abundance = 10**(abundances[d['element']] - 12.0)
    if 'wave' in d:
        e = d['atom'].getEmissivity(Te, Ne, wave=d['wave'], product=False)
    else:
        e = d['atom'].getEmissivity(Te, Ne, label=d['label'], product=False).T
    return abundance*Ne*Ne*e


if __name__ == '__main__':
    try:
        lineid = sys.argv[1]
        suffix = sys.argv[2] + '.fits'
    except IndexError:
        print('Usage: {} LINEID SUFFIX')
  
    hdu = fits.open('LineMaps/linesum-{}-excorr{}'.format(lineid, suffix))[iHDU]
    s = MUSE_CONSTANT*hdu.data
    em = find_emissivity_map(lineid, suffix)
    d_cm = 4.0*np.pi * s / em
    d_arcsec = d_cm / (D_M42_PC*AU_cm)

    fn = 'LineMaps/emissivity-' + lineid + suffix
    fits.PrimaryHDU(data=em,
                    header=hdu.header).writeto(fn, clobber=True)
    print(fn)

    fn = 'LineMaps/depth-' + lineid + suffix
    fits.PrimaryHDU(data=d_arcsec,
                    header=hdu.header).writeto(fn, clobber=True)
    print(fn)
