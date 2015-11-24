import pyneb
sii = pyneb.Atom('S', 2)
nii = pyneb.Atom('N', 2)
siii = pyneb.Atom('S', 3)
cliii = pyneb.Atom('Cl', 3)

def rsii_T_den(T, den):
    A = sii.getEmissivity(T, den, wave=6716)
    B = sii.getEmissivity(T, den, wave=6731)
    return A/B

def rnii_T_den(T, den):
    A = nii.getEmissivity(T, den, wave=5755)
    B = nii.getEmissivity(T, den, wave=6583)
    return A/B

def rsiii_T_den(T, den):
    A = siii.getEmissivity(T, den, wave=6312)
    B = siii.getEmissivity(T, den, wave=9069)
    return A/B

def rcliii_T_den(T, den):
    A = cliii.getEmissivity(T, den, wave=5538)
    B = cliii.getEmissivity(T, den, wave=5518)
    return A/B

# TODO: add in the 2-phase calculation
