r"""$Z$ partial widths in the SM.

Based on arXiv:1906.08815"""


from math import log
from scipy import constants
import flavio


# units: GeV=hbar=c=1
GeV = constants.giga * constants.eV
s = GeV / constants.hbar
m = s / constants.c
b = 1.e-28 * m**2
pb = constants.pico * b


# Table 1 of 1906.08815
adict = {
'Gammae,mu': [83.983, -0.1202, -0.06919, 0.00383, 0.0597, 0.8037, -0.015, -0.0195, 0.0032, -0.0956, -0.0078, -0.0095, 0.25, -1.08, 0.056, -0.37, 286],
'Gammatau': [83.793, -0.1200, -0.06905, 0.00382, 0.0596, 0.8023, -0.015, -0.0195, 0.0032, -0.0954, -0.0078, -0.0094, 0.25, -1.08, 0.056, -0.37, 285],
'Gammanu': [167.176, -0.1752, -0.1249, 0.00595, 0.1046, 1.253, -0.110, -0.0232, 0.0064, -0.187, -0.014, -0.014, 0.37, -0.085, 0.054, -0.30, 503],
'Gammau': [299.994, -0.6152, -0.2771, 0.0174, 0.2341, 4.051, -0.467, -0.0676, 0.017, 14.26, 1.6, -0.046, 1.82, -11.1, 0.16, -1.0, 1253],
'Gammac': [299.918, -0.6152, -0.2771, 0.0174, 0.2340, 4.051, -0.467, -0.0676, 0.017, 14.26, 1.6, -0.046, 1.82, -11.1, 0.16, -1.0, 1252],
'Gammad,s': [382.829, -0.6685, -0.3322, 0.0193, 0.2792, 3.792, -0.18, -0.0706, 0.020, 10.20, -2.4, -0.052, 0.71, -10.1, 0.16, -0.92, 1469],
'Gammab': [375.890, -0.6017, -0.3158, 0.0190, 0.227, -2.174, 0.042, -0.027, 0.021, 10.53, -2.4, -0.056, 1.2, -10.1, 0.15, -0.95, 1458],
'GammaZ': [2494.75, -4.055, -2.117, 0.122, 1.746, 19.68, -1.63, -0.432, 0.12, 58.61, -4.0, -0.32, 8.1, -56.1, 1.1, -6.8, 9267],
'Rl': [20751.6, -8.112, -1.174, 0.155, 0.16, -37.59, -10.9, 1.27, 0.29, 732.30, -44, -0.61, 5.7, -358, -4.7, 37, 11649],
'Rc': [17222.2, -4.049, -0.749, 0.0832, 1.08, 98.956, -15.1, -0.761, 0.080, 230.9, 125, 0.045, 36.9, -120, 1.2, -6.2, 3667],
'Rb': [21585.0, 4.904, 0.9149, -0.0535, -2.676, -292.21, 20.0, 1.97, -0.11, -131.9, -84, -0.27, 4.4, 71.9, -0.77, -4.4, -1790],
'sigma0had': [41489.6, 0.408, -0.320, 0.0424, 1.32, 60.17, 16.3, -2.31, -0.19, -579.58, 38, 0.010, 7.5, 85.2, 9.1, -68, -85957],
}

# Table 3 of 1906.08815
ddict = {
'sinthetaeffl2': [2314.64, 4.616, 0.539, -0.0737, 206, -25.71, 4.00, 0.288, 3.88, -6.49, -6560],
'sinthetaeffb2': [2327.04, 4.638, 0.558, -0.0700, 207, -9.554, 3.83, 0.179, 2.41, -8.24, -6630],
}

# Converting the table to appropriate powers of GeV
units = {
'Gammae,mu': 1e-3, # MeV -> GeV
'Gammatau': 1e-3, # MeV -> GeV
'Gammanu': 1e-3, # MeV -> GeV
'Gammau': 1e-3, # MeV -> GeV
'Gammac': 1e-3, # MeV -> GeV
'Gammad,s': 1e-3, # MeV -> GeV
'Gammab': 1e-3, # MeV -> GeV
'GammaZ': 1e-3, # MeV -> GeV
'Rl': 1e-3,
'Rc': 1e-5,
'Rb': 1e-5,
'sigma0had': pb, # pb
'sinthetaeffl2': 1e4,
'sinthetaeffb2': 1e4,
}

def Zobs(name, m_h, m_t, alpha_s, Dalpha, m_Z):
    r"""Expansion formula for $Z$ partial widths according to eq. (2.9) of
    arXiv:1906.08815.
    """
    flavio.citations.register("Dubovyk:2019szj")
    L_H = log(m_h / 125.7)
    D_H = m_h / 125.7 - 1
    D_t = (m_t / 173.2)**2 - 1
    D_alpha_s = alpha_s / 0.1184 - 1
    D_alpha = Dalpha / 0.059 - 1
    D_Z = m_Z / 91.1876 - 1
    a = adict[name]
    return (a[0] + a[1] * L_H + a[2] * L_H**2 + a[3] * L_H**4 + a[4] * D_H + a[5] * D_t + a[6] * D_t**2 + a[7] * D_t * L_H + a[8] * D_t * L_H**2 + a[9] * D_alpha_s + a[10] * D_alpha_s**2 + a[11] * D_alpha_s * D_H + a[12] * D_alpha_s * D_t + a[13] * D_alpha + a[14] * D_alpha * D_H + a[15] * D_alpha * D_t + a[16] * D_Z) * units[name]
    
def GammaZ_SM(par, f):
    if f in ['e', 'mu']:
        name = 'Gammae,mu'
    elif f in ['d', 's']:
        name = 'Gammad,s'
    elif 'nu' in f:
        name = 'Gammanu'
    else:
        name = 'Gamma' + f
    GSM = Zobs(name, par['m_h'], par['m_t'], par['alpha_s'], 0.059, par['m_Z'])
    return GSM + par['delta_' + name]
    
def sinthetaeff2(name, m_h, m_t, alpha_s, Dalpha, m_Z):
    r"""Expansion formula for $Z$ partial widths according to eq. (3.6) of
    arXiv:1906.08815.
    """
    flavio.citations.register("Dubovyk:2019szj")
    L_H = log(m_h / 125.7)
    D_H = m_h / 125.7 - 1
    D_t = (m_t / 173.2)**2 - 1
    D_alpha_s = alpha_s / 0.1184 - 1
    D_alpha = Dalpha / 0.059 - 1
    D_Z = m_Z / 91.1876 - 1
    d = ddict[name]
    return (d[0] + d[1] * L_H + d[2] * L_H**2 + d[3] * L_H**4 + d[4] * D_alpha + d[5] * D_t + d[6] * D_t**2 + d[7] * D_t * L_H + d[8] * D_alpha_s + d[9] * D_alpha_s * D_t + a[10] * D_Z) * units[name]

def Af_SM(par, f):
    if f in ['e', 'mu', 'tau']:
        name = 'sinthetaeffl2'
        Qf = 1
    elif f in ['u', 'c']:
        name = 'sinthetaeffl2'
        Qf = 2/3
    elif f in ['d', 's']:
        name = 'sinthetaeffl2'
        Qf = 1/3
    elif 'nu' in f:
        name = 'sinthetaeffl2'
        Qf = 0
    else:
        name = 'sinthetaeffb2'
        Qf = 1/3
    s2w_eff = sinthetaeff2(name, par['m_h'], par['m_t'], par['alpha_s'], 0.059, par['m_Z'])
    return (1 - 4 * Qf * s2w_eff) / (1 - 4 * Qf * s2w_eff + 8 * (Qf * s2w_eff)**2)

def AFBf_SM(par, f):
    A = Af_SM(par, f)
    Ae = Af_SM(par, 'e')
    return 3 / 4 * Ae * A
