r"""Functions for electron scattering."""

import flavio
from flavio.physics.zdecays.smeftew import d_gV, d_gA, gV_SM, gA_SM, _sinthetahat2
from math import pi, log, sqrt
from wilson import Wilson

def moller(wc_obj, par):
    GF, alpha_e, Q2 = par['GF'], 137, 0.026
    y = 0.5
    prefactor = (1 - y) / (1 + y**4 + (1 - y)**4)
    
    scale = flavio.config['renormalization scale']['zdecays']
    
    gV_SM_e = gV_SM('e', par)
    gA_SM_e = gA_SM('e', par)
    
    wc_dict = wc_obj.get_wcxf(sector='all', scale=scale, par=par, eft='SMEFT', basis='Warsaw')
    d_gV_e = d_gV('e', 'e', par, wc_dict)
    d_gA_e = d_gA('e', 'e', par, wc_dict)
    C_ll = wc_dict['ll_1111']
    C_ee = wc_dict['ee_1111']
    C_le = wc_dict['le_1111']
    sh2 = _sinthetahat2(par)
    
    L =  (1 / (sqrt(2) * GF) * C_ll - 2 * (gV_SM_e + gA_SM_e).real**2 - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_e + d_gA_e)).real) / 2
    R =  (1 / (sqrt(2) * GF) * C_ee - 2 * (gV_SM_e - gA_SM_e).real**2 - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_e - d_gA_e)).real) / 2
    LR =  (1 / (sqrt(2) * GF) * C_le - 4 * ((gV_SM_e + gA_SM_e) * (gV_SM_e - gA_SM_e)).real - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_e - d_gA_e)).real - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_e + d_gA_e)).real) / 2

    L_SM =  -1 * (gV_SM_e + gA_SM_e).real**2
    R_SM =  -1 * (gV_SM_e - gA_SM_e).real**2
    LR_SM =  -2 * ((gV_SM_e + gA_SM_e) * (gV_SM_e - gA_SM_e)).real
        
    w = Wilson({ 'CVLL_eeee': L, 'CVRR_eeee': R, 'CVLR_eeee': LR }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    dL_LE = wc_dict_WET['CVLL_eeee']
    dR_LE = wc_dict_WET['CVRR_eeee']

    w = Wilson({ 'CVLL_eeee': L_SM, 'CVRR_eeee': R_SM, 'CVLR_eeee': LR_SM }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    dL_LE -= wc_dict_WET['CVLL_eeee']
    dR_LE -= wc_dict_WET['CVRR_eeee']
    
    ALR_mol_NP = -2 * (-2 * dL_LE + 2 * dR_LE) * GF / (sqrt(2) * pi * alpha_e) * prefactor
    
    print('Match SMEFT onto LEFT at m_Z scale by hand')
    print('With running effect: {}'.format(ALR_mol_NP))
    
def moller_simple(wc_obj, par):
    GF, alpha_e, Q2 = par['GF'], 137, 0.026
    scale = flavio.config['renormalization scale']['zdecays']
    wc_dict = wc_obj.get_wcxf(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    dL_LE = wc_dict['CVLL_eeee']
    dR_LE = wc_dict['CVRR_eeee']
    
    y = 0.5
    prefactor = (1 - y) / (1 + y**4 + (1 - y)**4)
    
    ALR_mol_NP = -2 * (-2 * dL_LE + 2 * dR_LE) * GF / (sqrt(2) * pi * alpha_e) * prefactor
    
    print('Automatically caclulation by flavio')
    print(ALR_mol_NP)



'''
# Observable and Prediction instances
_tex = {'e': 'e', 'mu': r'\mu'}

for (l1, l2, l3) in [('mu', 'e', 'e'), ('mu', 'mu', 'mu'),
                     ('e', 'e', 'e'), ('e', 'mu', 'mu'),
                     ('e', 'mu', 'e'), ('mu', 'e', 'mu')]:
    _process_tex = r"\tau^-\to " + _tex[l1] + r"^-" + _tex[l2] + r"^+" + _tex[l3] + r"^-"
    _process_taxonomy = r'Process :: $\tau$ lepton decays :: LFV decays :: $\tau\to \ell^\prime\ell\ell$ :: $' + _process_tex + r"$"

    _obs_name = "BR(tau->" + l1 + l2 + l3 + ")"
    _obs = flavio.classes.Observable(_obs_name)
    _obs.set_description(r"Branching ratio of $" + _process_tex + r"$")
    _obs.tex = r"$\text{BR}(" + _process_tex + r")$"
    _obs.add_taxonomy(_process_taxonomy)
    flavio.classes.Prediction(_obs_name, br_taul1l2l3_fct(l1, l2, l3))
'''
