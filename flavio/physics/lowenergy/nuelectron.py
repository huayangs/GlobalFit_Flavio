r"""Functions for neutrino-electron scattering."""

import flavio
from flavio.physics.zdecays.smeftew import d_gV, d_gA, gV_SM, gA_SM, _sinthetahat2
from math import pi, log, sqrt
from wilson import Wilson

def nuelectron(wc_obj, par, Q2):
    GF, alpha_e = par['GF'], 137
    
    scale = flavio.config['renormalization scale']['zdecays']
    
    gV_SM_e = gV_SM('e', par)
    gA_SM_e = gA_SM('e', par)
    
    gV_SM_nu = gV_SM('numu', par)
    gA_SM_nu = gA_SM('numu', par)

    
    wc_dict = wc_obj.get_wcxf(sector='all', scale=scale, par=par, eft='SMEFT', basis='Warsaw')
    d_gV_e = d_gV('e', 'e', par, wc_dict)
    d_gA_e = d_gA('e', 'e', par, wc_dict)
    d_gV_nu = d_gV('numu', 'numu', par, wc_dict)
    d_gA_nu = d_gA('numu', 'numu', par, wc_dict)
    C_ll = wc_dict['ll_1122']
    C_le = wc_dict['le_2211']
    sh2 = _sinthetahat2(par)
    
    L =  (2 / (sqrt(2) * GF) * C_ll - 4 * ((gV_SM_e + gA_SM_e) * (gV_SM_nu + gA_SM_nu)).real - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_nu + d_gA_nu)).real - 4 * ((gV_SM_nu + gA_SM_nu) * (d_gV_e + d_gV_e)).real) / 2 # the factor in first parentheses is 1 or 2?
    LR =  (1 / (sqrt(2) * GF) * C_le
    - 4 * ((gV_SM_e - gA_SM_e) * (gV_SM_nu + gA_SM_nu)).real - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_nu + d_gA_nu)).real - 4 * ((gV_SM_nu + gA_SM_nu) * (d_gV_e - d_gV_e)).real) / 2
    
    
    L_SM =  -2 * ((gV_SM_e + gA_SM_e) * (gV_SM_nu + gA_SM_nu)).real
    LR_SM =  -2 * ((gV_SM_e - gA_SM_e) * (gV_SM_nu + gA_SM_nu)).real

        
    w = Wilson({ 'CVLL_numunumuee': L, 'CVLR_numunumuee': LR }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    dL_LE = wc_dict_WET['CVLL_numunumuee']
    dLR_LE = wc_dict_WET['CVLR_numunumuee']

    w = Wilson({ 'CVLL_numunumuee': L_SM, 'CVLR_numunumuee': LR_SM }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    dL_LE -= wc_dict_WET['CVLL_numunumuee']
    dLR_LE -= wc_dict_WET['CVLR_numunumuee']
    
    g_V_NP, g_A_NP = -1 * (dL_LE + dLR_LE), -1 * (dL_LE - dLR_LE)

    print('Match SMEFT onto LEFT at m_Z scale by hand')
    print('With running effect: {}, {}'.format(g_V_NP, g_A_NP))
    
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
