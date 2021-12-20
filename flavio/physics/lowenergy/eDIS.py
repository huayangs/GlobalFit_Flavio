r"""Functions for neutrino-electron scattering."""

import flavio
from flavio.physics.zdecays.smeftew import d_gV, d_gA, gV_SM, gA_SM, _sinthetahat2
from math import pi, log, sqrt
from wilson import Wilson

def d_gAV(wc_obj, par, Q2):
    GF, alpha_e = par['GF'], 137
    
    scale = flavio.config['renormalization scale']['zdecays']
    
    gV_SM_e = gV_SM('e', par)
    gA_SM_e = gA_SM('e', par)
    
    gV_SM_u = gV_SM('u', par)
    gA_SM_u = gA_SM('u', par)
    
    gV_SM_d = gV_SM('d', par)
    gA_SM_d = gA_SM('d', par)
    
    wc_dict = wc_obj.get_wcxf(sector='all', scale=scale, par=par, eft='SMEFT', basis='Warsaw')
    d_gV_e = d_gV('e', 'e', par, wc_dict)
    d_gA_e = d_gA('e', 'e', par, wc_dict)
    
    d_gV_u = d_gV('u', 'u', par, wc_dict)
    d_gA_u = d_gA('u', 'u', par, wc_dict)
    
    d_gV_d = d_gV('d', 'd', par, wc_dict)
    d_gA_d = d_gA('d', 'd', par, wc_dict)

    
    C_lq1 = wc_dict['lq1_1111']
    C_lq3 = wc_dict['lq3_1111']
    C_eu = wc_dict['eu_1111']
    C_ed = wc_dict['ed_1111']
    C_lu = wc_dict['lu_1111']
    C_ld = wc_dict['ld_1111']
    C_qe = wc_dict['qe_1111']

    sh2 = _sinthetahat2(par)
    
    L_u =  (1 / (sqrt(2) * GF) * (C_lq1 - C_lq3) - 4 * ((gV_SM_e + gA_SM_e) * (gV_SM_u + gA_SM_u)).real - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_u + d_gA_u)).real - 4 * ((gV_SM_u + gA_SM_u) * (d_gV_e + d_gV_e)).real) / 2
    R_u =  (1 / (sqrt(2) * GF) * C_eu - 4 * ((gV_SM_e - gA_SM_e) * (gV_SM_u - gA_SM_u)).real - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_u - d_gA_u)).real - 4 * ((gV_SM_u - gA_SM_u) * (d_gV_e - d_gV_e)).real) / 2
    LR_u =  (1 / (sqrt(2) * GF) * C_eu - 4 * ((gV_SM_e + gA_SM_e) * (gV_SM_u - gA_SM_u)).real - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_u - d_gA_u)).real - 4 * ((gV_SM_u - gA_SM_u) * (d_gV_e + d_gV_e)).real) / 2
    RL_u =  (1 / (sqrt(2) * GF) * C_qe - 4 * ((gV_SM_e - gA_SM_e) * (gV_SM_u + gA_SM_u)).real - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_u + d_gA_u)).real - 4 * ((gV_SM_u + gA_SM_u) * (d_gV_e - d_gV_e)).real) / 2

    L_d =  (1 / (sqrt(2) * GF) * (C_lq1 + C_lq3) - 4 * ((gV_SM_e + gA_SM_e) * (gV_SM_d + gA_SM_d)).real - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_d + d_gA_d)).real - 4 * ((gV_SM_d + gA_SM_d) * (d_gV_e + d_gV_e)).real) / 2
    R_d =  (1 / (sqrt(2) * GF) * C_ed - 4 * ((gV_SM_e - gA_SM_e) * (gV_SM_d - gA_SM_d)).real - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_d - d_gA_d)).real - 4 * ((gV_SM_d - gA_SM_d) * (d_gV_e - d_gV_e)).real) / 2
    LR_d =  (1 / (sqrt(2) * GF) * C_ld - 4 * ((gV_SM_e + gA_SM_e) * (gV_SM_d - gA_SM_d)).real - 4 * ((gV_SM_e + gA_SM_e) * (d_gV_d - d_gA_d)).real - 4 * ((gV_SM_d - gA_SM_d) * (d_gV_e + d_gV_e)).real) / 2
    RL_d =  (1 / (sqrt(2) * GF) * C_qe - 4 * ((gV_SM_e - gA_SM_e) * (gV_SM_d + gA_SM_d)).real - 4 * ((gV_SM_e - gA_SM_e) * (d_gV_d + d_gA_d)).real - 4 * ((gV_SM_d + gA_SM_d) * (d_gV_e - d_gV_e)).real) / 2


    L_SM_u =  -2 * ((gV_SM_e + gA_SM_e) * (gV_SM_u + gA_SM_u)).real
    R_SM_u =  -2 * ((gV_SM_e - gA_SM_e) * (gV_SM_u - gA_SM_u)).real
    LR_SM_u =  -2 * ((gV_SM_e + gA_SM_e) * (gV_SM_u - gA_SM_u)).real
    RL_SM_u =  -2 * ((gV_SM_e - gA_SM_e) * (gV_SM_u + gA_SM_u)).real
    
    L_SM_d =  -2 * ((gV_SM_e + gA_SM_e) * (gV_SM_d + gA_SM_d)).real
    R_SM_d =  -2 * ((gV_SM_e - gA_SM_e) * (gV_SM_d - gA_SM_d)).real
    LR_SM_d =  -2 * ((gV_SM_e + gA_SM_e) * (gV_SM_d - gA_SM_d)).real
    RL_SM_d =  -2 * ((gV_SM_e - gA_SM_e) * (gV_SM_d + gA_SM_d)).real
        
    w = Wilson({ 'CVLL_eeuu': L_u, 'CVRR_eeuu': R_u, 'CVLR_eeuu': LR_u, 'CVLR_uuee': RL_u, 'CVLL_eedd': L_d, 'CVRR_eedd': R_d, 'CVLR_eedd': LR_d, 'CVLR_ddee': RL_d }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    dL_LE_u = wc_dict_WET['CVLL_eeuu']
    dR_LE_u = wc_dict_WET['CVRR_eeuu']
    dLR_LE_u = wc_dict_WET['CVLR_eeuu']
    dRL_LE_u = wc_dict_WET['CVLR_uuee']
    dL_LE_d = wc_dict_WET['CVLL_eedd']
    dR_LE_d = wc_dict_WET['CVRR_eedd']
    dLR_LE_d = wc_dict_WET['CVLR_eedd']
    dRL_LE_d = wc_dict_WET['CVLR_ddee']

    w = Wilson({ 'CVLL_eeuu': L_SM_u, 'CVRR_eeuu': R_SM_u, 'CVLR_eeuu': LR_SM_u, 'CVLR_uuee': RL_SM_u, 'CVLL_eedd': L_SM_d, 'CVRR_eedd': R_SM_d, 'CVLR_eedd': LR_SM_d, 'CVLR_ddee': RL_SM_d }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=1, par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    dL_LE_u -= wc_dict_WET['CVLL_eeuu']
    dR_LE_u -= wc_dict_WET['CVRR_eeuu']
    dLR_LE_u -= wc_dict_WET['CVLR_eeuu']
    dRL_LE_u -= wc_dict_WET['CVLR_uuee']
    dL_LE_d -= wc_dict_WET['CVLL_eedd']
    dR_LE_d -= wc_dict_WET['CVRR_eedd']
    dLR_LE_d -= wc_dict_WET['CVLR_eedd']
    dRL_LE_d -= wc_dict_WET['CVLR_ddee']
    
    d_gAV_u, d_gAV_d, d_gVA_u, d_gVA_d = (dR_LE_u + dRL_LE_u - dL_LE_u - dLR_LE_u), (dR_LE_d + dRL_LE_d - dL_LE_d - dLR_LE_d), (dR_LE_u + dLR_LE_u - dL_LE_u - dRL_LE_u), (dR_LE_d + dLR_LE_d - dL_LE_d - dRL_LE_d)

    return d_gAV_u, d_gAV_d, d_gVA_u, d_gVA_d
    
def eDIS(wc_obj, par, Q2):
    GF, alpha_e = par['GF'], 137
    d_gAV_u, d_gAV_d, d_gVA_u, d_gVA_d = d_gAV(wc_obj, par, Q2)
    
    a1, a2 = d_gAV_u - d_gAV_d / 2, d_gVA_u - d_gVA_u / 2
    
    a1 *= 3 * GF / (5 * sqrt(2) * pi * alpha_e)
    a2 *= 3 * GF / (5 * sqrt(2) * pi * alpha_e)

    print('Match SMEFT onto LEFT at m_Z scale by hand')
    print('With running effect: {}, {}'.format(a1, a2))
   
def Qweak(wc_obj, par, Z, N):
    GF, alpha_e = par['GF'], 137
    d_gAV_u, d_gAV_d, d_gVA_u, d_gVA_d = d_gAV(wc_obj, par, 1)
    
    d_gAV_p, d_gAV_n = 2 * d_gAV_u + d_gAV_d, d_gAV_u + 2 * d_gAV_d
    
    return -2 * (Z * (d_gAV_p + 0.00005) + N * (d_gAV_n + 0.00006)) * (1 - alpha_e / (2 * pi))


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
