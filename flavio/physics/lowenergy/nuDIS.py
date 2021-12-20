r"""Functions for neutrino-electron scattering."""

import flavio
from flavio.physics.zdecays.smeftew import d_gV, d_gA, gV_SM, gA_SM, d_gW, _sinthetahat2
from flavio.physics.wdecays.gammaw import gWf_SM
from flavio.physics.wdecays.mw import dmW_SMEFT
from math import pi, log, sqrt
from wilson import Wilson

def couplings_eff(wc_obj, par, Q2, flav):
    GF, alpha_e = par['GF'], 137
    
    flav = 'mu'
    label = '22'
    
    scale = flavio.config['renormalization scale']['zdecays']
    
    gV_SM_u = gV_SM('u', par)
    gA_SM_u = gA_SM('u', par)
    gV_SM_d = gV_SM('d', par)
    gA_SM_d = gA_SM('d', par)
    gW_SM_q = gWf_SM('u', 'd', par)
        
    gV_SM_nu = gV_SM(flav, par)
    gA_SM_nu = gA_SM(flav, par)
    gW_SM_l = gWf_SM(flav, flav, par)
    
    wc_dict = wc_obj.get_wcxf(sector='all', scale=scale, par=par, eft='SMEFT', basis='Warsaw')
    d_gV_u = d_gV('u', 'u', par, wc_dict)
    d_gA_u = d_gA('u', 'u', par, wc_dict)
    d_gV_d = d_gV('d', 'd', par, wc_dict)
    d_gA_d = d_gA('d', 'd', par, wc_dict)
    d_gV_nu = d_gV(flav, flav, par, wc_dict)
    d_gA_nu = d_gA(flav, flav, par, wc_dict)

    d_gW_q = d_gW('u', 'd', par, wc_dict)
    d_gW_l = d_gW(flav, flav, par, wc_dict)
    
    dmW = dmW_SMEFT(par, wc_dict)
    
    C_lq1 = wc_dict['lq1_' + label + '11']
    C_lq3 = wc_dict['lq3_' + label + '11']
    
    C_lu = wc_dict['lu_' + label + '11']
    C_ld = wc_dict['ld_' + label + '11']

    sh2 = _sinthetahat2(par)
    
    L_u = (1 / (sqrt(2) * GF) * (C_lq1 + C_lq3) - 4 * ((gV_SM_nu + gA_SM_nu) * (gV_SM_u + gA_SM_u)).real - 4 * ((gV_SM_nu + gA_SM_nu) * (d_gV_u + d_gA_u)).real - 4 * ((gV_SM_u + gA_SM_u) * (d_gV_nu + d_gV_nu)).real) / 2
    LR_u = (1 / (sqrt(2) * GF) * C_lu
    - 4 * ((gV_SM_u - gA_SM_u) * (gV_SM_nu + gA_SM_nu)).real - 4 * ((gV_SM_u - gA_SM_u) * (d_gV_nu + d_gA_nu)).real - 4 * ((gV_SM_nu + gA_SM_nu) * (d_gV_u - d_gV_u)).real) / 2

    
    L_d = (1 / (sqrt(2) * GF) * (C_lq1 - C_lq3) - 4 * ((gV_SM_nu + gA_SM_nu) * (gV_SM_d + gA_SM_d)).real - 4 * ((gV_SM_nu + gA_SM_nu) * (d_gV_d + d_gA_d)).real - 4 * ((gV_SM_d + gA_SM_d) * (d_gV_nu + d_gV_nu)).real) / 2
    LR_d = (1 / (sqrt(2) * GF) * C_ld
    - 4 * ((gV_SM_d - gA_SM_d) * (gV_SM_nu + gA_SM_nu)).real - 4 * ((gV_SM_d - gA_SM_d) * (d_gV_nu + d_gA_nu)).real - 4 * ((gV_SM_nu + gA_SM_nu) * (d_gV_d - d_gV_u)).real) / 2

    L_nuedu = -1 * (2 / (sqrt(2) * GF) * C_lq3 - 4 * gW_SM_q * gW_SM_l * (1 - 2 * dmW) - 4 * gW_SM_q * d_gW_l - 4 * gW_SM_l * d_gW_q) / (2 * gW_SM_q) # check, especially the factor 2 in front of dmW
    
    #L_SM_u =  -2 * ((gV_SM_u + gA_SM_u) * (gV_SM_nu + gA_SM_nu)).real
    #LR_SM_u =  -2 * ((gV_SM_u - gA_SM_u) * (gV_SM_nu + gA_SM_nu)).real
    #L_SM_d =  -2 * ((gV_SM_d + gA_SM_d) * (gV_SM_nu + gA_SM_nu)).real
    #LR_SM_d =  -2 * ((gV_SM_d - gA_SM_d) * (gV_SM_nu + gA_SM_nu)).real

        
    w = Wilson({ 'CVLL_numunumuuu': L_u, 'CVLR_numunumuuu': LR_u, 'CVLL_numunumudd': L_d, 'CVLR_numunumudd': LR_d, 'CVL_dumunumu': L_nuedu }, scale=scale, eft='WET', basis='flavio')
    wc_obj_WET = flavio.WilsonCoefficients.from_wilson(w, par)
    wc_dict_WET = wc_obj_WET.get_wc(sector='all', scale=sqrt(Q2), par=par, eft='WET', basis='flavio')
    #print(wc_dict_WET['CVLL_eeee'])
    L_LE_u = -1 * wc_dict_WET['CVLL_numunumuuu']
    LR_LE_u = -1 * wc_dict_WET['CVLR_numunumuuu']
    L_LE_d = -1 * wc_dict_WET['CVLL_numunumudd']
    LR_LE_d = -1 * wc_dict_WET['CVLR_numunumudd']
    L_LE_nuedu = 2 * gW_SM_q * wc_dict_WET['CVLR_numunumudd']
    
    gL2 = (abs(L_LE_u)**2 + abs(L_LE_d)**2) / (abs(L_LE_nuedu)**2)
    gR2 = (abs(LR_LE_u)**2 + abs(LR_LE_d)**2) / (abs(L_LE_nuedu)**2)
    hL2 = (abs(L_LE_u)**2 - abs(L_LE_d)**2) / (abs(L_LE_nuedu)**2)
    hR2 = (abs(LR_LE_u)**2 - abs(LR_LE_d)**2) / (abs(L_LE_nuedu)**2)
    
    return gL2, gR2, hL2, hR2

def CCFR(wc_obj, par, Q2, flav):
    gL2, gR2, hL2, hR2 = couplings_eff(wc_obj, par, Q2, flav)
    kappa = 1.7897 * gL2 + 1.1479 * gR2 - 0.0916 * hL2 - 0.0782 * hR2
    print('With running effect: {}'.format(kappa))

def CHARM_CDHS(wc_obj, par, flav, label):
    if label=='CHARM1': r, Q2 = 0.456, 4**2
    if label=='CHARM2': r, Q2 = 0.429, 9**2
    if label=='CDHS': r, Q2 = 0.393, 10**2
    # check Q2 values
    gL2, gR2, hL2, hR2 = couplings_eff(wc_obj, par, Q2, flav)
    Rnu = gL2 + r * gR2
    Rnuba = gL2 + gR2 / 2
    print('With running effect, {} : Rv={}, Rvbar={}'.format(label, Rnu, Rnuba))


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
