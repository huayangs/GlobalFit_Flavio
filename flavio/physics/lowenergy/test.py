r"""Functions for electron scattering."""

import flavio
from wilson import Wilson
from moller import moller, moller_simple
from nuelectron import nuelectron
from eDIS import eDIS, Qweak
from nuDIS import CCFR, CHARM_CDHS

'''
print('Only O_ll_{1111}: 0.0001')
w = Wilson({ 'll_1111': 0.0001 }, scale=flavio.config['renormalization scale']['zdecays'], eft='SMEFT', basis='Warsaw')
par=flavio.default_parameters.get_central_all()
wc_obj=flavio.WilsonCoefficients.from_wilson(w, par)
moller(wc_obj, par)
moller_simple(wc_obj, par)

print('O_ll_{1111}: 0.0001, and O_phil1_{11}: 0.0001')
w = Wilson({ 'll_1111': 0.0001 , 'phil1_11': 0.0001 }, scale=flavio.config['renormalization scale']['zdecays'], eft='SMEFT', basis='Warsaw')
par=flavio.default_parameters.get_central_all()
wc_obj=flavio.WilsonCoefficients.from_wilson(w, par)
moller(wc_obj, par)
moller_simple(wc_obj, par)


print('Only O_ll_{1122}: 0.0001')
w = Wilson({ 'll_1122': 0.0001 }, scale=flavio.config['renormalization scale']['zdecays'], eft='SMEFT', basis='Warsaw')
par=flavio.default_parameters.get_central_all()
wc_obj=flavio.WilsonCoefficients.from_wilson(w, par)
nuelectron(wc_obj, par, 0)

print('Only O_lq1_{1111}: 0.0001')
w = Wilson({ 'lq1_1111': 0.0001 }, scale=flavio.config['renormalization scale']['zdecays'], eft='SMEFT', basis='Warsaw')
par=flavio.default_parameters.get_central_all()
wc_obj=flavio.WilsonCoefficients.from_wilson(w, par)
eDIS(wc_obj, par, 0.1)
Qweak(wc_obj, par, 100, 80)
'''

print('Only O_lq1_{1111}: 0.0001')
w = Wilson({ 'lq1_1111': 0.00000001 }, scale=flavio.config['renormalization scale']['zdecays'], eft='SMEFT', basis='Warsaw')
par=flavio.default_parameters.get_central_all()
wc_obj=flavio.WilsonCoefficients.from_wilson(w, par)
CCFR(wc_obj, par, 4**2, 'mu')
CHARM_CDHS(wc_obj, par, 'mu', 'CHARM1')
