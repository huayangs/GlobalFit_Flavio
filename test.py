import flavio
import flavio.plots

errors = flavio.sm_error_budget('BR(Bs->mumu)')
flavio.plots.error_budget_pie(errors)
