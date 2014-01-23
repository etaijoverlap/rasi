###############################################################################
#
#  RASI ... Reliability Analysis of Semiconductor Interfaces
#
#  (c) 2013-2014 Franz Schanovsky 
#  
#  This project is dedicated to the loving memory of Margarete and Johann 
#  Mittermayr.
#
###############################################################################
#
#    This software is licensed under the EUPL V 1.1
#
#    This software is provided "as is" without warranty of any kind, see the 
#    respective section in the EUPL. USE AT YOUR OWN RISK.
#
###############################################################################


from rasi.base import BasicCalculator

class RateCalculator4State(BasicCalculator):
    def __init__(self, **kwargs):
        self.init_input_variables(
                recharge_primary = None,
                recharge_secondary = None,
                nu = 0.,
                energies = {},
                temperature = None
                )

        self.init_output_variables(
                rates = {},
                timeconstants = {}
                )

        self.set_variables(kwargs)

    def do_update(self):
        from scipy.constants import k as kB
        from scipy import exp
        changed = False
        rates = self.rates
        if self.recharge_primary.update() or self.changed_recharge_primary:
            rates["1->2'"] = self.recharge_primary.oxidation_rate
            rates["2'->1"] = self.recharge_primary.reduction_rate
            changed = True
        if self.recharge_secondary.update() or self.changed_recharge_secondary:
            rates["1'->2"] = self.recharge_secondary.oxidation_rate
            rates["2->1'"] = self.recharge_secondary.reduction_rate
            changed = True
        if self.changed_energies or self.changed_nu or self.changed_temperature:
            E = self.energies ; nu = self.nu ; T = self.temperature
            for transition,Eb in E.iteritems():
                rates[transition] = nu * exp( -Eb/(kB*T) )
            changed = True

        if changed:
            for key,rate in rates.iteritems():
                self.timeconstants[key] = rate**-1
        return changed

