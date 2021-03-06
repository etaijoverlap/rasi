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
"""
The rates module collects classes that implement different defect models. Currently, 
there is only the *RateCalculator4State*, which implements the current "workhorse", the
four state defect. The purpose of the calculators in this module is to collect together
the rates for all microscopic transitions which are then used to compute measurable
quantities or the transient behavior of defects.
"""

from rasi.base import BasicCalculator

class RateCalculator4State(BasicCalculator):
    """
        Calculates the microscopic rates for a four state defect. 

        **Input**

        :var recharge_primary: NMP rate calculator for the primary configuration
        :var recharge_secondary: NMP rate calculator for the secondary configuration
        :var nu: Attempt frequency
        :var energies: Energy barriers
        :var temperature: Temperature

        **Output**

        :var rates: A dictionary containing the rates of the four state defect model
        :var timeconstants: A dictionary containing the time constants (inverse of the rates)

    """
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

