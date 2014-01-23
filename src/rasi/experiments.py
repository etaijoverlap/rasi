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

class StaticTDDS4State(BasicCalculator):
    def __init__(self,**kwargs):
        self.init_variables(
                inputs = {
                    "microscopic_rate_calculator": None
                },
                outputs = {
                    "tauc":None,
                    "taue":None,
                    "taue1":None,
                    "taue2":None
                })
        self.set_variables(kwargs)

    def do_update(self):
        if self.microscopic_rate_calculator == None:
            raise RuntimeError("No microscopic rate calculator set.")
        changed = self.microscopic_rate_calculator.update() or self.changed

        if changed:
            tau = self.microscopic_rate_calculator.timeconstants

            self.internal_tauc =  tau["1->2'"] + tau["2'->2"]*(1+tau["1->2'"]/tau["2'->1"])
            self.internal_taue1 = tau["2->2'"] + tau["2'->1"]*(1+tau["2->2'"]/tau["2'->2"])
            self.internal_taue2 = tau["2->1'"] + tau["1'->1"]*(1+tau["2->1'"]/tau["1'->2"])
            self.internal_taue = ((self.taue1)**-1+(self.taue2)**-1)**-1
            return True #My State has changed
        return False #Nothing happened

