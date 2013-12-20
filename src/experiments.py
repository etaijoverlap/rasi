from base import Calculator


class StaticTDDS4State(object):
    def __init__(self,microscopic_rate_calculator = None):
        Calculator.__init__(self)
        self.microscopic_rate_calculator = microscopic_rate_calculator

    def update(self):
        if self.microscopic_rate_calculator == None:
            raise RuntimeError("No microscopic rate calculator set.")
        changed = self.microscopic_rate_calculator.update()

        if changed:
            tau = self.microscopic_rate_calculator.timeconstants

            self.tauc =  tau["1->2'"] + tau["2'->2"]*(1+tau["1->2'"]/tau["2'->1"])
            self.taue1 = tau["2->2'"] + tau["2'->1"]*(1+tau["2->2'"]/tau["2'->2"])
            self.taue2 = tau["2->1'"] + tau["1'->1"]*(1+tau["2->1'"]/tau["1'->2"])
            self.taue = ((self.taue1)**-1+(self.taue2)**-1)**-1
            return True #My State has changed
        return False #Nothing happened

