class StaticTDDS4State(object):
    def __init__(self,microscopic_rate_calculator = None):
        self.__microscopic_rate_calculator = None
        self.__changed = False

        if microscopic_rate_calculator != None: self.microscopic_rate_calculator = microscopic_rate_calculator

    def update(self):
        if self.microscopic_rate_calculator == None:
            raise RuntimeError("No microscopic rate calculator set.")
        changed = self.microscopic_rate_calculator.update() or self.__changed

        if changed:
            tau = self.microscopic_rate_calculator.timeconstants

            self.tauc =  tau["1->2'"] + tau["2'->2"]*(1+tau["1->2'"]/tau["2'->1"])
            self.taue1 = tau["2->2'"] + tau["2'->1"]*(1+tau["2->2'"]/tau["2'->2"])
            self.taue2 = tau["2->1'"] + tau["1'->1"]*(1+tau["2->1'"]/tau["1'->2"])
            self.taue = ((self.taue1)**-1+(self.taue2)**-1)**-1
            self.__changed = False
            return True #My State has changed
        return False #Nothing happened

    def get_microscopic_rate_calculator(self):
        return self.__microscopic_rate_calculator
    def set_microscopic_rate_calculator(self,c):
        self.__microscopic_rate_calculator = c
        self.__changed = True
    microscopic_rate_calculator = property(get_microscopic_rate_calculator,set_microscopic_rate_calculator)
