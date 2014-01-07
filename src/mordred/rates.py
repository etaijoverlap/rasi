class RateCalculator4State(object):
    def __init__(self):
        self.recharge_primary = None
        self.recharge_secondary = None
        self.__nu = 0.
        self.__energies = {}
        self.__temperature = 1.
        self.__changed = True
        self.rates = {}
        self.timeconstants = {}

    def update(self):
        from scipy.constants import k as kB
        changed = False
        rates = self.rates
        if self.recharge_primary.update():
            rates["1->2'"] = self.recharge_primary.oxidation_rate
            rates["2'->1"] = self.recharge_primary.reduction_rate
            changed = True
        if self.recharge_secondary.update():
            rates["1'->2"] = self.recharge_secondary.oxidation_rate
            rates["2->1'"] = self.recharge_secondary.reduction_rate
            changed = True
        if self.__changed:
            E = self.__energies ; nu = self.__nu ; T = self.__temperature
            for transition,Eb in self.__energies:
                rates[transition] = nu * exp( -Eb[transition]/(kB*T) )
            changed = True

        if changed:
            for key,rate in rates.iteritems():
                self.timeconstants[key] = rate**-1
        self.__changed = False
        return changed

    def get_nu(self):
        return self.__nu
    def set_nu(self,nu):
        self.__nu = nu
        self.__changed = True
    nu = property(get_nu,set_nu)

    def get_energies(self):
        return self.__energies
    def set_energies(self,energies):
        self.__energies = energies
        self.__changed = True
    energies = property(get_energies,set_energies)

    def get_temperature(self):
        return self.__energies
    def set_temperature(self,temperature):
        self.__temperature = temperature
        self.__changed = True
    temperature = property(get_temperature,set_temperature)
