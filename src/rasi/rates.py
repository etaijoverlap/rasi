class RateCalculator4State(object):
    def __init__(self, recharge_primary = None, recharge_secondary = None, nu = None, energies = None, temperature = None):
        self.__recharge_primary = None
        self.__recharge_secondary = None
        self.__nu = 0.
        self.__energies = {}
        self.__temperature = None
        self.__primary_changed = False
        self.__secondary_changed = False
        self.__changed = True
        self.rates = {}
        self.timeconstants = {}

        if recharge_primary   != None: self.recharge_primary   = recharge_primary
        if recharge_secondary != None: self.recharge_secondary = recharge_secondary
        if nu                 != None: self.nu                 = nu
        if energies           != None: self.energies           = energies
        if temperature        != None: self.temperature        = temperature

    def update(self):
        from scipy.constants import k as kB
        from scipy import exp
        changed = False
        rates = self.rates
        if self.recharge_primary.update() or self.__primary_changed:
            rates["1->2'"] = self.recharge_primary.oxidation_rate
            rates["2'->1"] = self.recharge_primary.reduction_rate
            self.__primary_changed = False
            changed = True
        if self.recharge_secondary.update() or self.__secondary_changed:
            rates["1'->2"] = self.recharge_secondary.oxidation_rate
            rates["2->1'"] = self.recharge_secondary.reduction_rate
            self.__secondary_changed = False
            changed = True
        if self.__changed:
            E = self.energies ; nu = self.nu ; T = self.temperature
            for transition,Eb in E.iteritems():
                rates[transition] = nu * exp( -Eb/(kB*T) )
            self.__changed = False
            changed = True

        if changed:
            for key,rate in rates.iteritems():
                self.timeconstants[key] = rate**-1
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
        return self.__temperature
    def set_temperature(self,temperature):
        self.__temperature = temperature
        self.__changed = True
    temperature = property(get_temperature,set_temperature)

    def get_recharge_primary(self):
        return self.__recharge_primary
    def set_recharge_primary(self,r):
        self.__recharge_primary = r
        self.__primary_changed = True
    recharge_primary = property(get_recharge_primary,set_recharge_primary)

    def get_recharge_secondary(self):
        return self.__recharge_secondary
    def set_recharge_secondary(self,s):
        self.__recharge_secondary = s
        self.__secondary_changed = True
    recharge_secondary = property(get_recharge_secondary,set_recharge_secondary)
