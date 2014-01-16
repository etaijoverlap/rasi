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



class ClassicalLineShape(object):
    def __init__(self, k_occupied=None, k_unoccupied=None, equilibrium_shift=None, thermodynamic_level=None, temperature=None):
        """
            Properties:
            k_occupied    ... Spring constant for the mode when the electron state 
                        is occupied
            k_unoccupied  ... Spring constant for the model when the electron state
                        is unoccupied
            equilibrium_shift   ... Configuration shift along the modal coordinate induced
                        by the transition
            thermodynamic_level     ... `Thermodynamic level` (i.e. energetic shift between the 
                         parabolas referenced to some defined energy) for the transition
            temperature     ... Temperature (optional), used as default in later calculations
        """
        self.__k_occupied          = None
        self.__k_unoccupied        = None
        self.__equilibrium_shift   = None
        self.__thermodynamic_level = None
        self.__temperature         = None
        self.__changed             = False

        if k_occupied          != None: self.k_occupied          = k_occupied
        if k_unoccupied        != None: self.k_unoccupied        = k_unoccupied
        if equilibrium_shift   != None: self.equilibrium_shift   = equilibrium_shift
        if thermodynamic_level != None: self.thermodynamic_level = thermodynamic_level
        if temperature         != None: self.temperature         = temperature

    @staticmethod
    def __partitionfunction(Momega2,T):
        from scipy.constants import pi
        from scipy.constants import k as kB
        from scipy import sqrt
        # (4.79) in my thesis
        return sqrt(2*pi*kB*T/Momega2)

    def update(self):
        from numpy import vectorize
        if self.__changed:
            self.oxidation = vectorize(self._oxidation)
            self.reduction = vectorize(self._reduction)
            self.__changed = False
            return True
        return False

    def crossings(self,E):
        from scipy import sqrt
        """ Calculation of the intersection coordinates for the
            parabolas. (4.82) in my thesis. 
        """
        # Typo alert equation (4.82): Q' in the square root should be Q'**2
        k0 = self.k_occupied; kp = self.k_unoccupied ; shift = self.equilibrium_shift
        lvl = self.thermodynamic_level

        # The expressions are of the form ( a +/- sqrt(b) ) / D
        a = kp*shift
        b = k0*kp*shift**2 +(k0-kp)*(E-lvl)
        D = kp-k0

        if b < 0:
            raise ValueError("Parabolas have no real-valued crossings")
        return ( ( a + sqrt(b) ) / D, (a - sqrt(b) ) / D )

    def get_k_occupied(self):
        return self.__k_occupied
    def set_k_occupied(self,k):
        self.__k_occupied = k
        self.__changed = True
    k_occupied = property(get_k_occupied,set_k_occupied)

    def get_k_unoccupied(self):
        return self.__k_unoccupied 
    def set_k_unoccupied(self,k):
        self.__k_unoccupied  = k
        self.__changed = True
    k_unoccupied = property(get_k_unoccupied,set_k_unoccupied)

    def get_equilibrium_shift(self):
        return self.__equilibrium_shift
    def set_equilibrium_shift(self,s):
        self.__equilibrium_shift = s
        self.__changed = True
    equilibrium_shift = property(get_equilibrium_shift,set_equilibrium_shift)

    def get_thermodynamic_level(self):
        return self.__thermodynamic_level
    def set_thermodynamic_level(self,t):
        self.__thermodynamic_level = t
        self.__changed = True
    thermodynamic_level = property(get_thermodynamic_level,set_thermodynamic_level)

    def get_temperature(self):
        return self.__temperature
    def set_temperature(self,T):
        self.__temperature = T
        self.__changed     = True
    temperature = property(get_temperature,set_temperature)


        
    # To all Java/C++/... programmers: remember that those methods are DETACHABLE!
    def _oxidation(self,E):
        from scipy import exp,sqrt,pi
        from scipy.constants import k as kB
        T = self.temperature
        k0 = self.k_occupied     ; kp = self.k_unoccupied
        shift = self.equilibrium_shift ; lvl = self.thermodynamic_level
        Momega20 = 2*k0    ; Momega2P = 2*kp
        Z0 = self.__partitionfunction(Momega20,T)

        try:
            Q1,Q2 = self.crossings(E)
        except ValueError:
            return 0.
        denom1 = abs(Momega20*Q1+Momega2P*(shift-Q1))
        denom2 = abs(Momega20*Q2+Momega2P*(shift-Q2))

        # (4.83) Beware of another typo!
        return Z0**-1*(exp(-k0*Q1**2/(kB*T))/denom1+exp(-k0*Q2**2/(kB*T))/denom2)
    
    def _reduction(self,E):
        from scipy import exp,sqrt,pi
        from scipy.constants import k as kB
        T = self.temperature
        k0 = self.k_occupied     ; kp = self.k_unoccupied
        shift = self.equilibrium_shift ; lvl = self.thermodynamic_level
        Momega20 = 2*k0    ; Momega2P = 2*kp
        # First the partition function (4.80) in my thesis
        ZP = self.__partitionfunction(Momega2P,T)

        try:
            Q1,Q2 = self.crossings(E)
        except ValueError:
            return 0.

        denom1 = abs(Momega20*Q1+Momega2P*(shift-Q1))
        denom2 = abs(Momega20*Q2+Momega2P*(shift-Q2))

        # (4.84) Beware of another typo!
        return ZP**-1*(exp(-kp*(Q1-shift)**2/(kB*T))/denom1+exp(-kp*(Q2-shift)**2/(kB*T))/denom2)

