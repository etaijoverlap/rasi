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

class ClassicalLineShape(BasicCalculator):
    def __init__(self, **kwargs):
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

        self.init_input_variables(
                   k_occupied          = None,
                   k_unoccupied        = None,
                   equilibrium_shift   = None,
                   thermodynamic_level = None,
                   temperature         = None
                   )
        self.init_output_variables(
                   oxidation = None,
                   reduction = None
                   )           

    @staticmethod
    def __partitionfunction(Momega2,T):
        from scipy.constants import pi
        from scipy.constants import k as kB
        from scipy import sqrt
        # (4.79) in my thesis
        return sqrt(2*pi*kB*T/Momega2)

    def do_update(self):
        from numpy import vectorize
        if self.changed:
            self.internal_oxidation = vectorize(self._oxidation)
            self.internal_reduction = vectorize(self._reduction)
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

