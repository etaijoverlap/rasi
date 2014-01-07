
def calc_modal_vector(atoms1,atoms2):
    """
        Calculate the 'modal vector', i.e. the difference vector between the two configurations.
        The minimum image convention is applied!
    """
    from scipy.linalg import inv
    from scipy        import array,dot
    from scipy        import sign,floor
    cell1 = atoms1.get_cell()
    cell2 = atoms2.get_cell()

    # The cells need to be the same (otherwise the whole process won't make sense)
    if (cell1 != cell2).any():
        raise ValueError("Encountered different cells in atoms1 and atoms2. Those need to be the same.")
    cell = cell1

    icell = inv(cell)
                                            
    frac1 = atoms1.get_scaled_positions()
    frac2 = atoms2.get_scaled_positions()
    modal_vector_frac = frac1 - frac2
    for i in range(modal_vector_frac.shape[0]):
        for j in range(modal_vector_frac.shape[1]):
            if abs(modal_vector_frac[i,j]) > .5:
                value = modal_vector_frac[i,j]
                vsign = sign(modal_vector_frac[i,j])
                absvalue = abs(value)
                modal_vector_frac[i,j] = value - vsign*floor(absvalue+.5)
    return dot(modal_vector_frac,cell)


class ClassicalLineShape(object):
    def numpy_aware_ls(fn):
        import numpy
        def numpy_aware_version(self,E,T=None):
            if isinstance(E,numpy.ndarray):
                result = numpy.ndarray(E.shape)
                for i in xrange(E.shape[0]):
                    result[i] = numpy_aware_version(self,E[i],T)
            else:
                result = fn(self,E,T)
            return result
        return numpy_aware_version

    def __init__(self, kocc,kunocc,shift,lvl,T=None):
        """
            kocc    ... Spring constant for the mode when the electron state 
                        is occupied
            kunocc  ... Spring constant for the model when the electron state
                        is unoccupied
            shift   ... Configuration shift along the modal coordinate induced
                        by the transition
            lvl     ... `Thermodynamic level` (i.e. energetic shift between the 
                         parabolas referenced to some defined energy) for the transition
              T     ... Temperature (optional), used as default in later calculations
        """
        self.kocc   = kocc
        self.kunocc = kunocc
        self.shift  = shift
        self.lvl    = lvl
        self.T      = T
    def __check_T(self,T):
        if T == None:
            if self.T != None:
                T = self.T
            else:
                raise ValueError("No temperature given, no default temperature set for the lineshape.")
        return T

    def crossings(self,E):
        from scipy import sqrt
        """ Calculation of the intersection coordinates for the
            parabolas. (4.82) in my thesis. 
        """
        # Typo alert equation (4.82): Q' in the square root should be Q'**2
        k0 = self.kocc; kp = self.kunocc ; shift = self.shift
        lvl = self.lvl

        # The expressions are of the form ( a +/- sqrt(b) ) / D
        a = kp*shift
        b = k0*kp*shift**2 +(k0-kp)*(E-lvl)
        D = kp-k0

        if b < 0:
            raise ValueError("Parabolas have no real-valued crossings")
        return ( ( a + sqrt(b) ) / D, (a - sqrt(b) ) / D )

    def __partitionfunction(self,Momega2,T):
        from scipy.constants import pi
        from scipy.constants import k as kB
        from scipy import sqrt
        # (4.79) in my thesis
        return sqrt(2*pi*kB*T/Momega2)

        
    @numpy_aware_ls
    def oxidation(self,E,T=None):
        from scipy import exp,sqrt,pi
        from scipy.constants import k as kB
        T = self.__check_T(T)
        k0 = self.kocc     ; kp = self.kunocc
        shift = self.shift ; lvl = self.lvl
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
    
    @numpy_aware_ls
    def reduction(self,E,T=None):
        from scipy import exp,sqrt,pi
        from scipy.constants import k as kB
        T = self.__check_T(T)
        k0 = self.kocc     ; kp = self.kunocc
        shift = self.shift ; lvl = self.lvl
        Momega20 = 2*k0    ; Momega2P = 2*kp
        # First the partition function (4.80) in my thesis
        ZP = sqrt(2*pi*kB*T/Momega2P)

        try:
            Q1,Q2 = self.crossings(E)
        except ValueError:
            return 0.

        denom1 = abs(Momega20*Q1+Momega2P*(shift-Q1))
        denom2 = abs(Momega20*Q2+Momega2P*(shift-Q2))

        # (4.84) Beware of another typo!
        return ZP**-1*(exp(-kp*(Q1-shift)**2/(kB*T))/denom1+exp(-kp*(Q2-shift)**2/(kB*T))/denom2)


