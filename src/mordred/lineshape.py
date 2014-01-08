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
    def __init__(self):
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

class LineShapeBuffer(object):
    def __init__(self):
        self.__lineshape   = None
        self.__energy_grid = None
        self.__changed = False

    def update(self):
        from scipy.interpolate import interp1d
        ls_changed = self.lineshape.update()

        changed = self.__changed or ls_changed
        self.__changed = False

        if changed:
            self.oxidation_values = self.lineshape.oxidation(self.energy_grid)
            self.reduction_values = self.lineshape.reduction(self.energy_grid)

            self.oxidation = interp1d(self.energy_grid,self.oxidation_values,bounds_error=False,fill_value=0.)  
            self.reduction = interp1d(self.energy_grid,self.reduction_values,bounds_error=False,fill_value=0.)  
            return True
        return False

    def get_lineshape(self):
        return self.__lineshape
    def set_lineshape(self,l):
        self.__lineshape = l
        self.__changed = True
    lineshape = property(get_lineshape,set_lineshape)

    def get_energy_grid(self):
        return self.__energy_grid
    def set_energy_grid(self,E):
        self.__energy_grid = E
    energy_grid = property(get_energy_grid,set_energy_grid)

class LineShapeEnergyCorrector(object):
    def __init__(self):
        self.__lineshape         = None
        self.__correction_energy = 0.
        self.__changed = False

    def update(self):
        return self.lineshape.update()

    def oxidation(self,E):
        E_c = self.correction_energy
        return self.lineshape.oxidation(E + E_c)

    def reduction(self,E):
        E_c = self.correction_energy
        return self.lineshape.reduction(E + E_c)

    def get_lineshape(self):
        return self.__lineshape
    def set_lineshape(self,l):
        self.__lineshape = l
        self.__changed = True
    lineshape = property(get_lineshape,set_lineshape)

    def set_correction_energy(self,E):
        self.__correction_energy = E
        self.__changed = True
    def get_correction_energy(self):
        return self.__correction_energy
    correction_energy = property(get_correction_energy,set_correction_energy)

