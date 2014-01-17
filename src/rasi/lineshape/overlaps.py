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

def gaussian(sigma,x0,x):
    from math  import pi
    from scipy import exp
    return ((2.*pi*sigma**2)**-.5)*exp(-(x-x0)**2/(2*sigma**2)) 

class SmearedLineShape(object):
    def __init__(self,discrete_lineshape=None, smearing=None, accuracy=None):
        self.__discrete_lineshape = None
        self.__smearing           = None
        self.__accuracy           = 1e-12 # Has no effect yet!!

        self.__changed = False

        if discrete_lineshape != None: self.discrete_lineshape = discrete_lineshape
        if smearing           != None: self.smearing           = smearing
        if accuracy           != None: self.accuracy           = accuracy

    def update(self):
        ls_changed =  self.discrete_lineshape.update()
        self_changed = self.__changed
        self.__changed = False

        return self_changed or ls_changed

    def __smeared(self,energies,weights,E):
        sigma = self.smearing
        value = 0.
        for E0,weight in zip(energies,weights):
            value += weight*gaussian(sigma,E0,E)
        return value

    def oxidation(self,E):
        energies = self.discrete_lineshape.oxidation_energies
        weights  = self.discrete_lineshape.oxidation_weights
        return self.__smeared(energies,weights,E)

    def reduction(self,E):
        energies = self.discrete_lineshape.reduction_energies
        weights  = self.discrete_lineshape.reduction_weights
        return self.__smeared(energies,weights,E)

class DiscreteLineShape(object):
    def __init__(self, overlaps=None, thermodynamic_level=None, temperature = None):
        self.__overlaps            = None
        self.__temperature         = None
        self.__thermodynamic_level = None

        self.__changed = False

        if overlaps            != None: self.overlaps            = overlaps
        if temperature         != None: self.temperature         = temperature
        if thermodynamic_level != None: self.thermodynamic_level = thermodynamic_level

    def update(self):
        from scipy.constants import k as kB
        from scipy import exp,array
        self_changed     = self.__changed
        overlaps_changed = self.overlaps.update()
        self.__changed = False

        if self_changed or overlaps_changed:
            overlaps = self.overlaps
            T        = self.temperature
            S        = overlaps.overlap_matrix
            ET       = self.thermodynamic_level

            E_occupied   = overlaps.energies_occupied
            E_unoccupied = overlaps.energies_unoccupied

            p_occupied   = exp(-E_occupied/(kB*T))
            p_unoccupied = exp(-E_unoccupied/(kB*T))

            if hasattr(overlaps,"partition_function_occupied"):
                Z_occupied = overlaps.partition_function_occupied(T)
            else:
                Z_occupied = p_occupied.sum()

            if hasattr(overlaps,"partition_function_unoccupied"):
                Z_unoccupied = overlaps.partition_function_unoccupied(T)
            else:
                Z_unoccupied = p_unoccupied.sum()

            p_occupied   /= Z_occupied
            p_unoccupied /= Z_unoccupied

            oxidation = []
            reduction = []
            for i in xrange(len(p_occupied)):
                for j in xrange(len(p_unoccupied)):
                    oxidation.append((ET+E_occupied[i]-E_unoccupied[j] , p_occupied[i]   * S[i,j]**2))
                    reduction.append((ET+E_occupied[i]-E_unoccupied[j] , p_unoccupied[j] * S[j,i]**2))
            oxidation.sort(key=(lambda (E,w): E))
            reduction.sort(key=(lambda (E,w): E))

            self.oxidation_weights  = array([ w for E,w in oxidation ])
            self.oxidation_energies = array([ E for E,w in oxidation ])

            self.reduction_weights  = array([ w for E,w in reduction ])
            self.reduction_energies = array([ E for E,w in reduction ])
            return True
        return False

    def get_overlaps(self):
        return self.__overlaps
    def set_overlaps(self,o):
        self.__overlaps = o
        self.__changed = True
    overlaps = property(get_overlaps,set_overlaps)

    def get_temperature(self):
        return self.__temperature
    def set_temperature(self,T):
        self.__temperature = T
        self.__changed = True
    temperature = property(get_temperature,set_temperature)

    def get_thermodynamic_level(self):
        return self.__thermodynamic_level
    def set_thermodynamic_level(self,ET):
        self.__thermodynamic_level = ET
        self.__changed = True
    thermodynamic_level = property(get_thermodynamic_level,set_thermodynamic_level)
