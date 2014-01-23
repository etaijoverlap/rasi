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

from rasi.base import BasicCalculator

class SmearedLineShape(BasicCalculator):
    def __init__(self,**kwargs):
        self.init_input_variables(
                                   discrete_lineshape = None,
                                   smearing           = None,
                                 )
        self.init_output_variables()
        self.set_variables(kwargs)

    def do_update(self):
        ls_changed =  self.discrete_lineshape.update()

        return self.changed or ls_changed

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

class DiscreteLineShape(BasicCalculator):
    def __init__(self, **kwargs):
        self.init_input_variables(
                                   overlaps            = None,
                                   temperature         = None,
                                   thermodynamic_level = None
                                 )
        self.init_output_variables(
                                   oxidation_energies = None,
                                   oxidation_weights  = None,
                                   reduction_energies = None,
                                   reduction_weights  = None
                                  )
        self.set_variables(kwargs)

    def update(self):
        from scipy.constants import k as kB
        from scipy import exp,array
        overlaps_changed = self.overlaps.update()

        if self.changed or overlaps_changed:
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
                    reduction.append((ET+E_occupied[i]-E_unoccupied[j] , p_unoccupied[j] * S[i,j]**2))
            oxidation.sort(key=(lambda (E,w): E))
            reduction.sort(key=(lambda (E,w): E))

            self.internal_oxidation_weights  = array([ w for E,w in oxidation ])
            self.internal_oxidation_energies = array([ E for E,w in oxidation ])

            self.internal_reduction_weights  = array([ w for E,w in reduction ])
            self.internal_reduction_energies = array([ E for E,w in reduction ])
            return True
        return False

