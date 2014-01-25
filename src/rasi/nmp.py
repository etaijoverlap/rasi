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

class FullNMPTransition(BasicCalculator):
    def __init__(self, **kwargs):
        self.init_input_variables(
                    lineshape = None,
                    electronic_matrix_element = None,
                    mlambda                   = None
                    )
        self.init_output_variables(
                oxidation_rate = None,
                reduction_rate = None
                )
        self.set_variables(kwargs)
        

    def do_update(self):
        from scipy.integrate import trapz
        ls_changed = self.lineshape.update()
        em_changed = self.electronic_matrix_element.update()
        changed = (self.changed or ls_changed or em_changed)
        if changed:
            eme = self.electronic_matrix_element
            oxidation_lsf = self.lineshape.oxidation
            reduction_lsf = self.lineshape.reduction
            mlambda       = self.mlambda
            oxidation_rate = 0.
            reduction_rate = 0.
            Ev = eme.Ev
            for E,d in eme.oxidation_reservoir.itervalues():
                oxidation_rate += mlambda*trapz(d*oxidation_lsf(E-Ev),E)
            for E,d in eme.reduction_reservoir.itervalues():
                reduction_rate += mlambda*trapz(d*reduction_lsf(E-Ev),E)
                
            self.internal_oxidation_rate = oxidation_rate
            self.internal_reduction_rate = reduction_rate
            return True
        return False
    
    
class ColdCarrierNumerical(object):
    def __init__(self):
        self.__CBE_interface              = None
        self.__VBE_interface              = None
        self.__VBE_defect                 = None
        self.__thermodynamic_level        = None
        self.__temperature                = None
        self.__electron_density_interface = None
        self.__hole_density_interface     = None
        self.__tunneling_factor           = None
        self.__lineshape                  = None

        self.__changed = False

    def do_update(self):
        from scipy.constants import k as kB
        from math import exp

        changed = self.__changed
        self.__changed = False
        ls_changed = self.lineshape.update()
        
        if changed or ls_changed:
            Eci = self.CBE_interface
            Evi = self.VBE_interface
            Evd = self.VBE_defect
            ET  = self.thermodynamic_level
            T   = self.temperature
            n   = self.electron_density_interface
            p   = self.hole_density_interface
            TF  = self.tunneling_factor

            lsf = self.lineshape

            c_p = p*TF*lsf.oxidation(Evi-Evd)
            c_n = n*TF*lsf.reduction(Eci-Evd)
            #Detailed balance
            self.oxidation_rate = c_p + c_n * exp(-(Eci-ET)/(kB*T))
            self.reduction_rate = c_n + c_p * exp(-(ET-Evi)/(kB*T))
            return True
        return False

    def get_CBE_interface(self):
        return self.__CBE_interface
    def set_CBE_interface(self,v):
        self.__CBE_interface = v
        self.__changed = True
    CBE_interface = property(get_CBE_interface,set_CBE_interface)

    def get_VBE_interface(self):
        return self.__VBE_interface
    def set_VBE_interface(self,v):
        self.__VBE_interface = v
        self.__changed = True
    VBE_interface = property(get_VBE_interface,set_VBE_interface)

    def get_VBE_defect(self):
        return self.__VBE_defect
    def set_VBE_defect(self,v):
        self.__VBE_defect = v
        self.__changed = True
    VBE_defect = property(get_VBE_defect,set_VBE_defect)

    def get_thermodynamic_level(self):
        return self.__thermodynamic_level
    def set_thermodynamic_level(self,v):
        self.__thermodynamic_level = v
        self.__changed = True
    thermodynamic_level = property(get_thermodynamic_level,set_thermodynamic_level)

    def get_temperature(self):
        return self.__temperature
    def set_temperature(self,v):
        self.__temperature = v
        self.__changed = True
    temperature = property(get_temperature,set_temperature)

    def get_electron_density_interface(self):
        return self.__electron_density_interface
    def set_electron_density_interface(self,v):
        self.__electron_density_interface = v
        self.__changed = True
    electron_density_interface = property(get_electron_density_interface,set_electron_density_interface)

    def get_hole_density_interface(self):
        return self.__hole_density_interface
    def set_hole_density_interface(self,v):
        self.__hole_density_interface = v
        self.__changed = True
    hole_density_interface = property(get_hole_density_interface,set_hole_density_interface)

    def get_tunneling_factor(self):
        return self.__tunneling_factor
    def set_tunneling_factor(self,v):
        self.__tunneling_factor = v
        self.__changed = True
    tunneling_factor = property(get_tunneling_factor,set_tunneling_factor)

    def get_lineshape(self):
        return self.__lineshape
    def set_lineshape(self,v):
        self.__lineshape = v
        print self.__lineshape
        self.__changed = True
    lineshape = property(get_lineshape,set_lineshape)



#    def get_XXX(self):
#        return self.__XXX
#    def set_XXX(self,v):
#        self.__XXX = v
#        self.__changed = True
#    XXX = property(get_XXX,set_XXX)

