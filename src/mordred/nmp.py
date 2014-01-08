class FullNMPTransition(object):
    def __init__(self):
        self.__lineshape                 = None
        self.__electronic_matrix_element = None
        self.__mlambda                   = None
        self.__changed                   = False
        
    def update(self):
        from scipy.integrate import trapz
        ls_changed = self.lineshape.update()
        em_changed = False #self.electronic_matrix_element.update()
        changed = (self.__changed or ls_changed or em_changed)
        self.__changed = False
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
                
            self.oxidation_rate = oxidation_rate
            self.reduction_rate = reduction_rate
            return True
        return False
        
    
    def set_lineshape(self,ls):
        self.__lineshape = ls
        self.__changed = True
    def get_lineshape(self):
        return self.__lineshape
    lineshape = property(get_lineshape,set_lineshape)

    def set_electronic_matrix_element(self,d):
        self.__electronic_matrix_element = d
        self.__changed = True
    def get_electronic_matrix_element(self):
        return self.__electronic_matrix_element
    electronic_matrix_element = property(get_electronic_matrix_element,set_electronic_matrix_element)

    def set_mlambda(self,mlambda):
        self.__mlambda = mlambda
        self.__changed = True
    def get_mlambda(self):
        return self.__mlambda
    mlambda = property(get_mlambda,set_mlambda)
    
