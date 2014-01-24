"""
    Reliability Analysis of Semiconductor Interfaces -- ElectronicMatrixElement package
    -----------------------------------------------------------------------------------

    In this package, all the calculator objects for electronic matrix elements are defined.
    Electronic matrix elements are necessary to calculate capture and emission rates within
    the non-radiative multi-phonon transition theory... and possibly also in other theories.

"""


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

def lininterpolate(x,x1,y1,x2,y2):
    return (x-x1)/(x2-x1) * y1 + (x2-x)/(x2-x1) * y2

def loginterpolate(x,x1,y1,x2,y2):
    from math import log,exp
    if y1 == 0 or y2 == 0:
        return 0.0
    return exp(lininterpolate(x,x1,log(y1),x2,log(y2)))

def loginterpolatearray(x,x1,y1,x2,y2):
    from numpy import zeros
    result = zeros(y1.shape)
    for i in xrange(len(y1)):
        result[i] = loginterpolate(x,x1,y1[i],x2,y2[i])
    return result


class EMF1DPositionInterpolator(BasicCalculator):
    def __init__(self, **kwargs):
        self.init_input_variables(
                 emf      = None,
                 position = None,
                 )
        self.init_output_variables(
                 oxidation_reservoir = None,
                 reduction_reservoir = None,
                 Ec = None,
                 Ev = None,
                 phi = None
                 )
        self.set_variables(kwargs)


    def do_update(self):
        from numpy import vectorize
        if self.changed:
            x = self.position

            defects = [ defect for defect in self.emf._defects ]
            #defects.sort(key=(lambda defect: float(defect._info)))

            positions = [ float(defect._info) for defect in self.emf._defects ]

            if not (positions[0] < x and x < positions[-1]):
                raise ValueError("Defect position %g outside of simulated interval (%g,%g)"%(x,positions[0],positions[-1]))

            i_lower  = None
            i_higher = None
            for i,position in enumerate(positions):
                if x < position:
                    i_higher = i
                    i_lower  = i-1
                    break

            if defects[i_lower].oxidation_reservoir.keys() != defects[i_higher].oxidation_reservoir.keys():
                raise ValueError("Oxidation reservoirs of EMF defect objects contain different reservoirs")
            if defects[i_lower].reduction_reservoir.keys() != defects[i_higher].reduction_reservoir.keys():
                raise ValueError("Reduction reservoirs of EMF defect objects contain different reservoirs")

            oxidation_reservoir = {}
            for name in defects[i_lower].oxidation_reservoir.iterkeys():
                E_lower ,d_lower   = defects[i_lower].oxidation_reservoir[name]
                E_higher,d_higher  = defects[i_higher].oxidation_reservoir[name]
                x_lower  = positions[i_lower]
                x_higher = positions[i_higher] 
                if (E_lower != E_higher).any():
                    raise ValueError("Energy grids of EMF defect objects don't match")
                oxidation_reservoir[name] = E_lower,lininterpolate(x,x_lower,d_lower,x_higher,d_higher)

            reduction_reservoir = {}
            for name in defects[i_lower].reduction_reservoir.iterkeys():
                E_lower ,d_lower   = defects[i_lower].reduction_reservoir[name]
                E_higher,d_higher  = defects[i_higher].reduction_reservoir[name]
                x_lower  = positions[i_lower]
                x_higher = positions[i_higher] 
                if (E_lower != E_higher).any():
                    raise ValueError("Energy grids of EMF defect objects don't match")
                reduction_reservoir[name] = E_lower,lininterpolate(x,x_lower,d_lower,x_higher,d_higher)

            self.internal_oxidation_reservoir = oxidation_reservoir
            self.internal_reduction_reservoir = reduction_reservoir

            self.internal_Ec  = lininterpolate(x,x_lower,defects[i_lower].Ec ,x_higher,defects[i_higher].Ec )
            self.internal_Ev  = lininterpolate(x,x_lower,defects[i_lower].Ev ,x_higher,defects[i_higher].Ev )
            self.internal_phi = lininterpolate(x,x_lower,defects[i_lower].phi,x_higher,defects[i_higher].phi)
            return True
        return False

