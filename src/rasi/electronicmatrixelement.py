def lininterpolate(x,x1,y1,x2,y2):
    return (x-x1)/(x2-x1) * y1 + (x2-x)/(x2-x1) * y2

def loginterpolate(x,x1,y1,x2,y2):
    return exp(lininterpolate(x,x1,log(y1),x2,log(y2)))


class EMF1DPositionInterpolator(object):
    def __init__(self, emf=None, position=None):
        self.__emf      = None
        self.__position = None
        self.__changed  = False

        if emf      != None: self.emf      = emf
        if position != None: self.position = position

    def update(self):
        if self.__changed:
            x = self.position

            defects = [ defect for defect in self.emf._defects ]
            defects.sort(key=(lambda defect: float(defect._info)))

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
                if E_lower != E_higher:
                    raise ValueError("Energy grids of EMF defect objects don't match")
                oxidation_reservoir[name] = E_lower,loginterpolate(x,x_lower,d_lower,x_higher,d_higher)

            reduction_reservoir = {}
            for name in defects[i_lower].reduction_reservoir.iterkeys():
                E_lower ,d_lower   = defects[i_lower].reduction_reservoir[name]
                E_higher,d_higher  = defects[i_higher].reduction_reservoir[name]
                x_lower  = positions[i_lower]
                x_higher = positions[i_higher] 
                if E_lower != E_higher:
                    raise ValueError("Energy grids of EMF defect objects don't match")
                reduction_reservoir[name] = E_lower,loginterpolate(x,x_lower,d_lower,x_higher,d_higher)

            self.oxidation_reservoir = oxidation_reservoir
            self.reduction_reservoir = reduction_reservoir
            self.__changed = False
            return True
        return False

    def get_emf(self):
        return self.__emf
    def set_emf(self,emf):
        self.__emf = emf
        self.__changed = True
    emf = property(get_emf,set_emf)

    def get_position(self):
        return self.__position
    def set_position(self,d):
        self.__position = d
        self.__changed = True
    position = property(get_position,set_position)
