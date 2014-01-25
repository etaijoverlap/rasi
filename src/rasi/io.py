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


class EMFv1Defect:
    def __init__(self,iterable=None):
        self._info = None
        self.Ec = None
        self.Ev = None
        self.phi = None
        self.oxidation_reservoir = {}
        self.reduction_reservoir = {}
        if iterable != None:
            self.read_from_stream(iterable)

    def __read_reservoir_set(self,emffile):
        from scipy import array
        setname = emffile.next().strip()
        nlines  = int(emffile.next())
        Ed = []
        for i in xrange(nlines):
            [Eval,dval] = [ float(val) for val in emffile.next().split() ]
            Ed.append((Eval,dval))
        Ed.sort(key=(lambda (Eval,dval):Eval))
        E = []; d = []
        for Eval,dval in Ed:
            E.append(Eval); d.append(dval)
        return setname,array(E),array(d)
        
    def read_from_stream(self,emffile):
        self._info = emffile.next()
        [self.phi,self.Ec,self.Ev] = [float(val) for val in emffile.next().split() ]
        [noxi,nred] = [ int(val) for val in emffile.next().split() ]
        for i in range(noxi): 
            setname,E,d = self.__read_reservoir_set(emffile)
            self.oxidation_reservoir[setname] = (E,d)
        for i in range(nred): 
            setname,E,d = self.__read_reservoir_set(emffile)
            self.reduction_reservoir[setname] = (E,d)

class EMFv1:
    def __init__(self):
        self._filename = None
        self._dimensionality = None
        self._devicesimulator = None
        self._theorylevel = None
        self._type = None
        self._defects = []
        
    def read(self,emffile):
        [self._dimensionality,self._devicesimulator,self._theorylevel,self._type] = emffile.next().split()
        ndefects = int(emffile.next())
        for i in range(ndefects):
            defect = EMFv1Defect()
            defect.read_from_stream(emffile)
            self._defects.append(defect)

def EMF(emffile):
    if type(emffile) == str:
        emffile = open(emffile)
    title = emffile.next()
    if not title.startswith("EMFILE"):
        raise ValueError("EMF file %s: invalid format"%emffile.name)
    if int(title.split()[1][1:]) == 1:
        emf = EMFv1()
        emf._filename = emffile.name
        emf.read(emffile)
        return emf
    else:
        raise ValueError("EMF version not recognized.")
