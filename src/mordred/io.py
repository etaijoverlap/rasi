class EMFDefect:
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
        

class EMF:
    def __init__(self,filename=None):
        self._filename = None
        self._version  = None
        self._dimensionality = None
        self._devicesimulator = None
        self._theorylevel = None
        self._type = None
        self._defects = []
        if filename != None:
            self.read(filename)
        
    def read(self,emffile):
        if type(emffile) == str:
            emffile = open(emffile)
        self._filename = emffile.name
        title = emffile.next()
        if not title.startswith("EMFILE"):
            raise ValueError("EMF file %s: invalid format"%emffile.name)
        self._version = int(title.split()[1][1:])
        [self._dimensionality,self._devicesimulator,self._theorylevel,self._type] = emffile.next().split()
        ndefects = int(emffile.next())
        for i in range(ndefects):
            defect = EMFDefect()
            defect.read_from_stream(emffile)
            self._defects.append(defect)
