class OutOfBoundsError(Exception):
    def __init__(self,msg,above):
        Exception.__init__(self,msg)
        self.above = above

class FitVariable(object):
    def __init__(self,name, value = None, limits = None):
        self.name              = name
        self.__base_value      = None
        self._limits           = None
        self.__value           = None
        self.__optimizer_value = None

        if value == None and limits == None:
            raise ValueError("Either a value or range has to be given!")

        if value  != None: self.base_value = value
        if limits != None: self.limits     = limits

        self.value = value

    def check_limits(self):
        if self.limits != None:
            if self.value > self.limits[1]:
                raise OutOfBoundsError("Value %g above allowed region (%g,%g)"%(self.value,self.limits[0],self.limits[1]),above=True)
            if self.value < self.limits[0]:
                raise OutOfBoundsError("Value %g below allowed region (%g,%g)"%(self.value,self.limits[0],self.limits[1]),above=False)

class LinearFitVariable(FitVariable):
    def get_base_value(self):
        return self.__base_value
    def set_base_value(self,v):
        self.__base_value = v
        if self.limits == None:
            self.__shift = 0.0
            self.__scale = v**-1
    base_value = property(get_base_value,set_base_value)

    def get_limits(self):
        return self._limits
    def set_limits(self,l):
        self._limits = l
        if l != None:
            (lower,upper) = l
            self.__scale = (upper-lower)**-1
            self.__shift = -.5*(upper+lower)
        else:
            self.__scale = self.__base_value**-1
            self.__shift = 0.0
    limits = property(get_limits,set_limits)

    def get_value(self):
        return self.__value
    def set_value(self,v):
        self.__value = v
        self.__optimizer_value = self.__scale*(v+self.__shift)
        self.check_limits()
    value = property(get_value,set_value)

    def get_optimizer_value(self):
        return self.__optimizer_value
    def set_optimizer_value(self,v):
        self.__optimizer_value = v
        self.__value = v/self.__scale - self.__shift
        self.check_limits()
    optimizer_value = property(get_optimizer_value,set_optimizer_value)
    

class LogarithmicFitVariable(FitVariable):
    def get_base_value(self):
        return self.__base_value
    def set_base_value(self,v):
        from scipy import log
        self.__base_value = v
        if self.limits == None:
            self.__shift = 0.0
            self.__scale = log(v)**-1
    base_value = property(get_base_value,set_base_value)

    def get_limits(self):
        return self._limits
    def set_limits(self,l):
        from scipy import log
        self._limits =l
        if l != None:
            (lower,upper) = l
            self.__scale = (log(upper)-log(lower))**-1
            self.__shift = -.5*(log(upper)+log(lower))
        else:
            self.__scale = log(self.__base_value)**-1
            self.__shift = 0.0
    limits = property(get_limits,set_limits)

    def get_value(self):
        return self.__value
    def set_value(self,v):
        from math import log
        self.__value = v
        self.__optimizer_value = self.__scale*(log(v)+self.__shift)
        self.check_limits()
    value = property(get_value,set_value)

    def get_optimizer_value(self):
        return self.__optimizer_value
    def set_optimizer_value(self,v):
        from math import exp
        self.__optimizer_value = v
        self.__value = exp(v/self.__scale - self.__shift)
        self.check_limits()
    optimizer_value = property(get_optimizer_value,set_optimizer_value)
