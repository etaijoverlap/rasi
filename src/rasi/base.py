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

class BaseCalculator(object):
    def __init__(self):
        self.__dict__["_input_variables"]  = {}
        self.__dict__["_output_variables"]  = {}
        self.__dict__["_changes"]  = {}

    @staticmethod
    def check_reserved(name):
        if name.startswith("internal_"):
            raise KeyError("Variable names starting with 'internal_' are reserved.")
        if name.startswith("changed_"):
            raise KeyError("Variable names starting with 'changed_' are reserved.")
        if name == "changed":
            raise KeyError("Variable name 'changed' is reserved.")

        

    def set_input_variable(self,name,value):
        BaseCalculator.check_reserved(name)
        self.__dict__["_input_variables"][name] = value
        self.__dict__["_changes"][name] = False
       

    def set_output_variable(self,name,value):
        BaseCalculator.check_reserved(name)
        self.__dict__["_output_variables"][name] = value

    def init_variables(self,inputs=None,outputs = None):
        if inputs == None:
            inputs = {}
        if outputs == None:
            outputs = {}
        self.__dict__["_input_variables"]  = {}
        self.__dict__["_output_variables"]  = {}
        self.__dict__["_changes"]  = {}
        for key,value in inputs.iteritems():
            self.set_input_variable(key,value)
        for key,value in outputs.iteritems():
            self.set_output_variable(key,value)

    def init_input_variables(self,**kwargs):
        if not "_input_variables" in self.__dict__:
            self.__dict__["_input_variables"] = {}
            self.__dict__["_changes"] = {}
        for key,value in kwargs.iteritems():
            self.set_input_variable(key,value)

    def init_output_variables(self,**kwargs):
        if not "_output_variables" in self.__dict__:
            self.__dict__["_output_variables"] = {}
        for key,value in kwargs.iteritems():
            self.set_output_variable(key,value)

    def set_variables(self,kwargs):
        for key,value in kwargs.iteritems():
            if not key in self._input_variables.keys():
                raise AttributeError("%s is not an input variable."%key)
            self._input_variables[key] = value
            

    def __getattr__(self,var):
        if var.startswith("__"):
            raise AttributeError("Attribute %s not found."%var)
        if var.startswith("changed_"):
            varname = var[8:]
            if varname in self._changes.keys():
                return self._changes[varname]
            else:
                raise AttributeError("%s is not an output parameter."%varname)
        if var == "changed":
            return True in self._changes.values()

        if var in self._input_variables.keys():
            return self._input_variables[var]
        if var in self._output_variables.keys():
            return self._output_variables[var]
        raise AttributeError("Parameter %s not found."%var)

    def __setattr__(self,var,val):
        if var.startswith("__"):
            raise AttributeError("Attribute %s not found."%var)
        if var.startswith("internal_"):
            varname = var[9:]
            if varname in self._output_variables.keys():
                self._output_variables[varname]=val
            else:
                raise AttributeError("%s is not an output parameter."%varname)
            return

        if var in self._input_variables.keys():
            self._input_variables[var] = val
            self._changes[var] = True
            return
        if var in self._output_variables.keys():
            raise AttributeError("Direct setting of output variables not permitted. Use the internal_* alias for this purpose.")
        raise AttributeError("Parameter %s not found."%var)

    def __delattr__(self,var):
        raise AttributeError("Deletion of parameters is not permitted (yet).")

    def __dir__(self):
        return self._input_variables.keys() + self._output_variables.keys()

    def update(self):
        changed = self.do_update()
        for key in self._changes:
            self._changes[key]=False
        return changed
