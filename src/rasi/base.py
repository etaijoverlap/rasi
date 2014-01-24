"""
    Reliability Analysis of Semiconductor Interfaces -- BASE package
    ----------------------------------------------------------------

    This package defines basic objects and functions that are usually invisible
    to the user. If you just plan to use this library rather than extending
    it by writing new objects or methods, this is most likely the worst place
    to start. Go away, there is nothing here for you to find.
    If you are actually planning to add some functionality to RASI, be welcome
    and make yourself a home.
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

class BasicCalculator(object):
    """
        BasicCalculator:
        ---------------

        The calculation of experimentally observable quantities in this library
        usually proceeds within a hierachy of objects, where each object receives
        input from the user and the child objects and provides output to either
        the user or other objects which are in the next level of the hierachy.

        A basic concept of this computational tree is that the output parameters
        of each node are only valid after an "update()" signal has been processed.
        On receiving this signal, the node updates its internal state and at the
        same time propagates the signal to its child nodes.

        The BasicCalculator base class implements some functionality to make the
        development of node objects easier. A node object should always be derived
        from BasicCalculator, and have a constructor of the form:

        def __init__(self,**kwargs):
            self.init_input_variables(
                                      <input variable 1> = <initial value 1>,
                                      <input variable 2> = <initial value 2>,
                                                 ...
                                     )
            self.init_output_variables(
                                      <output variable 1> = <initial value 3>,
                                      <output variable 2> = <initial value 4>,
                                                 ...
                                      )
            self.set_variables(kwargs)

        Further, objects which are derived from BasicCalculator should NOT implement
        the update method directly (otherwise bad things will happen), but instead
        implement a 

            def do_update(self):
                 <update code here>

        method, which is called by the update method of BasicCalculator.

        After input and output variables are registered as in the example code above,
        they can in principle be used like normal variables of the object. The
        BasicCalculator class keeps track of changes to the input variables, which
        can be checked via the *_changed parameters. These flags are automatically
        reset after the update method has been processed.

        If for example "foo" is an input parameter of the node object "Bar", the following code holds:
            
        bar = Bar()

        print bar.foo_changed  # prints "False"

        bar.foo = 4

        print bar.foo_changed # prints "True"

        bar.update()

        print bar.foo_changed # prints "False"

        There is also a global flag "bar.changed" which signals if ANY of the input variables
        has changed.

        The output parameters can be directly read but not directly written. The idea is 
        to avoid accidential overwriting of output data by the user. Writing to output
        parameters has to proceed via the internal_* aliases. 

        If for example "spam" is an output parameter of the node object "Eggs", the following code holds:

        eggs = Eggs()

        eggs.spam = 4 # throws AttributeError

        eggs.internal_spam = 4 # Works, but should only be used in the update method

        print eggs.spam # prints "4"
    """
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
        BasicCalculator.check_reserved(name)
        self.__dict__["_input_variables"][name] = value
        self.__dict__["_changes"][name] = False
       

    def set_output_variable(self,name,value):
        BasicCalculator.check_reserved(name)
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
            self.__setattr__(key,value)
            

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
