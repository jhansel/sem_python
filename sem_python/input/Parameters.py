from ..utilities.conversion_utilities import stringToBool, stringToInt, stringToFloat, stringToFloatList, stringToFunction, stringToStringList
from ..utilities.error_utilities import error


## Class for declaring and retrieving input parameters
class Parameters(object):

    def __init__(self, factory):
        self.factory = factory

        self.descriptions = dict()
        self.values = dict()
        self.types = dict()
        self.string_selections = dict()

    ## Retrieves an input parameter
    # @param name  the name of the input parameter to retrieve
    # @return the value of the parameter
    def get(self, name):
        # first check if parameter was manually set for this object
        if name in self.values:
            return self.values[name]
        else:
            # check if parameter was registered
            if name in self.descriptions:
                # check if parameter exists in factory
                if self.factory.hasParameter(name):
                    return self.factory.getParameter(name)
                else:
                    error("The parameter '" + name + "' was not set and does not have a default.")
            else:
                error("'" + name + "' is not a registered parameter.")

    ## Returns whether or not a parameter is available
    # @param name  the name of the parameter for which to check existence
    # @return whether or not the parameter is available
    def has(self, name):
        if name in self.values:
            return True
        else:
            return False

    ## Returns whether or not a parameter is registered
    # @param name  the name of the parameter for which to check registration
    # @return whether or not the parameter is registered
    def hasRegisteredParam(self, name):
        if name in self.descriptions:
            return True
        else:
            return False

    ## Sets an input parameter
    # @param name  name of the parameter to set
    # @param value  value of the parameter as a string
    def set(self, name, value):
        if name in self.descriptions:
            if self.types[name] == "string_list":
                self.values[name] = stringToStringList(value)
            elif self.types[name] == "bool":
                self.values[name] = stringToBool(value)
            elif self.types[name] == "int":
                self.values[name] = stringToInt(value)
            elif self.types[name] == "float":
                self.values[name] = stringToFloat(value)
            elif self.types[name] == "float_list":
                self.values[name] = stringToFloatList(value)
            elif self.types[name] == "function":
                self.values[name] = value
            elif self.types[name] == "general":
                self.values[name] = value
            elif self.types[name] == "parsed_function":
                self.values[name] = stringToFunction(value)
            elif self.types[name] == "string":
                self.values[name] = value
            elif self.types[name] == "StringSelection":
                if value in self.string_selections[name]:
                    self.values[name] = value
                else:
                    selection_list_string = ""
                    for acceptable_value in self.string_selections[name]:
                        selection_list_string += "\n  " + acceptable_value
                    error(
                        "The parameter '" + name + "' must take one of the following values:" +
                        selection_list_string)
        else:
            error("The parameter '" + name + "' is not registered.")

    def registerParameterInternal(self, name, param_type, description, default):
        self.types[name] = param_type
        if name == "type":
            error("'type' is a reserved name and cannot be used for a parameter name.")
        if name in self.descriptions:
            error("The parameter '" + name + "' has already been registered.")
        else:
            self.descriptions[name] = description
            if default != None:
                self.values[name] = default

    def registerStringListParameter(self, name, description, default=None):
        self.registerParameter(name, description, default)
        self.types[name] = "string_list"

    def registerBoolParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "bool", description, default)

    def registerIntParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "int", description, default)

    def registerFloatParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "float", description, default)

    def registerFloatListParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "float_list", description, default)

    def registerFunctionParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "function", description, default)

    def registerParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "general", description, default)

    def registerParsedFunctionParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "parsed_function", description, default)

    def registerStringParameter(self, name, description, default=None):
        self.registerParameterInternal(name, "string", description, default)

    ## Registers a string parameter that has a selection of values
    # @param name  name of the parameter
    # @param selection  list of acceptable values
    # @param description  description of the parameter
    # @param default  optional default value
    def registerStringSelectionParameter(self, name, selection, description, default=None):
        self.registerParameterInternal(name, "StringSelection", description, default)
        self.string_selections[name] = selection
