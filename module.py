from part import Part

class Module:
    """Class represents the most module componenet of a plasmid:
        Data attributes:
        Name: String
        Unique ID: String
        Module Type: String (either Cargo, Marker, Origin or Invariant)
        Structure: List of Part objects
        Description: String, add some information about the part
        Fetched: True or False, allows for different structure construction if sequence is imported"""
       
    def __init__(self, name, unique_id, module_type, description="",structure =[], fetched = False):
        self.name = name
        self.unique_id = unique_id
        self.module_type = module_type
        self.description = description
        self.fetched = fetched
        self.structure = structure



    @property
    def name(self):
        """gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """gets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """sets the unique_id attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def module_type(self):
        """gets the module_type attribute"""
        return self._module_type

    @module_type.setter
    def module_type(self, module):
        """sets the unique_id attribute"""
        if module not in ["cargo", "marker", "origin", "invariant"]:
            raise ValueError("Module must be either a cargo, marker, origin, invariant")
        self._module_type = module

    @property
    def structure(self):
        """gets the structure attribute"""
        return self._structure

    @structure.setter
    def structure(self, structure_list):
        """Function which sets the structure
        Input: list of part(s) objects
        Depending on module_type a predefined structure with corresponding flanking RE sites are assigned"""
        

        if not isinstance(structure_list, list):
            raise TypeError("Input structure must be of Type: list")
        for parts in structure_list:
            if not isinstance(parts, Part):
                raise TypeError("List should contain Part types only")
        
        if self.module_type == "cargo": #Creates Part objects for flanking RE sites dynamically
            self.cargo = None
            self.PacI = Part("PacI","paci","TTAATTAA","re") 
            self.SpeI = Part("SpeI","spei","ACTAGT","re")
            self._structure = [self.PacI,self.cargo,self.SpeI]   
            self._structure[1:2] = structure_list    #The input part list is inserted by a slice to remove the placeholder


        if self.module_type == "origin":
            self.replication = None
            self.FseI = Part("FseI","fsei","GGCCGGCC","re")
            self.AscI = Part("AscI","asci","GGCGCGCC","re")
            self._structure = [self.FseI,self.replication,self.AscI]
            self._structure[1:2] = structure_list #The input part list is inserted by a slice to remove the placeholder


        if self.module_type == "marker" and self.fetched == False:
            self.marker = None
            self.PshAI = Part("pshai","1","GACGTC","re") ##Need to sort this logic out a bit more
            self.SwaI = Part("swai","1","ATTTAAAT","re")
            self._structure = [self.SwaI,self.marker,self.PshAI]
            self._structure[1:2] = structure_list #The input part list is inserted by a slice to remove the placeholder

        if self.module_type == "marker" and self.fetched == True:
            self._structure = structure_list #If sequence is imported I want to assign PshAI dynamically as certain SEVA plasmids contain different sequence

        if self.module_type =="invariant":
            self._structure = structure_list





    def get_sequence(self):
        """Obtains a sequence by concatenating sequence of parts together"""
        sequence = ""
        for i in self.structure:
            sequence += i.sequence
        return sequence