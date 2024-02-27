import re

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


class Part:
    """Class represents the most fundamental componenet of a plasmid:
    Data attributes:
    Name: String
    Unique ID: String
    Sequence: String (only ACTG and no spaces)
    Role: String, gives function for further info and serialiastion
    Description: String, add some information about the part"""

    def __init__(self, name, unique_id, sequence, role, description=""):
        self.name = name
        self.unique_id = unique_id
        self.sequence = sequence
        self.role = role
        self.description = description

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        if not isinstance(sequence, str):
            raise TypeError("Name must be a string")
        if not all(c in "ACTG" for c in sequence) or " " in sequence:
            raise ValueError("Sequence must only contain A,C,T,G with no spaces")
        self._sequence = sequence

    @property
    def role(self):
        return self._role

    @role.setter
    def role(self, role):
        if not isinstance(role, str):
            raise TypeError("Role must be a string")
        self._role = role

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, description):
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description


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
        return self._name

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def module_type(self):
        return self._module_type

    @module_type.setter
    def module_type(self, module):
        if module not in ["cargo", "marker", "origin", "invariant"]:
            raise ValueError("Module must be either a cargo, marker, origin, invariant")
        self._module_type = module

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, structure_list):
        """Function which sets the structure depending on the module type"""

        if not isinstance(structure_list, list):
            raise TypeError("Input structure must be of Type: list")
        for parts in structure_list:
            if not isinstance(parts, Part):
                raise TypeError("List should contain Part types only")
        
        if self.module_type == "cargo":
            self.cargo = None
            self.PacI = Part("PacI","1","TTAATTAA","misc_feature")
            self.SpeI = Part("SpeI","1","ACTAGT","misc_feature")
            self._structure = [self.PacI,self.cargo,self.SpeI]
            self._structure[1:2] = structure_list


        if self.module_type == "origin":
            self.replication = None
            self.FseI = Part("FseI","1","GGCCGGCC","misc_feature")
            self.AscI = Part("AscI","1","GGCGCGCC","misc_feature")
            self._structure = [self.FseI,self.replication,self.AscI]
            self._structure[1:2] = structure_list


        if self.module_type == "marker" and self.fetched == False:
            self.marker = None
            self.PshAI = Part("PshAI","1","GACGTC","misc_feature") ##Need to sort this logic out a bit more
            self.SwaI = Part("SwaI","1","ATTTAAAT","misc_feature")
            self._structure = [self.SwaI,self.marker,self.PshAI]
            self._structure[1:2] = structure_list

        if self.module_type == "marker" and self.fetched == True:
            self._structure = structure_list

        if self.module_type =="invariant":
            self._structure = structure_list





    def get_sequence(self):
        """Obtains a sequence by dynamicaly concatenating sequence of parts together"""
        sequence = ""
        for i in self.structure:
            sequence += i.sequence
        return sequence


class Plasmid:
    """Class which builds final plasmid:
        Data attributes:
        Name: String
        Unique ID: String
        Module Type: String (either Cargo, Marker, Origin or Invariant)
        Structure: List of Module objects
        Description: String, add some information about the part
        """
    
    def __init__(self, name, unique_id, description=""):
        self.name = name
        self.unique_id = unique_id
        self.description = description

        self.cargo_module = None
        self.marker_module = None
        self.origin_module = None


        self.t1_part = Part("t1_part","1","CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT","terminator")
        self.t1 = Module("t1","1","invariant", structure =[self.t1_part])

        self.t0_part = Part("t0_part","1","CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG","terminator")
        self.SanDI = Part("SanDI_part","1","GGGTCCC","misc_feature")
        self.t0 = Module("t0","1","invariant", structure =[self.t0_part,self.SanDI])

        self.orit_part = Part("orit_part","1","CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA","oriT")
        self.orit = Module("orit","1","invariant", structure =[self.orit_part])

        self._structure = [self.cargo_module,self.t0,self.marker_module,self.orit,self.origin_module,self.t1] # On initilisation I want to skip first validation

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id


    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, structure):
        if not isinstance(structure, list):
            raise TypeError("Input structure must be of Type: list")
        for modules in structure:
            if not isinstance(modules, Module):
                raise TypeError("List should contain Module types only")
        
        
        correct_val = ["cargo","marker","origin"]
        order = []
        for i in structure:
            order.append(i.module_type) #Creates a list containing the module types and compares it with predefined structure to see if its correct
        if correct_val != order:
            raise ValueError("Please ensure the modules are in the correct order: [Cargo,Marker,Origin]")



        iterVal = iter(structure)
        for parts, element in enumerate(self._structure):
            if element is None: #Finds which variables are "None" in structure AKA placeholders and then swaps them for each corresponding module
                self._structure[parts] = next(iterVal)
            

    def get_sequence(self): #it iterates through each module in structure and then calls their get_sequence function and concatenates all together
        plasmid_sequence = ""
        for module in self.structure:
            plasmid_sequence += module.get_sequence()
        return plasmid_sequence


    def insert_scar(self, scar, module = None, position = None): #Method to insert scars 
        module_dict = {
        "T1": self.t1,
        "T0": self.t0,
        "OriT": self.orit
    }
    
        module = module_dict.get(module) # Get the module object based on the module name
    
        if module: #Inserts the scar depending on if it comes before or after the module
           if position == "before":
               module.structure.insert(0, scar)
           elif position == "after":
               module.structure.append(scar)



    def gather_data(self): #Method which gathers the data of each part in the plasmid and puts it all into a dictionary with key being name and values being other data attributes and locations
        sequence = self.get_sequence()
        positions = []
        end_point = []
        my_dict = {}


        for module in self.structure:
            for part in module.structure:
                position = sequence.find(part.sequence) #Finds where the part occurs in the sequence and adds it to dictionary along with other data attributes
                positions.append(position) 
                my_dict[part.name] = [position]
                my_dict[part.name].append(part.unique_id)
                my_dict[part.name].append(part.role)
                my_dict[part.name].append(part.description)

        for i in range(len(positions)-1):     #Looks at position of where next part is and then -1 as that will be the endpoint
            end_point.append(positions[i+1] - 1)
        end_point.append(len(sequence))

        for part, (key, value) in enumerate(my_dict.items()):
            # Update each key in my_dict with its corresponding end point
            if part < len(end_point):
                my_dict[key].append(end_point[part])

        return my_dict


    def serialise_genbank(self,path): #Serialises in genbank fomrat
        sequence = self.get_sequence()
        sequence_r = SeqRecord(Seq(sequence), id = self.unique_id, name = self.name, description = self.description )
        print(sequence_r)
        data_dict = self.gather_data() #Obtains data dictionary and uses this to add features of each part
        sequence_r.annotations["molecule_type"] = "DNA"

        for key, values in data_dict.items():
            feature = SeqFeature(
                FeatureLocation(values[0], values[4]), type= values[2],
                qualifiers={
                    "gene": [key],  # Gene name
                    "locus_tag": [values[1]],  # A unique identifier for the gene
                    "note": [values[3]]  # Any additional notes
                }
            )
            sequence_r.features.append(feature)

        
        # Write the SeqRecord to a GenBank file
        with open(path, "w") as output_handle:
            SeqIO.write(sequence_r, output_handle, "genbank")

        print(f"GenBank file saved at {path}")
        print(sequence_r)




def sort_sites(plasmid_sequence):  #Sorts the plasmid using the non-variable sites , my class logic is a bit different now so i am going to change this a bit as some code is redundant 
    predefined_sites = {
        "OriT": "CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA",
        "T0":"CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG",
        "T1": "CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT",
        "PshAI": "GACNNNNGTC",  
        "SwaI": "ATTTAAAT",
        "AscI": "GGCGCGCC",
        "FseI": "GGCCGGCC",
        "PacI": "TTAATTAA",
        "SpeI": "ACTAGT",
        "SanDI": "GGGTCCC",     
             }
        
    correct_order = ["PacI", "SpeI", "T0", "SanDI", "SwaI", "PshAI", "OriT", "FseI", "AscI", "T1"]
    site_positions = {}
    first_occurance = []

    for site_name, site_sequence in predefined_sites.items():
        pattern = re.sub('N', '[ACGT]', site_sequence)  #Replaces "N" for either ACTG
        positions = []

        for i in range(len(plasmid_sequence)):
            if re.match(pattern, plasmid_sequence[i:]):
                positions.append(i+1) #Appends the position but +1 as starts from 0
                
        if positions:
            site_positions[site_name] = positions #Stores position at the current site_name
            first_occurance.append((site_name, positions[0])) #Sorts by where they first occur 

        first_occurance.sort(key=lambda x: x[1]) #Order each site into ascending postions
    
    sorted_sites = {site: site_positions[site] for site, _ in first_occurance} #Create a new dictionary which has sites ordered
    
    for i in sorted_sites:
        if len(sorted_sites[i]) > 1:
            raise ValueError(f"The plasmid contains more than one: {i}") #Checks if there is only one site
        
    plasmid_order = list(sorted_sites.keys())
    if plasmid_order != correct_order:
        raise ValueError("Please make sure plasmid is of valid SEVA Format") #Checks order is correct
    
    return sorted_sites
    
        

def get_module(plasmid_sequence, name, ID, module_type, description =""):
    """Creates a new module object by slicing the plasmid sequence depending on which module_type you want"""
    sorted_sites = sort_sites(plasmid_sequence) #Sorts the plasmid so that it can be sliced easily

    if module_type == "cargo": #Slices the plasmid between restriction enzyme site to get Cargo module
        cargo_sequence = (plasmid_sequence[(sorted_sites["PacI"][0]+7):sorted_sites["SpeI"][0]-1])
        cargo_part = Part("cargo_part", "1", cargo_sequence, "misc_feature")
        cargo_module = Module(name, ID, "cargo", description, [cargo_part])
    
        return cargo_module

    if module_type == "marker": #Slices the plasmid between restriction enzyme site to get Marker module
        marker_sequence = (plasmid_sequence[(sorted_sites["SwaI"][0]+7):sorted_sites["PshAI"][0]-1])
        marker_part = Part("module_part", "1", marker_sequence, "misc_feature")

        marker_module = Module(name, ID, "marker", description, fetched = True)
        marker_module.SwaI = Part("SwaI","1","ATTTAAAT","misc_feature")
        pshai = (plasmid_sequence[(sorted_sites["PshAI"][0]-1):sorted_sites["PshAI"][0]+9]) #As the PshAI contains NNN this needs to be created 
        marker_module.PshAI = Part("PshAI","1",pshai,"misc_feature")
        marker_module.structure = [marker_module.SwaI,marker_part,marker_module.PshAI]

        return marker_module
    
    if module_type == "origin": #Slices the plasmid between restriction enzyme site to get Origin module
        origin_sequence = (plasmid_sequence[(sorted_sites["FseI"][0]+7):sorted_sites["AscI"][0]-1])
        origin_part = Part("origin_part", "1", origin_sequence, "rep_origin")
        origin_module = Module(name, ID, "origin", description, [origin_part])
        
        return origin_module

    

def get_scars(plasmid_sequence): #Slices plasmid between certain sites to check for any potential scars as each SEVA plasmid has certain scars, some don't
    sorted_sites = sort_sites(plasmid_sequence)
    scars = {}
    check_scars = {"Scar after T1": (plasmid_sequence[0:sorted_sites["PacI"][0]-1]),
                   "Scar before T0": (plasmid_sequence[(sorted_sites["SpeI"][0]+5):sorted_sites["T0"][0]-1]),
                   "Scar after T0": (plasmid_sequence[(sorted_sites["SanDI"][0]+6):sorted_sites["SwaI"][0]-1]),
                   "Scar before OriT": (plasmid_sequence[(sorted_sites["PshAI"][0]+9):sorted_sites["OriT"][0]-1]),
                   "Scar after OriT": (plasmid_sequence[(sorted_sites["OriT"][0]+245):sorted_sites["FseI"][0]-1]),
                   "Scar before T1": (plasmid_sequence[(sorted_sites["AscI"][0]+7):sorted_sites["SwaI"][0]-1]),
                   }
    
    for key, value in check_scars.items():
        if value:
            scars[key] = value
    print(scars)
    return scars




