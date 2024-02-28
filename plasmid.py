from module import Module
from part import Part
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO


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
        self.marker_module = None  #Three placeholder variables to just provide template structure
        self.origin_module = None

        #Creating the non-variable modules on construction of an object
        self.t1_part = Part("t1_part","t1_part","CAGCTGTCTAGGGCGGCGGATTTGTCCTACTCAGGAGAGCGTTCACCGACAAACAACAGATAAAACGAAAGGCCCAGTCTTTCGACTGAGCCTTTCGTTTTATTTGATGCCT","terminator")
        self.t1 = Module("t1","t1_module","invariant", structure =[self.t1_part])

        self.t0_part = Part("t0_part","t0_part","CTTGGACTCCTGTTGATAGATCCAGTAATGACCTCAGAACTCCATCTGGATTTGTTCAGAACGCTCGGTTGCCGCCGGGCGTTTTTTATTGGTGAGAATCCAG","terminator")
        self.SanDI = Part("SanDI_part","sandi","GGGTCCC","re")
        self.t0 = Module("t0","t0_module","invariant", structure =[self.t0_part,self.SanDI])

        self.orit_part = Part("orit_part","orit_part","CTTTTCCGCTGCATAACCCTGCTTCGGGGTCATTATAGCGATTTTTTCGGTATATCCATCCTTTTTCGCACGATATACAGGATTTTGCCAAAGGGTTCGTGTAGACTTTCCTTGGTGTATCCAACGGCGTCAGCCGGGCAGGATAGGTGAAGTAGGCCCACCCGCGAGCGGGTGTTCCTTCTTCACTGTCCCTTATTCGCACCTGGCGGTGCTCAACGGGAATCCTGCTCTGCGAGGCTGGCCGTA","oriT")
        self.orit = Module("orit","orit_module","invariant", structure =[self.orit_part])
        
        # Skip first validation on initilisation so ._
        self._structure = [self.cargo_module,self.t0,
                           self.marker_module,self.orit,
                           self.origin_module,self.t1
                           ] 

    @property
    def name(self):
        """Gets the name attribute"""
        return self._name

    @name.setter
    def name(self, name):
        """Sets the name attribute, type: String"""
        if not isinstance(name, str):
            raise TypeError("Name must be a string")
        self._name = name

    @property
    def unique_id(self):
        """Sets the unique_id attribute"""
        return self._unique_id

    @unique_id.setter
    def unique_id(self, id):
        """Sets the unique_id attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id


    @property
    def structure(self):
        """Gets the structure attribute"""
        return self._structure

    @structure.setter
    def structure(self, structure):
        """Sets the name attribute
            Input: A list of three individual seperate modules, type Cargo,Marker,Origin"""
        if not isinstance(structure, list):
            raise TypeError("Input structure must be of Type: list")
        for modules in structure:
            if not isinstance(modules, Module):
                raise TypeError("List should contain Module types only")
        
         #Creates a list containing the module types and compares it with predefined structure to see if its correct
        correct_val = ["cargo","marker","origin"]
        order = []
        for i in structure:
            order.append(i.module_type)
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


    def insert_scar(self, scar, module = None, position = None):
        """"Method to insert scars""" 
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



    def gather_data(self): 
        """Method which gathers the data of each part in the plasmid and puts it all into a dictionary with key being name and values being other data attributes and locations"""
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


    def serialise_genbank(self,path):
        """Serialises in genbank format: Need input file path"""
        sequence = self.get_sequence()
        sequence_r = SeqRecord(Seq(sequence), id = self.unique_id, name = self.name, description = self.description )
        print(sequence_r)
        data_dict = self.gather_data() #Obtains data dictionary and uses this to add features of each part

        conversion = {
                    're': 'misc_feature',  #Dictionary to convert the role into the GenBank accepted role
                    'cargo': 'misc_feature',
                    'abr': 'misc_feature',
                     'ori': 'rep_origin'
        }


        for part , attribute in data_dict.items(): # Role Conversion 
            if attribute[2] in conversion:
                attribute[2] = conversion[attribute[2]]
       
            
        sequence_r.annotations["molecule_type"] = "DNA" #Annotate the sequence to show that its DNA

        #Going through the data dictionary an using it to annotate each part of the sequence
        for part, attribute in data_dict.items():
            feature = SeqFeature(
                FeatureLocation(attribute[0], attribute[4]), type= attribute[2],
                qualifiers={
                    "gene": [part],  # Gene name
                    "locus_tag": [attribute[1]],  # A unique identifier for the gene
                    "note": [attribute[3]]  # Any additional notes
                }
            )
            sequence_r.features.append(feature)

        
        # Write the SeqRecord to a GenBank file
        with open(path, "w") as output_handle:
            SeqIO.write(sequence_r, output_handle, "genbank")

        print(f"GenBank file saved at {path}")
        print(sequence_r) #Help show that the serialisation has worked