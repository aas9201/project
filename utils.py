from part import Part
from module import Module
import re



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
        cargo_part = Part("cargo_part", "cargo", cargo_sequence, "cargo")
        cargo_module = Module(name, ID, "cargo", description, [cargo_part])
    
        return cargo_module

    if module_type == "marker": #Slices the plasmid between restriction enzyme site to get Marker module
        marker_sequence = (plasmid_sequence[(sorted_sites["SwaI"][0]+7):sorted_sites["PshAI"][0]-1])
        marker_part = Part("module_part", "module_part", marker_sequence, "abr")

        marker_module = Module(name, ID, "marker", description, fetched = True)
        marker_module.SwaI = Part("SwaI","swai","ATTTAAAT","re")
        pshai = (plasmid_sequence[(sorted_sites["PshAI"][0]-1):sorted_sites["PshAI"][0]+9]) #As the PshAI contains NNN this needs to be created 
        marker_module.PshAI = Part("PshAI","pshai",pshai,"re")
        marker_module.structure = [marker_module.SwaI,marker_part,marker_module.PshAI]

        return marker_module
    
    if module_type == "origin": #Slices the plasmid between restriction enzyme site to get Origin module
        origin_sequence = (plasmid_sequence[(sorted_sites["FseI"][0]+7):sorted_sites["AscI"][0]-1])
        origin_part = Part("origin_part", "origin_part", origin_sequence, "ori")
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