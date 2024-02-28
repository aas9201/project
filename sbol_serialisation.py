from sbol2 import ComponentDefinition, Sequence, Document, setHomespace
from sbol2 import Document
from plasmid import Plasmid

setHomespace('http://project7.com')
doc = Document()
doc.clear()


def adding_component(displayid = str, ontology = str ):
    """Creates SBOL component"""
    comp = ComponentDefinition(displayid)
    comp.roles = f'http://identifiers.org/so/SO:{ontology}'
    doc.addComponentDefinition(comp)
    return comp

def assignsequence(id, code = str):
    """Assigns the sequence to SBOL Component"""
    seq = Sequence(id, code)
    return seq



def serialise_sbol(instance,path):
    """Function which serialises in SBOL format
    Input: Plasmid object and file path"""
    if not isinstance(instance, Plasmid):
        raise TypeError("Please input a plasmid object")
    
    doc.clear()
    primary_structure = []
    data_dict = instance.gather_data()
    conversion = {
                    're': '0001687',  # Dictionary to convert roles into SBOL ontology codes
                    'cargo': '0000704',
                    'abr': '0000001',
                    'terminator': '0000141',
                    'oriT': '0000724',
                    'ori': '0000296',
                    'scar': '0001953',
        }

    # Apply conversion
    for part, attribute in data_dict.items():
        if attribute[2] in conversion:
            attribute[2] = conversion[attribute[2]]
    #Add sequence to dictionary as need to add component and sequence in one loop for it to work
    for module in instance.structure:
        for part in module.structure:
            data_dict[part.name].append(part.sequence)

    for part, attibute in data_dict.items():
        component = adding_component(attibute[1], attibute[2]) #Add Sbol component
        component.sequence = assignsequence(component.displayId,attibute[5]) #Add Sequence to component
        primary_structure.append(component) #Add the component to primary structure 

    instance_built = adding_component(instance.name,"0000155") #Create final plasmid component for hierarchy
    instance_built.assemblePrimaryStructure(primary_structure) #Assign primary structure to plasmid component
    instance_built.compile()
    doc.write(path)
