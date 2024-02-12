class Part:
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
    def name(self, value):
        if not isinstance(value, str):
            raise TypeError("Name must be a string")
        self._name = value

    @property
    def unique_id(self):
        return self._unique_id

    @unique_id.setter
    def unique_id(self, value):
        if not isinstance(value, str):
            raise TypeError("ID must be a string")
        self._unique_id = value

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if not isinstance(value, str):
            raise TypeError("Name must be a string")
        if not all(c in "ACTG" for c in value) or " " in value:
            raise ValueError("Sequence must only contain A,C,T,G with no spaces")
        self._sequence = value

    @property
    def role(self):
        return self._role

    @role.setter
    def role(self, value):
        if not isinstance(value, str):
            raise TypeError("Role must be a string")
        self._role = value

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, value):
        if not isinstance(value, str):
            raise TypeError("Description must be a string")
        self._description = value


class Module:

    def __init__(self, name, unique_id, module_type, description=""):
        self.name = name
        self.unique_id = unique_id
        self.module_type = module_type
        self.description = description
        self.structure = []

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, str):
            raise TypeError("Name must be a string")
        self._name = value

    @property
    def unique_id(self):
        return self._unique_id

    @unique_id.setter
    def unique_id(self, value):
        if not isinstance(value, str):
            raise TypeError("ID must be a string")
        self._unique_id = value

    @property
    def module_type(self):
        return self._module_type

    @module_type.setter
    def module_type(self, value):
        if value not in ["cargo", "marker", "origin"]:
            raise ValueError("Module must be either a cargo, marker or origin")
        self._module_type = value

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, value):
        if not isinstance(value, list):
            raise TypeError("Input structure must be of Type: list")
        for i in value:
            if not isinstance(i, Part):
                raise TypeError("List should contain Part types only")
        self._structure = value

    def get_sequence(self):
        sequence = ""
        for i in self.structure:
            sequence += i.sequence
        return sequence


class Plasmid:

    def __init__(self, name, unique_id, description=""):
        self.name = name
        self.unique_id = unique_id
        self.description = description
        self.structure = []

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        if not isinstance(value, str):
            raise TypeError("Name must be a string")
        self._name = value

    @property
    def unique_id(self):
        return self._unique_id

    @unique_id.setter
    def unique_id(self, value):
        if not isinstance(value, str):
            raise TypeError("ID must be a string")
        self._unique_id = value

    @property
    def module_type(self):
        return self._module_type

    @module_type.setter
    def module_type(self, value):
        if value not in ["cargo", "marker", "origin"]:
            raise ValueError("Module must be either cargo, marker or origin")
        self._module_type = value

    @property
    def structure(self):
        return self._structure

    @structure.setter
    def structure(self, value):
        if not isinstance(value, list):
            raise TypeError("Input structure must be of Type: list")
        for i in value:
            if not isinstance(i, Module):
                raise TypeError("List should contain Module types only")
        self._structure = value

    def get_sequence(self):
        plasmid_sequence = ""
        for i in self.structure:
            plasmid_sequence += i.get_sequence()
        return plasmid_sequence


promoter = Part("promoter", "1", "AAA", "pro")
cds = Part("cds", "2", "TTT", "cds")
terminator = Part("terminator", "3", "GGG", "ter")

cargo = Module("alpha", "4", "cargo")
cargo_list = [promoter, cds, terminator]
cargo.structure = cargo_list
sequence = cargo.get_sequence()
print(sequence)

promoter1 = Part("promoter", "1", "GGG", "pro")
cds1 = Part("cds", "2", "TTT", "cds")
terminator1 = Part("terminator", "3", "AAA", "ter")

abr = Module("alpha", "4", "marker")
abr_list = [promoter1, cds1, terminator1]
abr.structure = abr_list
sequence1 = abr.get_sequence()
print(sequence1)

Plasmid1 = Plasmid("plasmid", "2")
plasmid_list = [cargo, abr]
Plasmid1.structure = plasmid_list
plasmid_sequence = Plasmid1.get_sequence()
print(plasmid_sequence)
