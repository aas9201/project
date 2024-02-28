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
        """gets the name attribute, type: String"""
        if not isinstance(id, str):
            raise TypeError("ID must be a string")
        self._unique_id = id

    @property
    def sequence(self):
        """gets the sequence attribute"""
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """sets the sequence attribute, type: String, Only contains A,C,T or G with no spaces"""
        if not isinstance(sequence, str):
            raise TypeError("Name must be a string")
        if not all(c in "ACTG" for c in sequence) or " " in sequence:
            raise ValueError("Sequence must only contain A,C,T,G with no spaces")
        self._sequence = sequence

    @property
    def role(self):
        """gets the role attribute"""
        return self._role

    @role.setter
    def role(self, role):
        """sets the role attribute, type: String, must only be certain roles also"""
        if not isinstance(role, str):
            raise TypeError("Role must be a string")
        if role not in ["gene", "abr", "ori", "terminator","promoter","scar","re","oriT", "cargo"]:
            raise ValueError("Role must be either gene, abr, ori, terminator, promoter, scar, re, oriT, cargo")


        self._role = role

    @property
    def description(self):
        """gets the description attribute"""
        return self._description

    @description.setter
    def description(self, description):
        """sets the description attribute, type: String"""
        if not isinstance(description, str):
            raise TypeError("Description must be a string")
        self._description = description