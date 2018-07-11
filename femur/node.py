class Node:
    def __init__(self, id, x, y, z):
        self.id = int(id) - 1
        self.x = x
        self.y = y
        self.z = z
        self.is_noslip = False
        self.is_input = False
        self.is_output = False

    def get_position(self):
        return [self.x, self.y, self.z]

    def __eq__(self, other):
        return self.get_position() == other.get_position()
