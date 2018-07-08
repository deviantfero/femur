class Node:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.is_noslip = False
        self.is_input = False
        self.is_output = False

    def get_position(self):
        return [self.x, self.y, self.z]
