class ProblemData:
    def __init__(self, success, time, coms=None, moving_foot_pos=None, all_feet_pos=None, surface_indices=None):
        self.success = success
        self.time = time
        self.coms = coms
        self.moving_foot_pos = moving_foot_pos
        self.all_feet_pos = all_feet_pos
        self.surface_indices = surface_indices

    def __str__(self):
        string = "ProblemData: "
        string += "\n \t Success: " + str(self.success)
        string += "\n \t Time: " + str(self.time)
        return string

    def __repr__(self):
        return self.__str__()
