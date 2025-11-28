import time as tt
import numpy as np

class time:
    def __init__(self, execution):
        # ----- Public attributes ----- #
        self.name = "TIME" # Name is harcoded for time objects
        self.control = execution['time']
        self.path = execution['paths']['executionPath']
        self.stepping = execution['time'].get('stepping', 'automatic')
        self.start = execution['time']['startTime']
        self.end = execution['time']['endTime']
        self.delta = execution['time']['deltaT']
        self.value = self.start  # Current time
        self.startDate = tt.perf_counter()
        self.span = np.array([self.start, self.start + self.delta])
        self.pf = 0
        # Output variable mapping
        self.varMap = {
            "time":       "value"
        }

    def elapsed(self):
        return tt.perf_counter() - self.startDate

    def advance(self, dt=None):
        if dt is not None:
            self.delta = dt
        self.value += self.delta
        self.span += self.delta
        self.pf = (self.value - self.start) / (self.end - self.start)
