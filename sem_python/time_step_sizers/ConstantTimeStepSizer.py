from .TimeStepSizer import TimeStepSizer, TimeStepSizerParameters


class ConstantTimeStepSizerParameters(TimeStepSizerParameters):

    def __init__(self, factory):
        TimeStepSizerParameters.__init__(self, factory)
        self.registerFloatParameter("dt", "Time step size")


class ConstantTimeStepSizer(TimeStepSizer):

    def __init__(self, params):
        TimeStepSizer.__init__(self, params)
        self.dt = params.get("dt")

    def getTimeStepSizeInternal(self, U):
        return self.dt
