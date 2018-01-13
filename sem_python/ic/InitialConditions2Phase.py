from ..input.Parameters import Parameters
from ..utilities.error_utilities import error


class InitialConditions2PhaseParameters(Parameters):

    def __init__(self):
        Parameters.__init__(self)

        def one(x):
            return 1

        self.registerStringParameter("mesh_name", "Name of the mesh to which this corresponds")
        self.registerParsedFunctionParameter("A", "Cross-sectional area of flow channel", one)
        self.registerParsedFunctionParameter("vf1", "Volume fraction of phase 1")
        self.registerParsedFunctionParameter("rho1", "Density of phase 1")
        self.registerParsedFunctionParameter("p1", "Pressure of phase 1")
        self.registerParsedFunctionParameter("T1", "Temperature of phase 1")
        self.registerParsedFunctionParameter("u1", "Velocity of phase 1")
        self.registerParsedFunctionParameter("rho2", "Density of phase 2")
        self.registerParsedFunctionParameter("p2", "Pressure of phase 2")
        self.registerParsedFunctionParameter("T2", "Temperature of phase 2")
        self.registerParsedFunctionParameter("u2", "Velocity of phase 2")


class InitialConditions2Phase(object):

    def __init__(self, params):
        self.mesh_name = params.get("mesh_name")
        self.A = params.get("A")
        self.vf1 = params.get("vf1")
        self.p = [params.get("p1"), params.get("p2")]
        self.u = [params.get("u1"), params.get("u2")]

        # one may supply either rho or T, but not both
        if params.has("rho1") and params.has("rho2"):
            has_rho = True
        else:
            has_rho = False
        if params.has("T1") and params.has("T2"):
            has_T = True
        else:
            has_T = False
        if (params.has("T1") and params.has("rho2")) or (params.has("rho1") and params.has("T2")):
            error("ICs cannot supply a mix of T and rho between the phases.")
        elif has_rho and has_T:
            error("ICs cannot specify both T and rho.")
        elif has_rho:
            self.rho = [params.get("rho1"), params.get("rho2")]
            self.specified_rho = True
        elif has_T:
            self.T = [params.get("T1"), params.get("T2")]
            self.specified_rho = False
        else:
            error("Either 'rho' or 'T' for each phase must be specified for IC.")
