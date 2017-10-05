import numpy as np

from ..base.enums import ModelType
from .DoFHandler import DoFHandler, DoFHandlerParameters

class DoFHandler2PhaseNonInteractingParameters(DoFHandlerParameters):
  def __init__(self):
    DoFHandlerParameters.__init__(self)

class DoFHandler2PhaseNonInteracting(DoFHandler):
  def __init__(self, params):
    DoFHandler.__init__(self, params)
    self.model_type = ModelType.TwoPhaseNonInteracting
    self.n_phases = 2
    self.n_vf_equations = 0
    self.setup()

    # create array for volume fraction
    self.vf1 = np.zeros(self.n_node)
    for ic in self.ics:
      # get corresponding mesh
      mesh_name = ic.mesh_name
      i_mesh = self.mesh_name_to_mesh_index[mesh_name]
      mesh = self.meshes[i_mesh]

      # compute volume fraction for each node on mesh
      vf1 = ic.vf1
      for k_mesh in xrange(mesh.n_node):
        k = self.k_from_k_mesh(k_mesh, i_mesh)
        self.vf1[k] = vf1(mesh.x[k_mesh])

  def aA1(self, U, k):
    return self.vf1[k] * self.A[k]
