import unittest

from sem_python.base.enums import ModelType
from sem_python.base.Factory import Factory


class FEMDoFHandlerTester(unittest.TestCase):
    ## Tests the DoF indices assigned to junction constraints
    #
    # 1-phase: 3 variables
    #
    # mesh1: 3 cells
    # mesh2: 5 cells
    # mesh3: 4 cells
    # mesh4: 1 cell
    # junction1(mesh1, mesh2): 2 constraints
    # junction2(mesh2, mesh3): 6 constraints
    # junction3(mesh1, mesh2): 3 constraints
    # junction4(mesh3, mesh4): 1 constraint
    #
    # The resulting ordering is [mesh1, junction1, junction3, mesh2, junction2, mesh3, junction4, mesh4]
    #
    # With all of this information, the correct DoF indices for each junction are:
    # junction1: [12, 13]
    # junction2: [35, 36, 37, 38, 39, 40]
    # junction3: [14, 15, 16]
    # junction4: [56]
    #
    def testJunctionConstraintDoFIndices(self):
        # create the factory
        factory = Factory()

        # create the meshes
        n_cells_list = [3, 5, 4, 1]
        meshes = list()
        for i, n_cells in enumerate(n_cells_list):
            name = "mesh" + str(i + 1)
            params = {"name": name, "n_cell": n_cells}
            meshes.append(factory.createObject("UniformMesh", params))

        # ICs
        ics = [
            factory.createObject(
                "InitialConditions1Phase", {
                    "mesh_name": meshes[0].name,
                    "A": "1",
                    "rho": "1",
                    "u": "1",
                    "p": "1"
                })
        ]

        # create the DoF handler
        params = {"meshes": meshes, "ics": ics, "model_type": ModelType.OnePhase}
        dof_handler = factory.createObject("FEMDoFHandler", params)

        # create an EoS
        eos_list = [factory.createObject("TestEoS")]

        # create the junctions
        n_constraints_list = [2, 6, 3, 1]
        n_junctions = len(n_constraints_list)
        meshes_list = [
            ["mesh1", "mesh2"], ["mesh2", "mesh3"], ["mesh1", "mesh2"], ["mesh3", "mesh4"]
        ]
        sides_list = [["right", "left"]] * n_junctions
        junctions = list()
        for i in range(n_junctions):
            params = {
                "mesh_names": meshes_list[i],
                "mesh_sides": sides_list[i],
                "dof_handler": dof_handler,
                "eos_list": eos_list,
                "n_constraints": n_constraints_list[i]
            }
            junctions.append(factory.createObject("TestJunction", params))

        # update the DoF handler with the junction constraints
        dof_handler.updateWithJunctionConstraints(junctions)

        # check all of the constraint DoF indices are the expected
        expected_constraint_dof_indices = [[12, 13], [35, 36, 37, 38, 39, 40], [14, 15, 16], [56]]
        for i, junction in enumerate(junctions):
            self.assertEqual(junction.i_constraint, expected_constraint_dof_indices[i])
