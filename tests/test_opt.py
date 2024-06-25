import pytest

from libsyn_tools.opt.formulation_milp import SolverMILP
from libsyn_tools.opt.schema import SchedulerInput, random_functional_modules
from .test_chem_schema import test_reaction_network1, OperationNetwork, test_reaction_network2, StateOfMatter


@pytest.fixture
def operation_network1(test_reaction_network1):
    test_reaction_network1.dummy_quantify()
    network1 = OperationNetwork.from_reaction_network(test_reaction_network1)
    return network1


@pytest.fixture
def operation_network2(test_reaction_network2):
    test_reaction_network2.dummy_quantify()
    # LIQUID fix
    for r in test_reaction_network2.chemical_reactions:
        if len(r.reactants + r.reagents) == 1:
            r.reactants[0].state_of_matter = StateOfMatter.LIQUID
    network2 = OperationNetwork.from_reaction_network(test_reaction_network2)
    return network2


@pytest.fixture
def functional_modules1():
    return random_functional_modules(random_seed=42, max_capacity=1, max_module_number_per_type=1)


@pytest.fixture
def functional_modules2():
    return random_functional_modules(random_seed=19, max_capacity=2, max_module_number_per_type=2)


class TestSchedulerInput:
    def test_build_input(self, operation_network1, functional_modules1, operation_network2, functional_modules2):
        for fms in [functional_modules1, functional_modules2]:
            for on in [operation_network1, operation_network2]:
                si = SchedulerInput.build_input(on, fms)
                # print(si.summary)
                assert si.summary['# of operations'] in [13, 41]
                assert si.summary['# of precedence relations'] in [8, 28]


class TestFormulationMILP:

    def test_solver1(self, operation_network1, functional_modules1):
        required_operation_types = set([o.type for o in operation_network1.operations])
        functional_modules = [fm for fm in functional_modules1 if
                              set(fm.can_process).issubset(required_operation_types)]
        si = SchedulerInput.build_input(operation_network1, functional_modules)
        solver = SolverMILP(input=si)
        solver.solve()

    def test_solver2(self, operation_network2, functional_modules2):
        required_operation_types = set([o.type for o in operation_network2.operations])
        functional_modules = [fm for fm in functional_modules2 if
                              set(fm.can_process).issubset(required_operation_types)]
        si = SchedulerInput.build_input(operation_network2, functional_modules)
        solver = SolverMILP(input=si)
        solver.solve()
