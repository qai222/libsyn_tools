Library synthesis tools
---

Software tools developed for scheduling and simulating synthesis campaigns.

For extending the simulator with custom Operations or spawners, see
`docs/sim/EXTENDING.md` and the runnable example `examples/sim_minimal.py`.


### Schedule optimization

#### basic usage
Scheduler input and output are defined in [libsyn_tools/opt/schema.py](libsyn_tools/opt/schema.py).
Detailed definitions of the variables/parameters can be found in the `Methods` section of 
the article `Schedule Optimization for Chemical Library Synthesis`. 
In addition to manual construction, scheduler input instances can be automatically constructed from given
[`OperationNetwork` instance](libsyn_tools.chem_schema.network_operation.OperationNetwork) instances.
```python
from libsyn_tools.opt import SchedulerInput
from libsyn_tools.utils import json_load
operation_network = json_load("operation_network.json")
scheduler_input = SchedulerInput.build_input(operation_network=operation_network, ...)
```

Two methods are provided to solve the scheduling problem, the baseline greedy algorithm and the MILP exact formulation
using `gurobi` (a free licence is available through the distributor for academic users).
```python
from libsyn_tools.opt import SolverMILP, SolverBaseline
solver_baseline = SolverBaseline(input=scheduler_input, ...)
solver_milp = SolverMILP(input=scheduler_input, ...)
solver_baseline.solve()
solver_milp.solve()
# outputs are available at solver_x.output
```

#### reproduce paper results
Chemical libraries are stored in folders in [workplace_opt/LIBS.zip](workplace_opt/LIBS.zip).
Each folder contains an `operation_network.json` and a `reaction_network.json` 
that can be loaded as objects using `libsyn_tools.utils.json_load`.
Given the folder path containing the aforementioned two JSON files, 
the script [workplace_opt/run_scheduler.py](workplace_opt/run_scheduler.py) automatically constructs scheduling instance
and solve them using the baseline and the MILP formulation.
- Sample output folders can be found in [lst_app/routes](lst_app/routes). 
This path also allows registering scheduler output, 
so it can be visualized using the app at [lst_app/lst_app.py](lst_app/lst_app.py).
- Output folders of the 720 scheduling instances from the paper are stored in [workplace_opt/RUNS.7z](workplace_opt/RUNS.7z).






[//]: # (### Chemistry Platform Simulator)
