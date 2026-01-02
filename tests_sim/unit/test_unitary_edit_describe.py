# ### THIS IS THE START OF CONTENT OF tests_sim/unit/test_unitary_edit_describe.py ###
from __future__ import annotations

from libsyn_tools.sim.operation.unitary_edit import (
    AddObjectProperty,
    ChangeDataProperty,
    Create,
)


def test_unitary_edit_describe_is_stable() -> None:
    create = Create(instance_1_iri="ex:thing")
    assert create.describe() == "CREATE ex:thing"

    set_prop = ChangeDataProperty(instance_1_iri="ex:thing", property_iri="ex:prop", data_value=1)
    assert set_prop.describe() == "SET ex:thing ex:prop = 1"

    add_rel = AddObjectProperty(
        instance_1_iri="ex:thing",
        instance_2_iri="ex:other",
        property_iri="ex:rel",
    )
    assert add_rel.describe() == "ADD ex:thing ex:rel ex:other"
# ### THIS IS THE END OF CONTENT OF tests_sim/unit/test_unitary_edit_describe.py ###
