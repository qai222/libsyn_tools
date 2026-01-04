import pytest


def test_apply_tx_rolls_back_on_unexpected_exception(monkeypatch):
    """If an edit.apply() throws unexpectedly, the engine must restore snapshots and abort."""
    import simpy

    from libsyn_tools.sim.effect_engine import EffectEngine, EngineMechanicalError
    from libsyn_tools.sim.knowledge_graph import MaterialContainer
    from libsyn_tools.sim.operation.unitary_edit import Create

    env = simpy.Environment()
    # Create a present=false object that exists in KG lookup
    c = MaterialContainer(identifier="c")
    c.is_present = {False}

    edit = Create(instance_1_iri=c.identifier)

    # Force an unexpected exception inside apply()
    def boom():
        raise RuntimeError("boom")

    monkeypatch.setattr(edit, "apply", boom)

    eng = EffectEngine(shapes_graph=None)

    with pytest.raises(EngineMechanicalError):
        eng.apply_tx([edit], env, operation_id="op")

    # Must have been rolled back
    assert c.is_present == {False}
