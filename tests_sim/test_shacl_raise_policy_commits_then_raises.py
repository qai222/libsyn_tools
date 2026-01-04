import pytest


def _overflow_shape(target_iri: str):
    """A minimal node shape: currentVolume must be <= 9.9."""
    from rdflib import Graph, Namespace, URIRef, BNode, Literal
    from rdflib.namespace import SH, XSD

    LIB = Namespace("https://libsyn-sim/kg/")
    g = Graph()
    shape = URIRef("urn:libsyn:shape:overflow")
    p = BNode()
    g.add((shape, SH.targetNode, URIRef(target_iri)))
    g.add((shape, SH.property, p))
    g.add((p, SH.path, LIB.currentVolume))
    g.add((p, SH.maxInclusive, Literal(9.9, datatype=XSD.double)))
    return g


def _world():
    """Minimal world: reservoir has 12 mL; destination vial starts empty."""
    from libsyn_tools.chem_schema import Chemical
    from libsyn_tools.sim.knowledge_graph import MaterialContainer, PortionOfMaterial

    res = MaterialContainer(identifier="reservoir")
    v1 = MaterialContainer(identifier="v1")
    pip = MaterialContainer(identifier="pipette")

    for o in (res, v1, pip):
        o.is_present = {True}
        o.has_pool_type = {"TEST"}

    # reservoir contains one POM of 12 mL
    pom = PortionOfMaterial()
    chem = Chemical(smiles="O", density=1.0, molecular_weight=18.0, mass=12.0)
    pom.add_chemical(chem)
    pom.is_present = {True}
    pom.is_directly_contained_by = {res}

    # capacity set just so we have a meaningful reference (not used by the shape)
    v1.has_capacity = {9.9}
    return v1, res, pip


def test_shacl_raise_policy_commits_then_raises():
    """With raise_shacl=True, SHACL violations should raise after commit."""
    from libsyn_tools.sim.simulation import Simulation
    from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByPortionSize
    from libsyn_tools.sim.effect_engine import SHACLValidationError

    v1, res, pip = _world()

    op = TransferMaterialByPortionSize(
        identifier="fill",
        participant_source=res.identifier,
        participant_destination=v1.identifier,
        participant_device=pip.identifier,
        portion_size=1.0,
        temporal_cost=0.0,
    )

    sim = Simulation([op], shacl_shapes=_overflow_shape(v1.identifier))
    sim.effect_engine.raise_shacl = True

    with pytest.raises(SHACLValidationError):
        sim.run()
