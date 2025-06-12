from twa.data_model.base_ontology import KnowledgeGraph

from libsyn_tools.chem_schema import ChemicalReaction
from libsyn_tools.sim.knowledge_graph import LabObject, PortionOfMaterial
from .core import Action, UnitaryEdit, UnitaryEditType


# TODO presumptions

class Evaporate(Action):
    """ heat container """

    container_iri: str
    """ the iri of the container """

    evaporation_module_iri: str
    """ the iri of the module """

    reaction: ChemicalReaction

    def get_resources(self) -> list[str]:
        """ a list of iris of the resources, they are assumed to be `LabObject` instances """
        resources = [self.container_iri, self.evaporation_module_iri]
        for iri in resources:
            resource = KnowledgeGraph.get_object_from_lookup(iri=iri)
            resources += [part.instance_iri for part in LabObject.get_all_parts(resource)]
        return resources

    def get_action_effects(self) -> list[UnitaryEdit]:
        unitary_edits = []

        container = KnowledgeGraph.get_object_from_lookup(self.container_iri)
        # heating_module = KnowledgeGraph.get_object_from_lookup(self.heating_module_iri)
        poms_before = LabObject.get_directly_contained_individuals(container, PortionOfMaterial)
        for pom_before in poms_before:
            annihilate_pom_before = UnitaryEdit(
                type=UnitaryEditType.ANNIHILATE,
                instance_1_iri=pom_before.instance_iri,
            )
            unitary_edits.append(annihilate_pom_before)
        pom_after = PortionOfMaterial()
        for c in self.reaction.products:
            pom_after.has_ingredient.add(c.model_dump_json())
        create_pom_after = UnitaryEdit(
            type=UnitaryEditType.CREATE,
            instance_1_iri=pom_after.instance_iri,
        )
        unitary_edits.append(create_pom_after)
        return unitary_edits
