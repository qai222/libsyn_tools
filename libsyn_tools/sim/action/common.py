from .core import Action

"""
stub for now
"""


class Weigh(Action):
    container_iri: str
    balance_iri: str

    def get_resources(self):
        return [self.container_iri, self.balance_iri]

    def get_action_effects(self):
        # No graph mutation – measurement logged via event record
        return []


class Stir(Action):
    container_iri: str
    stirrer_iri: str
    duration_s: float | None = None

    def get_resources(self):
        return [self.container_iri, self.stirrer_iri]

    def get_action_effects(self):
        return []


class Shake(Action):
    container_iri: str
    shaker_iri: str
    duration_s: float | None = None

    def get_resources(self):
        return [self.container_iri, self.shaker_iri]

    def get_action_effects(self):
        return []


class Filter(Action):
    mixture_iri: str
    filter_media_iri: str
    filtrate_container_iri: str

    def get_resources(self):
        return [self.mixture_iri, self.filter_media_iri, self.filtrate_container_iri]

    def get_action_effects(self):
        # Stub – implement solid/liquid separation edits
        return []


class Dry(Action):
    container_iri: str
    dryer_iri: str

    def get_resources(self):
        return [self.container_iri, self.dryer_iri]

    def get_action_effects(self):
        return []
