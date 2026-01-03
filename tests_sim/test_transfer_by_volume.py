from __future__ import annotations

import pytest

from libsyn_tools.sim.knowledge_graph import MaterialContainer
from libsyn_tools.sim.operation_preset.transfer import TransferMaterialByVolume


def _make_container(*, capacity: float = 100.0, present: bool = True) -> MaterialContainer:
    c = MaterialContainer()
    c.has_capacity.add(capacity)
    c.is_present = {present}
    return c


def test_transfer_by_volume_empty_source_raises_clear_error() -> None:
    src = _make_container(capacity=100.0)
    dst = _make_container(capacity=100.0)
    dev = _make_container(capacity=100.0)

    op = TransferMaterialByVolume(
        participant_source=src.identifier,
        participant_destination=dst.identifier,
        participant_device=dev.identifier,
        transfer_volume=1.0,
    )

    with pytest.raises(RuntimeError) as exc:
        op.get_operation_effects()

    # Avoid ZeroDivisionError / cryptic failures.
    msg = str(exc.value).lower()
    assert "no material" in msg or "empty" in msg or "no pom volume" in msg
