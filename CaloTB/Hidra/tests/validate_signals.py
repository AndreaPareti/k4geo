"""Check the one-electron run exercises both channels and preserves physical CellIDs."""

import math
import sys

from podio import reading


def main():
    module_count = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    reader = reading.get_reader(sys.argv[1])
    events = list(reader.get("events"))
    if len(events) != 1:
        raise RuntimeError(f"Expected one event, found {len(events)}")
    # Bit positions are the readout contract, checked independently of the SD decoder.
    for name, kind in (("DRBTScin", 0), ("DRBTCher", 1)):
        hits = events[0].get(name)
        if len(hits) == 0:
            raise RuntimeError(f"No {name} hits: the signal action/optical setup is not exercised")
        total = 0.0
        seen = set()
        hit_modules = set()
        for hit in hits:
            cell = hit.getCellID()
            module_id = (cell >> 15) & 255
            hit_modules.add(module_id)
            column = (cell >> 29) & 0xFFFF
            row = (cell >> 45) & 0xFFFF
            if not (
                cell & 31 == 1
                and (cell >> 5) & 1023 == 0
                and 0 <= module_id < module_count
                and (cell >> 23) & 63 == 63
                and 0 <= column < 64
                and 0 <= row < 16
                and row % 2 == kind
                and (cell >> 61) & 3 == 3
                and (cell >> 63) & 1 == kind
            ):
                raise RuntimeError(f"Unexpected {name} CellID: {cell:#018x}")
            if cell in seen:
                raise RuntimeError("Duplicate fibre hit")
            seen.add(cell)
            signal = hit.getEnergy()  # Existing IDEA convention: this field contains PE counts.
            if not math.isfinite(signal) or signal <= 0:
                raise RuntimeError("Invalid IDEA PE signal")
            contributions = list(hit.getContributions())
            if not contributions:
                raise RuntimeError("Missing arrival-time contributions")
            if not math.isclose(sum(c.getEnergy() for c in contributions), signal, rel_tol=1e-5):
                raise RuntimeError("Arrival-bin sum does not match fibre signal")
            for contribution in contributions:
                time = contribution.getTime()
                if (
                    not math.isfinite(time)
                    or time < 0
                    or not math.isclose(time / 0.1, round(time / 0.1), abs_tol=1e-3)
                ):
                    raise RuntimeError("Invalid IDEA arrival-time bin")
            total += signal
        if module_count > 1 and len(hit_modules) < 2:
            raise RuntimeError(
                f"The TB26 electron run did not exercise multiple modules in {name}"
            )
        print(
            f"{name}: {len(hits)} fibres in {len(hit_modules)} modules, {total:g} PE; IDs and bins OK."
        )


if __name__ == "__main__":
    main()
