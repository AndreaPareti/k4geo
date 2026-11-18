"""Validate the TB26 reference map, physical neighbours and unchanged fibre hierarchy."""

import math
from pathlib import Path
import sys

import DDG4
import ROOT


# Independent reference from the block explicitly labelled TB26 in HidraSim.
REFERENCE = (
    (-1, -1, 0, -1, -1),
    (-1, 1, 2, 3, -1),
    (-1, 4, 5, 6, -1),
    (-1, 7, 8, 9, -1),
    (10, 11, 12, 13, 14),
    (15, 16, 17, 18, 19),
    (20, 21, 22, 23, 24),
    (25, 26, 27, 28, 29),
    (30, 31, 32, 33, 34),
    (35, 36, 37, 38, 39),
    (40, 41, 42, 43, 44),
    (45, 46, 47, 48, 49),
    (50, 51, 52, 53, 54),
    (55, 56, 57, 58, 59),
    (60, 61, 62, 63, 64),
    (65, 66, 67, 68, 69),
    (-1, 70, 71, 72, -1),
    (-1, 73, 74, 75, -1),
    (-1, 76, 77, 78, -1),
    (-1, -1, 79, -1, -1),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def save_map(path, positions, outline, mm):
    """Write a front-view SVG using actual placed polygon vertices and module IDs."""
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="800" '
        'viewBox="-335 -292 670 584">',
        '<rect x="-335" y="-292" width="670" height="584" fill="white"/>',
    ]
    for module_id, (x, y) in sorted(positions.items()):
        points = " ".join(f"{(x + vx) / mm:.5f},{-(y + vy) / mm:.5f}" for vx, vy in outline)
        fill = "#e6f1ff" if module_id % 2 else "#eef5e8"
        parts.append(
            f'<polygon points="{points}" fill="{fill}" stroke="#35516b" stroke-width="0.3"/>'
        )
        parts.append(
            f'<text x="{x / mm:.5f}" y="{-y / mm:.5f}" text-anchor="middle" '
            f'dominant-baseline="middle" font-family="sans-serif" font-size="9">{module_id}</text>'
        )
    parts.append("</svg>")
    Path(path).write_text("\n".join(parts) + "\n")


def main():
    ROOT.gROOT.SetBatch(True)
    kernel = DDG4.Kernel()
    kernel.loadGeometry("file:" + str(Path(sys.argv[1]).resolve()))
    detector = kernel.detectorDescription()
    manager = detector.manager()
    mm = ROOT.dd4hep.mm
    pitch_x, pitch_y = 128 * mm, 16 * math.sqrt(3) * mm
    expected = {
        module: ((2 - column) * pitch_x, (row - 9.5) * pitch_y)
        for row, entries in enumerate(REFERENCE)
        for column, module in enumerate(entries)
        if module >= 0
    }
    envelope = manager.GetVolume("HidraPrototype_envelope")
    require(envelope.GetNdaughters() == 1, "Expected one air stave")
    stave = envelope.GetNode(0).GetVolume()
    require(stave.GetMaterial().GetName() == "Air", "Stave introduces support material")
    require(stave.GetNdaughters() == 80, "TB26 must contain exactly 80 modules")
    positions = {}
    module_air = manager.GetVolume("HidraPrototype_module_air")
    require(module_air.GetNdaughters() == 1024, "Wrong fibre count per module")
    for index in range(stave.GetNdaughters()):
        node = stave.GetNode(index)
        module_id = node.GetNumber()
        require(
            module_id in expected and module_id not in positions, "Invalid or duplicate module ID"
        )
        translation = node.GetMatrix().GetTranslation()
        x, y, z = translation[0], translation[1], translation[2]
        positions[module_id] = (x, y)
        require(
            math.isclose(x, expected[module_id][0], abs_tol=1e-9 * mm)
            and math.isclose(y, expected[module_id][1], abs_tol=1e-9 * mm)
            and math.isclose(z, 0, abs_tol=1e-9 * mm),
            "TB26 placement/orientation mismatch",
        )
        rotation = node.GetMatrix().GetRotationMatrix()
        require(
            all(
                math.isclose(rotation[i], 1 if i in (0, 4, 8) else 0, abs_tol=1e-12)
                for i in range(9)
            ),
            "TB26 must be unrotated",
        )
        tower = node.GetVolume()
        require(
            tower.GetMaterial().GetName() == "Air", "A module border contains support material"
        )
        require(tower.GetNode(0).GetNumber() == 63, "IDEA air copy-number depth changed")
        require(
            ROOT.addressof(tower.GetNode(0).GetVolume()) == ROOT.addressof(module_air),
            "Modules do not reuse the same validated tube layout",
        )

        # Navigate to an S and C core in every module and check the real ancestor
        # depths that the unchanged IDEA sensitive action relies on.
        for row in (0, 1):
            tube_node = module_air.GetNode(row)
            local = tube_node.GetMatrix().GetTranslation()
            current = manager.FindNode(x + local[0], y + local[1], 100 * mm)
            require(current.GetNumber() == (2 | row), "Wrong S/C core copy number")
            nav = manager.GetCurrentNavigator()
            require(
                nav.GetMother(1).GetNumber() == 1
                and nav.GetMother(2).GetNumber() == row
                and nav.GetMother(3).GetNumber() == 63
                and nav.GetMother(4).GetNumber() == module_id
                and nav.GetMother(5).GetNumber() == 0,
                "IDEA hierarchy/module addressing changed",
            )
    require(set(positions) == set(range(80)), "Module IDs are not exactly 0 through 79")

    # Check steel-to-steel contacts across real X and Y module boundaries,
    # sampling a small distance inside each tube on either side of the contact.
    radius = mm
    half_height = ((16 - 1) * math.sqrt(3) + 4 / math.sqrt(3)) * mm / 2
    epsilon = 1e-3 * mm

    def steel_module(x, y):
        node = manager.FindNode(x, y, 100 * mm)
        require(
            node.GetVolume().GetMaterial().GetName() == "HidraSteel304",
            "A module contact has an air gap or added border material",
        )
        return manager.GetCurrentNavigator().GetMother(3).GetNumber()

    x, y = positions[37]
    contact_x = x + 64.5 * mm
    contact_y = y - half_height + 2 * radius / math.sqrt(3)
    require(
        steel_module(contact_x - epsilon, contact_y) == 37
        and steel_module(contact_x + epsilon, contact_y) == 36,
        "Horizontal module contact is not shared by the expected modules",
    )
    lower = (x + 62.5 * mm, y - half_height + 15 * math.sqrt(3) * mm + 2 * mm / math.sqrt(3))
    ux, uy = positions[42]
    upper = (ux + 63.5 * mm, uy - half_height + 2 * mm / math.sqrt(3))
    contact_x, contact_y = (lower[0] + upper[0]) / 2, (lower[1] + upper[1]) / 2
    require(
        steel_module(contact_x - epsilon / 2, contact_y - epsilon * math.sqrt(3) / 2) == 37
        and steel_module(contact_x + epsilon / 2, contact_y + epsilon * math.sqrt(3) / 2) == 42,
        "Vertical module contact is not shared by the expected modules",
    )

    manager.CheckOverlaps(1e-5 * mm)
    require(manager.GetListOfOverlaps().GetEntries() == 0, "TB26 has overlaps or extrusions")
    if len(sys.argv) > 2:
        shape = module_air.GetShape()
        outline = [(shape.GetX(i), shape.GetY(i)) for i in range(shape.GetNvert())]
        save_map(sys.argv[2], positions, outline, mm)
    print(
        "TB26 geometry OK: 80 unique modules, 81,920 fibres, exact reference placements, "
        "unchanged scoring depths, matching steel contacts and no overlaps."
    )
    kernel.terminate()


if __name__ == "__main__":
    main()
