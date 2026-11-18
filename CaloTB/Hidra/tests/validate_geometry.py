"""Validate the built module, hex-cell tiling and IDEA hierarchy using ROOT geometry."""

import json
import math
import sys
from array import array

import DDG4
import ROOT


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    ROOT.gROOT.SetBatch(True)
    kernel = DDG4.Kernel()
    kernel.loadGeometry("file:" + sys.argv[1])
    detector = kernel.detectorDescription()
    manager = detector.manager()
    radius = detector.constantAsDouble("HidraTubeOuterRadius")
    length = detector.constantAsDouble("HidraModuleLength")
    columns = detector.constantAsLong("HidraFiberColumns")
    rows = detector.constantAsLong("HidraFiberRows")
    pitch_x = detector.constantAsDouble("HidraModulePitchX")
    pitch_y = detector.constantAsDouble("HidraModulePitchY")
    module = manager.GetVolume("HidraPrototype_module_air")
    shape = module.GetShape()
    require(shape.IsA().GetName() == "TGeoXtru", "Module must have the hex-cell extruded outline")
    require(module.GetNdaughters() == rows * columns, "Wrong tube count")

    envelope = manager.GetVolume("HidraPrototype_envelope")
    stave_node = envelope.GetNode(0)
    tower_node = stave_node.GetVolume().GetNode(0)
    air_node = tower_node.GetVolume().GetNode(0)
    require(
        (stave_node.GetNumber(), tower_node.GetNumber(), air_node.GetNumber()) == (0, 0, 63),
        "Wrong stave/tower/air copy numbers",
    )
    for volume in (stave_node.GetVolume(), tower_node.GetVolume(), module):
        require(
            volume.GetMaterial().GetName() == "Air",
            "A module separator/support material was introduced",
        )
        require(
            ROOT.addressof(volume.GetShape()) == ROOT.addressof(shape),
            "Stave/tower envelopes must share the interlocking outline",
        )

    centres = {}
    types = {0: 0, 1: 0}
    for index in range(module.GetNdaughters()):
        tube_node = module.GetNode(index)
        copy = tube_node.GetNumber()
        column, row = copy >> 16, copy & 0xFFFF
        require((column, row) not in centres, "Duplicate fibre ID")
        require(0 <= column < columns and 0 <= row < rows, "Tube address is out of range")
        centre = tube_node.GetMatrix().GetTranslation()
        centres[column, row] = (centre[0], centre[1])
        tube = tube_node.GetVolume()
        require(tube.GetMaterial().GetName() == "Air", "Tube navigation container must be air")
        require(not tube.IsVisible(), "Navigation tube container must stay invisible")
        require(
            tube.GetNdaughters() == 2, "Tube must contain cladding and an explicit steel annulus"
        )
        require(math.isclose(tube.GetShape().GetRmax(), radius), "Tube radius mismatch")
        require(math.isclose(tube.GetShape().GetDz(), length / 2), "Tube length mismatch")
        clad_node = tube.GetNode(0)
        shell = tube.GetNode(1).GetVolume()
        require(
            shell.GetMaterial().GetName() == "HidraSteel304", "Capillary shell material mismatch"
        )
        require(shell.GetNdaughters() == 0, "Steel shell must be a drawable leaf volume")
        require(shell.IsVisible(), "Steel shell must be visible")
        require(
            math.isclose(shell.GetShape().GetRmin(), clad_node.GetVolume().GetShape().GetRmax()),
            "Shell/cladding boundary has a gap or overlap",
        )
        require(math.isclose(shell.GetShape().GetRmax(), radius), "Steel outer radius mismatch")
        require(math.isclose(shell.GetShape().GetDz(), length / 2), "Steel length mismatch")
        require(
            math.isclose(
                shell.GetShape().Capacity() + clad_node.GetVolume().GetShape().Capacity(),
                tube.GetShape().Capacity(),
                rel_tol=1e-12,
            ),
            "Steel and cladding do not completely fill the former tube material map",
        )
        core_node = clad_node.GetVolume().GetNode(0)
        kind = row % 2
        require(
            clad_node.GetNumber() == 1 and core_node.GetNumber() == (2 | kind),
            "Tube hierarchy cannot be decoded by IDEA's unchanged action",
        )
        expected = "HidraPolystyrene" if kind == 0 else "HidraPMMA"
        require(
            core_node.GetVolume().GetMaterial().GetName() == expected,
            "S/C material or row parity mismatch",
        )
        types[kind] += 1
        # Check the tube circumference against the actual module shape, including edge tubes.
        for angle_index in range(24):
            angle = angle_index * math.tau / 24
            inset_radius = radius * (1 - 1e-8)
            point = array(
                "d",
                (
                    centre[0] + inset_radius * math.cos(angle),
                    centre[1] + inset_radius * math.sin(angle),
                    0,
                ),
            )
            require(shape.Contains(point), f"Tube ({column}, {row}) protrudes outside the module")
    require(types[0] == types[1] == columns * rows // 2, "Unbalanced S/C tube population")

    # Check navigated materials inside each radial layer, away from interfaces.
    for row, core_material, clad_material in (
        (0, "HidraPolystyrene", "HidraPMMAScin"),
        (1, "HidraPMMA", "HidraFluorinatedPolymer"),
    ):
        x, y = centres[0, row]
        core_radius = detector.constantAsDouble("HidraCoreOuterRadius")
        clad_radius = detector.constantAsDouble("HidraCladdingOuterRadius")
        for radial, expected in (
            (core_radius / 2, core_material),
            ((core_radius + clad_radius) / 2, clad_material),
            ((clad_radius + radius) / 2, "HidraSteel304"),
        ):
            node = manager.FindNode(x + radial, y, length / 4)
            require(
                node.GetVolume().GetMaterial().GetName() == expected,
                f"Radial material map changed in row {row}: expected {expected}",
            )

    # Every neighbouring row pair must touch at the common capillary radius.
    for (column, row), (x, y) in centres.items():
        if row + 1 < rows:
            nx, ny = centres[column, row + 1]
            require(
                math.isclose(math.hypot(nx - x, ny - y), 2 * radius, rel_tol=1e-11),
                "Adjacent-row tubes have a gap/overlap",
            )

    points = [(shape.GetX(i), shape.GetY(i)) for i in range(shape.GetNvert())]
    area = abs(sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(points, points[1:] + points[:1]))) / 2
    require(
        math.isclose(area, pitch_x * pitch_y, rel_tol=1e-11),
        "Module footprint does not fill one repeat cell",
    )

    def edge_set(dx=0, dy=0):
        # Normalize to radius to make tolerance independent of ROOT length units.
        return {
            tuple(
                sorted(
                    (
                        (round((a[0] + dx) / radius, 8), round((a[1] + dy) / radius, 8)),
                        (round((b[0] + dx) / radius, 8), round((b[1] + dy) / radius, 8)),
                    )
                )
            )
            for a, b in zip(points, points[1:] + points[:1])
        }

    edges = edge_set()
    require(
        len(edges & edge_set(pitch_x, 0)) == 2 * rows - 1,
        "Horizontal module neighbours do not match along their zigzag edges",
    )
    require(
        len(edges & edge_set(0, pitch_y)) == 2 * columns - 1,
        "Vertical module neighbours do not match along their zigzag edges",
    )
    # Validate actual cross-boundary tube spacing as well as envelope boundaries.
    for row in range(rows):
        left, right = centres[0, row], centres[columns - 1, row]
        require(
            math.isclose(right[0] + pitch_x - left[0], 2 * radius, rel_tol=1e-11),
            "Cross-module tubes do not touch in X",
        )
    bottom = centres[0, 0]
    top = centres[0, rows - 1]
    require(
        math.isclose(
            math.hypot(bottom[0] - top[0], bottom[1] + pitch_y - top[1]), 2 * radius, rel_tol=1e-11
        ),
        "Cross-module tubes do not touch in Y",
    )

    # Regress the blank display caused by stopping at hidden cladding layers.
    manager.SetVisLevel(8)
    manager.SetVisOption(0)
    canvas = ROOT.TCanvas("HidraVisibilityCheck", "HiDRa visibility", 800, 600)
    manager.GetTopVolume().Draw("")
    canvas.Update()
    visible = manager.GetGeomPainter().CountVisibleNodes()
    require(visible >= rows * columns, "Classic viewer produced incomplete or empty drawing data")
    # Check the actual volumes rendered, not only a nonzero total node count.
    for name in (
        "HidraPrototype_S_tube_steel",
        "HidraPrototype_C_tube_steel",
        "HidraPrototype_S_tube_core",
        "HidraPrototype_C_tube_core",
    ):
        volume = manager.GetVolume(name)
        require(
            volume.TestAttBit(ROOT.TGeoAtt.kVisOnScreen), f"Classic combined view omitted {name}"
        )
    canvas.Close()
    viewer = ROOT.RGeomViewer(manager)
    viewer.SetLimits()  # Match geoWebDisplay's actual default node/face limits.
    web = viewer.Description()
    web.SetVisLevel(8)
    web.ProduceDrawData()
    payload = str(web.GetDrawJson())
    require("{" in payload, "Web viewer produced no geometry drawing data")
    drawing = json.loads(payload[payload.index("{") :])
    # Inspect emitted scene objects: visibility flags alone missed this regression.
    nodes = {node["id"]: node for node in drawing["nodes"]}
    rendered = {"S_steel": 0, "C_steel": 0, "S_core": 0, "C_core": 0}
    for item in drawing["visibles"]:
        node = nodes[item["nodeid"]]
        name = node["name"]
        for kind in ("S", "C"):
            if f"HidraPrototype_{kind}_tube_steel" in name:
                require(
                    node["material"] == "HidraSteel304", "Rendered steel has the wrong material"
                )
                rendered[f"{kind}_steel"] += 1
            elif f"HidraPrototype_{kind}_tube_core" in name:
                rendered[f"{kind}_core"] += 1
    for kind, count in rendered.items():
        require(count == rows * columns // 2, f"Web renderer emitted {count} {kind} objects")
    print(f"Display scene OK: {rendered}; classic rendering includes shells and cores.")

    manager.CheckOverlaps(radius * 1e-5)
    require(
        manager.GetListOfOverlaps().GetEntries() == 0, "ROOT reports module overlaps/extrusions"
    )
    print(
        f"Geometry OK: {types[0]} S + {types[1]} C tubes; wall-free matching boundaries; no overlaps."
    )
    kernel.terminate()


if __name__ == "__main__":
    main()
