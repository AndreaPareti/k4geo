"""HiDRa viewer with explicit limits for all 80 TB26 modules."""

import argparse
import json
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--web", action="store_true", help="Use the web viewer (default)")
    mode.add_argument("--classic", action="store_true", help="Use native geoDisplay")
    parser.add_argument(
        "--single-module", action="store_true", help="Load the original one-module geometry"
    )
    parser.add_argument("--xml", type=Path, help="Override the compact file")
    content = parser.add_mutually_exclusive_group()
    content.add_argument(
        "--modules",
        action="store_true",
        help="Draw solid modules in orange/blue checkerboard colours; hide tubes and fibres",
    )
    content.add_argument(
        "--details", action="store_true", help="Draw all tubes and fibre cores (default)"
    )
    content.add_argument(
        "--module-borders",
        action="store_true",
        help="Draw only wireframe module outlines; hide all tubes and fibres",
    )
    parser.add_argument("-level", "--level", type=int, default=8)
    parser.add_argument("-visopt", "--visopt", type=int, default=0)
    parser.add_argument("-maxnodes", "--maxnodes", type=int, default=200000)
    parser.add_argument("--maxfaces", type=int, default=1000000)
    parser.add_argument(
        "--check-scene",
        action="store_true",
        help="Validate emitted display objects without opening a window",
    )
    args, extra = parser.parse_known_args()
    name = "HidraPrototype_v01.xml" if args.single_module else "HidraTB26_v01.xml"
    compact = args.xml or Path(__file__).resolve().parents[1] / "compact" / name
    compact = compact.resolve()
    if args.classic:
        if args.check_scene:
            parser.error("--check-scene uses the web scene generator")
        if args.modules or args.module_borders:
            parser.error("Module display modes use the web viewer; omit --classic")
        os.execvp(
            "geoDisplay",
            [
                "geoDisplay",
                str(compact),
                "-level",
                str(args.level),
                "-visopt",
                str(args.visopt),
                *extra,
            ],
        )
    if extra:
        parser.error("Unrecognized web viewer arguments: " + " ".join(extra))

    import ROOT
    import DDG4

    if args.check_scene:
        ROOT.gROOT.SetBatch(True)
    kernel = DDG4.Kernel()
    kernel.loadGeometry("file:" + str(compact))
    manager = kernel.detectorDescription().manager()
    manager.SetVisLevel(args.level)
    if args.modules or args.module_borders:
        # Display-only changes in this process: use the actual module polygons.
        # Hide the logical parents and stop traversal at each tower placement.
        for volume_name in (
            "HidraPrototype_envelope",
            "HidraPrototype_stave",
            "HidraPrototype_module_air",
            "HidraPrototype_S_tube",
            "HidraPrototype_C_tube",
            "HidraPrototype_S_tube_clad",
            "HidraPrototype_C_tube_clad",
            "HidraPrototype_S_tube_core",
            "HidraPrototype_C_tube_core",
            "HidraPrototype_S_tube_steel",
            "HidraPrototype_C_tube_steel",
        ):
            manager.GetVolume(volume_name).SetVisibility(False)
        manager.GetTopVolume().SetVisibility(False)
        tower = manager.GetVolume("HidraPrototype_tower")
        tower.SetVisibility(True)
        tower.SetTransparency(0)
        tower.SetLineColor(ROOT.kBlue + 1)
        tower.SetLineWidth(2)
        stave = manager.GetVolume("HidraPrototype_stave")
        for index in range(stave.GetNdaughters()):
            node = stave.GetNode(index)
            node.SetVisibility(True)
            node.SetVisDaughters(False)
        if args.modules:
            # Placements share one logical tower. Make two display-only copies
            # with the same shape and contents so each placement can have a colour.
            display_volumes = []
            for colour, rgb in (("blue", "#0072B2"), ("orange", "#E69F00")):
                volume = tower.CloneVolume()
                volume.SetName(f"HidraPrototype_tower_display_{colour}")
                volume.SetLineColor(ROOT.TColor.GetColor(rgb))
                display_volumes.append(volume)
            positions = {
                str(stave.GetNode(i).GetName()): tuple(
                    float(stave.GetNode(i).GetMatrix().GetTranslation()[axis]) for axis in (0, 1)
                )
                for i in range(stave.GetNdaughters())
            }
            columns = {x: i for i, x in enumerate(sorted({p[0] for p in positions.values()}))}
            rows = {y: i for i, y in enumerate(sorted({p[1] for p in positions.values()}))}
            module_grid = {name: (columns[x], rows[y]) for name, (x, y) in positions.items()}
            for index in range(stave.GetNdaughters()):
                node = stave.GetNode(index)
                column, row = module_grid[str(node.GetName())]
                node.SetVolume(display_volumes[(column + row) % 2])
    viewer = ROOT.RGeomViewer(manager)
    viewer.SetTitle("HiDRa TB26" if not args.single_module else "HiDRa single module")
    viewer.SetLimits(args.maxnodes, args.maxfaces)
    if args.modules:
        viewer.SetDrawOptions("")
        viewer.SetTitle("HiDRa solid modules")
    elif args.module_borders:
        viewer.SetDrawOptions("wire")
        viewer.SetTitle("HiDRa module borders")
    if args.check_scene:
        scene = viewer.Description()
        scene.ProduceDrawData()
        payload = str(scene.GetDrawJson())
        if "{" not in payload:
            raise RuntimeError("Viewer produced no scene data")
        drawing = json.loads(payload[payload.index("{") :])
        nodes = {node["id"]: node for node in drawing["nodes"]}
        if args.modules or args.module_borders:
            stave = manager.GetVolume("HidraPrototype_stave")
            expected = {str(stave.GetNode(i).GetName()) for i in range(stave.GetNdaughters())}
            drawn = [nodes[item["nodeid"]]["name"] for item in drawing["visibles"]]
            if len(drawn) != len(expected) or set(drawn) != expected:
                raise RuntimeError(f"Module scene contains incorrect objects: {drawn}")
            expected_option = "wire" if args.module_borders else ""
            if drawing["cfg"].get("drawopt") != expected_option:
                raise RuntimeError("Incorrect solid/wireframe module rendering mode")
            if args.modules:
                colours = {
                    nodes[item["nodeid"]]["name"]: nodes[item["nodeid"]]["color"]
                    for item in drawing["visibles"]
                }
                palette = ("#0072b2", "#e69f00")
                expected_colours = {
                    name: palette[(column + row) % 2]
                    for name, (column, row) in module_grid.items()
                }
                if colours != expected_colours:
                    raise RuntimeError(
                        "Module scene does not have orange/blue checkerboard colours"
                    )
                print(
                    "Checkerboard colours OK: "
                    f"{sum(c == palette[0] for c in colours.values())} blue, "
                    f"{sum(c == palette[1] for c in colours.values())} orange."
                )
            style = "wireframe" if args.module_borders else "solid"
            print(f"Module scene OK: {len(drawn)} {style} modules, no tubes or fibres.")
            kernel.terminate()
            return
        counts = {"S_steel": 0, "C_steel": 0, "S_core": 0, "C_core": 0}
        for item in drawing["visibles"]:
            name = nodes[item["nodeid"]]["name"]
            for kind in ("S", "C"):
                if f"HidraPrototype_{kind}_tube_steel" in name:
                    counts[f"{kind}_steel"] += 1
                elif f"HidraPrototype_{kind}_tube_core" in name:
                    counts[f"{kind}_core"] += 1
        modules = manager.GetVolume("HidraPrototype_stave").GetNdaughters()
        fibres = manager.GetVolume("HidraPrototype_module_air").GetNdaughters()
        expected = modules * fibres // 2
        if any(value != expected for value in counts.values()):
            raise RuntimeError(
                f"Viewer limits/visibility truncate the geometry: {counts}, expected {expected} each"
            )
        print(f"Full scene OK: {counts}; no modules/shells/cores omitted.")
        kernel.terminate()
        return
    viewer.Show()
    print("Close the viewer and stop this process with Ctrl+C when finished.")
    ROOT.gApplication.Run()


if __name__ == "__main__":
    main()
