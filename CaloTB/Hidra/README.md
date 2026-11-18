# HiDRa test-beam geometries

Two compact descriptions use the same constructor, materials, fibre packing,
readout and unchanged IDEA tube sensitive action:

| Compact file | Geometry |
| --- | --- |
| `HidraTB26_v01.xml` | TB26: 80 modules, centred and unrotated; default for steering/display |
| `HidraPrototype_v01.xml` | Original single module, retained for regression and close inspection |

## TB26 layout

`compact/prototype_tb26.xml` contains the 5-column by 20-row slot map from the
block explicitly labelled **TB26** in standalone HidraSim's
`include/DREMTubesGeoPar.hh`. Empty slots are `-1`; the 80 active modules have
unique IDs 0 through 79. This reference has ID 11 in its proper slot, unlike the
duplicate ID in the later standalone active-map variant. The constructor rejects
duplicate IDs, malformed rows, mismatched module counts and out-of-range IDs.

Grid coordinates follow standalone HidraSim: x decreases as the grid column
increases, y increases with grid row, and all modules have z=0. For grid column
c and row r, the module centre is:

```text
x = (2 - c) * 128 mm
y = (r - 9.5) * (16 * sqrt(3)) mm
z = 0
```

No platform translation or rotation is applied. The array contains 81,920 fibres:
40,960 scintillation and 40,960 Cherenkov channels. It spans approximately
641 mm by 554.833 mm transversely and 2500 mm along the fibres. One common air
stave contains 80 placements of a shared tower/module volume. Each tower copy
number is the module ID; the stave copy number remains zero.

The full module map is the first TB26 step. Sensor packages, PMT/SiPM hardware
grouping, beamline/preshower/leakage devices and nonzero platform settings are
not yet included. The response remains IDEA's uniform per-fibre model.

## Module geometry and matching boundaries

Each module contains 64 columns by 16 staggered rows of capillaries: 512 S and
512 C fibres, each 2500 mm long. The standalone prototype supplies the radii
(tube 1 mm, cladding 0.5 mm, core 0.485 mm) and stainless steel (7.93 g/cm3;
Fe/Cr/Ni = 68/20/12%). The readout is at the positive local-z fibre end.

The module boundary follows the exposed edges of the outer hexagonal tube
cells, using HidraSim's `calcmod` construction. It is a zigzag polygon rather
than a rectangular support box. Tower and module-air envelopes share that
outline; the full array's common stave is an invisible Air box. There are no
module walls, coatings or inter-module clearances.

Packing uses exact `sqrt(3)`, replacing the standalone rounded value `1.733`.
The repeat pitches are 128 mm in x and approximately 27.712813 mm in y. The
module's bounding box is 129 mm by approximately 28.290163 mm: neighbouring
bounding boxes intersect while their actual zigzag polygons match exactly.
The tests navigate through touching steel on both sides of X and Y module
contacts, in addition to checking all module positions and ROOT overlaps.

## Explicit capillaries

Each tube has an invisible cylindrical Air container whose children are a
visible steel annulus and a cladding cylinder; the core remains inside the
cladding. The steel extends from the cladding outer radius (0.5 mm) to the tube
outer radius (1 mm). Steel and cladding completely fill the container, with no
added air gap and the same spatial material map as the original steel-parent
representation. The separate steel leaf makes capillaries and cores reliably
visible in both ROOT renderers. Steel is grey, S cores blue and C cores red.

## IDEA hierarchy and signal model

Actual volumes preserve the fibre-core touchable depths expected by the
unchanged IDEA barrel decoder:

| Depth from core | Volume | Copy number |
| --- | --- | --- |
| 0 | Core | 2 for S, 3 for C |
| 1 | Cladding | 1 |
| 2 | Invisible tube container | `(column << 16) | row` |
| 3 | Module air | 63 |
| 4 | Tower | TB26 module ID 0–79; 0 for the single module |
| 5 | Common air stave | 0 |

The envelope carries system ID 1. The `BARREL` flag selects IDEA's decoder;
the prototype itself has straight parallel fibres. The readout uses IDEA's
barrel field encoding with physical-fibre segmentation. The signed 8-bit tower
field accommodates all 80 IDs. Even rows are S and odd rows C; standalone
type-local fibre numbers are `column*(rows/2) + row/2`.

`plugins/DRTubesSDAction.cpp` and `plugins/DRTubesSglHpr.hh` are unchanged.
Steering attaches that action to the S/C cores and C cladding, with the same
Cherenkov/optical setup and no filter that rejects zero-deposit optical steps.
Scintillation conversion, Cherenkov conversion, Birks correction, attenuation,
timing and the 0.1 ns arrival-bin default remain IDEA's model. `DRBTScin` and
`DRBTCher` retain its PE-in-energy-field convention; `DRTubeEdep` is separate
truth information. No HidraSim raw-bin conversion or SimSiPM digitizer is added.

## Run

After setting up Key4hep and this k4geo build, run from a separate output directory:

```bash
ddsim --steeringFile /path/to/k4geo/CaloTB/Hidra/run_sim/steering.py
```

The default is one 10 GeV electron incident along positive z on TB26. To retain
the single-module comparison, override the compact file on the command line:

```bash
ddsim --steeringFile /path/to/k4geo/CaloTB/Hidra/run_sim/steering.py \
  --compactFile /path/to/k4geo/CaloTB/Hidra/compact/HidraPrototype_v01.xml \
  --outputFile HidraSingleModule.edm4hep.root
```

## Display

Stop an old display process before loading newly rebuilt geometry. The launcher
defaults to TB26 and the ROOT web viewer, with explicit limits of 200,000 visible
nodes and 1,000,000 faces. The complete array needs 163,840 shell/core objects;
native viewers' smaller defaults can otherwise silently truncate it.

```bash
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display.sh --modules
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display.sh --details
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display.sh
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display.sh --single-module
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display.sh --check-scene
```

`--modules` draws the full, opaque module envelopes at their actual placements,
hiding their contents. Blue and orange alternate by grid row and column in a
checkerboard pattern, including around empty slots. Use it to check the complete
calorimeter layout. Colours apply only in the viewer process.
`--details` draws all steel capillaries and S/C fibre cores; this remains the
default when no content option is given. These two options and `--module-borders`
are mutually exclusive. Add `--single-module` to any mode to inspect one module.
The module modes use the web viewer.

`--classic` invokes native `geoDisplay` with the XML filename first, depth 8 and
`-visopt 0`; its own rendering limits may need adjustment for the full array.
The supplied web launcher explicitly controls both node and face limits and
checks every emitted S/C steel/core object with `--check-scene`. `--level`,
`--maxnodes`, `--maxfaces` and `--xml` allow overrides.

For module borders alone, use the dedicated web launcher:

```bash
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display_modules.sh
bash /path/to/k4geo/CaloTB/Hidra/run_sim/display_modules.sh --single-module
```

This draws the actual interlocking module envelopes in blue wireframe and hides
all tubes and fibres. It is equivalent to `display.sh --module-borders`.
Visibility changes apply only inside the viewer process; no material separator
is added to the simulation. Shift + left drag or right drag pans the view,
left drag rotates it, and the mouse wheel zooms.

## Validation and layout map

CTest label `Hidra` runs eleven tests: the existing three single-module checks
and four TB26 checks (geometry/map, electron simulation, multi-module signals,
and complete display scene), plus solid-module and module-border scenes for both layouts.
They verify:

- The reference positions and all unique module IDs 0–79, with no rotations.
- 1024 fibres per module and the unchanged physical scoring hierarchy.
- Explicit steel-shell/core materials, matching boundaries, and zero overlaps.
- Positive S/C signals spanning multiple modules, valid CellIDs and arrival bins.
- All 40,960 objects of each S-steel, C-steel, S-core and C-core type in the full scene.
- Exactly 80 module wireframes (one for the single module) with no tube/fibre meshes.
- Exactly 80 solid module envelopes (one for the single module) with no tube/fibre meshes.

The TB26 geometry test writes `HidraTB26_module_map.svg` in the test build
directory. Its labelled polygons come from actual module placements and are
viewed along z with positive y upwards. These are geometry and integration
checks, not a validation of test-beam calibration or detector performance.
