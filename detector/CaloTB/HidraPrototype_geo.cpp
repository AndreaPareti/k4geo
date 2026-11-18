// HiDRa test-beam geometry: one module or a mapped array of hexagonal modules.
// Packing and the zigzag module outline are adapted from DRCalo/HidraSim's
// DREMTubesDetectorConstruction::calcmod. Nested tube construction follows
// k4geo's DRBarrelTubes::DRTubesconstructor::assert_tube_existence.

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Printout.h"
#include "XML/Utilities.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Outline {
  std::vector<double> x;
  std::vector<double> y;
};

// Trace the exposed edges of the hexagonal cells surrounding the outer tubes.
// This is an air envelope, not a support wall. Even row counts allow copies to
// tile with pitches 2*columns*radius and sqrt(3)*rows*radius, without a separator.
Outline moduleOutline(double radius, int rows, int columns) {
  const double sqrt3 = std::sqrt(3.);
  const double side = 2. * radius / sqrt3;
  const double width = (2. * columns + 1.) * radius;
  const double height = (rows - 1.) * sqrt3 * radius + 4. * radius / sqrt3;
  const int ny = 2 * rows + 1;
  const int nx = 2 * columns - 2;
  const int size = 2 * (ny + nx);
  Outline outline{std::vector<double>(size), std::vector<double>(size)};
  auto& x = outline.x;
  auto& y = outline.y;

  for (int i = 1; i < ny; ++i) {
    y[i] = y[i - 1] + ((i - 1) % 2 + 1) * 0.5 * side;
    x[i] = -((i + 1) / 2 % 2) * radius;
  }
  for (int i = ny; i < ny + nx; ++i) {
    const int j = i - ny;
    y[i] = y[ny - 1] + 0.5 * (i % 2) * side;
    x[i] = x[ny - 1] + radius * (j + 1);
  }
  for (int i = ny + nx; i < 2 * ny + nx; ++i) {
    const int j = i - ny - nx;
    y[i] = y[ny - j];
    x[i] = x[ny + nx - 1] + radius + ((j + 1) / 2 % 2) * radius;
  }
  for (int i = 2 * ny + nx; i < size; ++i) {
    const int j = i - 2 * ny - nx;
    y[i] = 0.5 * side * (i % 2);
    x[i] = x[ny + nx - j - 1];
  }
  for (int i = 0; i < size; ++i) {
    x[i] = -x[i] + width / 2. - radius;
    y[i] -= height / 2.;
  }
  return outline;
}

// The same two logical-volume trees are reused for all S/C tube placements.
dd4hep::Volume makeTube(dd4hep::Detector& description, dd4hep::SensitiveDetector sensitive, const std::string& name,
                        const xml_comp_t& module, bool cherenkov) {
  const xml_comp_t capillary_xml = module.child(_Unicode(capillary));
  const xml_comp_t cladding_xml = module.child(cherenkov ? _Unicode(cher_clad) : _Unicode(scin_clad));
  const xml_comp_t core_xml = module.child(cherenkov ? _Unicode(cher_core) : _Unicode(scin_core));
  const double half_length = module.attr<double>(_Unicode(length)) / 2.;
  const double radius = capillary_xml.outer_r();
  const double clad_radius = cladding_xml.outer_r();
  const double core_radius = core_xml.outer_r();
  if (!(half_length > 0. && core_radius > 0. && core_radius < clad_radius && clad_radius < radius)) {
    dd4hep::except(name.c_str(), "Require positive length and 0 < core radius < cladding radius < tube radius.");
  }

  // Partition the tube explicitly: steel annulus + cladding cylinder.
  // Their union fills the air container exactly, so this adds no air gap and
  // preserves the former steel-parent material map and the SD touchable depths.
  // A separate leaf shell is also drawn by ROOT's classic/web leaf renderers.
  dd4hep::Volume capillary(name, dd4hep::Tube(0., radius, half_length), description.air());
  dd4hep::Volume steel(name + "_steel", dd4hep::Tube(clad_radius, radius, half_length),
                       description.material(capillary_xml.materialStr()));
  dd4hep::Volume cladding(name + "_clad", dd4hep::Tube(0., clad_radius, half_length),
                          description.material(cladding_xml.materialStr()));
  dd4hep::Volume core(name + "_core", dd4hep::Tube(0., core_radius, half_length),
                      description.material(core_xml.materialStr()));
  capillary.setVisAttributes(description, "HidraAirVis");
  steel.setVisAttributes(description, capillary_xml.visStr());
  cladding.setVisAttributes(description, cladding_xml.visStr());
  core.setVisAttributes(description, core_xml.visStr());

  // Match IDEA: both cores score signal; C cladding receives the action's
  // optical-photon termination. S cladding is passive PMMA.
  core.setSensitiveDetector(sensitive);
  if (cherenkov)
    cladding.setSensitiveDetector(sensitive);

  auto clad_placement = capillary.placeVolume(cladding, 1);
  clad_placement.addPhysVolID("clad", 1).addPhysVolID("cherenkov", cherenkov ? 1 : 0);
  // DRTubesSDAction reads bits 0 (C/S) and 1 (core) of this copy number.
  auto core_placement = cladding.placeVolume(core, 2 | (cherenkov ? 1 : 0));
  core_placement.addPhysVolID("core", 1);
  // Keep cladding as child 0, then add the non-sensitive shell as child 1.
  capillary.placeVolume(steel, 0);
  return capillary;
}

struct ModulePlacement {
  int id;
  double x;
  double y;
};

// Row-major map as in standalone HidraSim: x decreases with grid column,
// y increases with grid row. Negative entries leave the corresponding slot empty.
std::vector<ModulePlacement> readLayout(const xml_det_t& detector_xml, double pitch_x, double pitch_y) {
  if (!detector_xml.hasChild(_Unicode(layout)))
    return {{0, 0., 0.}};

  const std::string name = detector_xml.nameStr();
  const xml_comp_t layout = detector_xml.child(_Unicode(layout));
  const int columns = layout.attr<int>(_Unicode(columns));
  const int rows = layout.attr<int>(_Unicode(rows));
  const int expected = layout.attr<int>(_Unicode(expected_modules));
  if (columns < 1 || rows < 1 || expected < 1 || expected > 128)
    dd4hep::except(name.c_str(), "Invalid layout dimensions or module count (tower IDs use signed 8-bit fields).");

  std::vector<ModulePlacement> placements;
  std::set<int> ids;
  int grid_row = 0;
  for (xml_coll_t row_xml(layout, _Unicode(row)); row_xml; ++row_xml, ++grid_row) {
    if (grid_row >= rows)
      dd4hep::except(name.c_str(), "Too many rows in the module layout.");
    const xml_comp_t row = row_xml;
    std::istringstream values(row.attr<std::string>(_Unicode(ids)));
    for (int column = 0; column < columns; ++column) {
      int id = -1;
      if (!(values >> id) || id < -1 || id > 127)
        dd4hep::except(name.c_str(), "Invalid or missing module ID at layout row %d column %d.", grid_row, column);
      if (id < 0)
        continue;
      if (!ids.insert(id).second)
        dd4hep::except(name.c_str(), "Duplicate module ID %d in layout.", id);
      placements.push_back({id, ((columns - 1.) / 2. - column) * pitch_x, (grid_row - (rows - 1.) / 2.) * pitch_y});
    }
    std::string extra;
    if (values >> extra)
      dd4hep::except(name.c_str(), "Too many entries in layout row %d.", grid_row);
  }
  if (grid_row != rows || static_cast<int>(placements.size()) != expected)
    dd4hep::except(name.c_str(), "Layout row/module counts do not match the declared values.");
  return placements;
}

dd4hep::Ref_t create_detector(dd4hep::Detector& description, xml_h element, dd4hep::SensitiveDetector sensitive) {
  const xml_det_t detector_xml = element;
  const std::string name = detector_xml.nameStr();
  dd4hep::DetElement detector(name, detector_xml.id());
  dd4hep::xml::setDetectorTypeFlag(element, detector);
  sensitive.setType("calorimeter");

  const xml_comp_t module_xml = detector_xml.child(_Unicode(module));
  const int rows = module_xml.attr<int>(_Unicode(rows));
  const int columns = module_xml.attr<int>(_Unicode(columns));
  if (rows < 2 || rows % 2 != 0 || rows > 65536 || columns < 1 || columns > 32768)
    dd4hep::except(name.c_str(), "Require an even row count in [2,65536] and columns in [1,32768].");
  const xml_comp_t capillary_xml = module_xml.child(_Unicode(capillary));
  const double radius = capillary_xml.outer_r();
  const double length = module_xml.attr<double>(_Unicode(length));
  const double width = (2. * columns + 1.) * radius;
  const double height = (rows - 1.) * std::sqrt(3.) * radius + 4. * radius / std::sqrt(3.);
  const auto placements = readLayout(detector_xml, 2. * columns * radius, std::sqrt(3.) * rows * radius);
  double half_x = 0., half_y = 0.;
  for (const auto& module : placements) {
    half_x = std::max(half_x, std::abs(module.x) + width / 2.);
    half_y = std::max(half_y, std::abs(module.y) + height / 2.);
  }
  const xml_dim_t dimensions = detector_xml.dimensions();
  if (dimensions.x() < 2. * half_x || dimensions.y() < 2. * half_y || dimensions.z() < length)
    dd4hep::except(name.c_str(), "The prototype envelope must contain the complete module layout.");

  dd4hep::Volume envelope(name + "_envelope",
                          dd4hep::Box(dimensions.x() / 2., dimensions.y() / 2., dimensions.z() / 2.),
                          description.air());
  envelope.setVisAttributes(description, detector_xml.visStr());
  const auto outline = moduleOutline(radius, rows, columns);
  dd4hep::ExtrudedPolygon module_shape(name + "_module_shape", outline.x, outline.y, {-length / 2., length / 2.},
                                       {0., 0.}, {0., 0.}, {1., 1.});

  // Preserve real volumes at SD depths 3/4/5. The full array uses one common
  // air stave containing interlocking towers; no rectangular box is placed
  // around an individual module and no support material is introduced.
  dd4hep::Solid stave_shape = module_shape;
  if (placements.size() > 1)
    stave_shape = dd4hep::Box(half_x, half_y, length / 2.);
  dd4hep::Volume stave_volume(name + "_stave", stave_shape, description.air());
  dd4hep::Volume tower_volume(name + "_tower", module_shape, description.air());
  dd4hep::Volume module_air(name + "_module_air", module_shape, description.air());
  for (auto volume : {stave_volume, tower_volume, module_air})
    volume.setVisAttributes(description, "HidraAirVis");

  const auto scintillation_tube = makeTube(description, sensitive, name + "_S_tube", module_xml, false);
  const auto cherenkov_tube = makeTube(description, sensitive, name + "_C_tube", module_xml, true);
  for (int column = 0; column < columns; ++column) {
    for (int row = 0; row < rows; ++row) {
      const bool cherenkov = row % 2 != 0;
      const double x = width / 2. - radius - 2. * radius * column - (cherenkov ? radius : 0.);
      const double y = -height / 2. + std::sqrt(3.) * radius * row + 2. * radius / std::sqrt(3.);
      auto placement = module_air.placeVolume(cherenkov ? cherenkov_tube : scintillation_tube, (column << 16) | row,
                                              dd4hep::Position(x, y, 0.));
      placement.addPhysVolID("col", column).addPhysVolID("row", row);
    }
  }

  auto air_placement = tower_volume.placeVolume(module_air, 63);
  air_placement.addPhysVolID("air", 63);
  auto stave_placement = envelope.placeVolume(stave_volume, 0);
  stave_placement.addPhysVolID("stave", 0);
  dd4hep::DetElement stave(detector, "stave_0", 0);
  stave.setPlacement(stave_placement);
  for (const auto& position : placements) {
    auto tower_placement =
        stave_volume.placeVolume(tower_volume, position.id, dd4hep::Position(position.x, position.y, 0.));
    tower_placement.addPhysVolID("tower", position.id);
    dd4hep::DetElement module(stave, "module_" + std::to_string(position.id), position.id);
    module.setPlacement(tower_placement);
  }
  auto mother = description.pickMotherVolume(detector);
  auto placement = mother.placeVolume(envelope);
  placement.addPhysVolID("system", detector_xml.id());
  detector.setPlacement(placement);
  dd4hep::printout(dd4hep::INFO, name, "Built %zu wall-free modules, each %d S + %d C fibres, in one air stave.",
                   placements.size(), columns * rows / 2, columns * rows / 2);
  return detector;
}

} // namespace

DECLARE_DETELEMENT(HidraPrototype_v01, create_detector)
