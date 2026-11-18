"""TB26 HiDRa geometry using IDEA option 2's unchanged tube sensitive action."""

from pathlib import Path

from DDSim.DD4hepSimulation import DD4hepSimulation
from g4units import GeV, mm

# DDSim 1.36 injects SIM but does not define __file__ when reading steering.
if "__file__" in globals():
    steering_path = Path(__file__).resolve()
else:
    steering_path = Path(SIM.steeringFile).resolve()

SIM = DD4hepSimulation()
SIM.compactFile = str(steering_path.parents[1] / "compact" / "HidraTB26_v01.xml")
SIM.runType = "batch"
SIM.numberOfEvents = 1
SIM.outputFile = "HidraTB26.edm4hep.root"
SIM.physicsList = "FTFP_BERT"
SIM.enableGun = True
SIM.gun.particle = "e-"
SIM.gun.energy = 10 * GeV
SIM.gun.position = (0, 0, -1300 * mm)
SIM.gun.direction = (0, 0, 1)
SIM.random.seed = 12345
SIM.part.userParticleHandler = ""
SIM.action.mapActions["HidraPrototype"] = "DRTubesSDAction"
# Optical steps can have zero energy deposition and must reach the action.
SIM.filter.calo = ""


def setup_cerenkov(kernel):
    """Same Cherenkov/optical setup as SteeringFile_IDEA_o2_v01.py."""
    from DDG4 import PhysicsList

    sequence = kernel.physicsList()
    cerenkov = PhysicsList(kernel, "Geant4CerenkovPhysics/CerenkovPhys")
    cerenkov.MaxNumPhotonsPerStep = 1000
    cerenkov.VerboseLevel = 0
    cerenkov.enableUI()
    sequence.adopt(cerenkov)
    optical = PhysicsList(kernel, "Geant4OpticalPhotonPhysics/OpticalGammaPhys")
    optical.addParticleConstructor("G4OpticalPhoton")
    optical.VerboseLevel = 0
    optical.enableUI()
    sequence.adopt(optical)


SIM.physics.setupUserPhysics(setup_cerenkov)
