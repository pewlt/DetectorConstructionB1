#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

namespace B1
{

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void MyDetectorConstruction::ConstructSDandField() override;
  private:
    //Methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
    G4int  fNofLayers = -1;

};

}


#endif
