#ifndef B1DetectorConstruction_h //A CHANGER ?? Je ne sais pas quoi mettre
#define B1DetectorConstruction_h 1

#include "G4VSensitiveDetector.hh"

#include "ScintHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

namespace B1
{

/// Scitillator sensitive detector class
///
/// In Initialize(), it creates one hit for each slab of the scintillator
/// hit for accounting the total quantities in all layers.
///
/// The values are accounted in hits in ProcessHits() function which is called
/// by Geant4 kernel at each step.

class ScintillatorSD : public G4VSensitiveDetector
{
  public:
    ScintillatorSD(const G4String& name,
                  const G4String& hitsCollectionName,
                  G4int nofCells);
    ~ScintillatorSD() override = default;

    // methods from base class
    void   Initialize(G4HCofThisEvent* hitCollection) override;
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    void   EndOfEvent(G4HCofThisEvent* hitCollection) override;

  private:
    SlabHitsCollection* fHitsCollection = nullptr;
    G4int fNofCells = 0;
};

}



#endif
