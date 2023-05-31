


/// \file B1bis/src/SlabHit.hh
/// \brief Implementation of the B1::SlabHit class



#ifndef SlabHit_h
#define SlabHit_h


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
using namespace std;

namespace B1
{
class SlabHit : public G4VHit
{
    public:
        SlabHit() = default;
        SlabHit(const SlabHit&) = default;
        ~SlabHit() override = default;

        
        G4double GetEdep() {return fEdep;};
        void SetEdep(double edep) {fEdep += edep;};
        
        G4double GetHitTime() {return hitTime;};
        void SetHitTime(double hit_time) {hitTime += hit_time;};
        
        vector<G4int> GetSlabVector() {return SlabVector;};
        void SetSlabVector(G4int Slab_number) {SlabVector.push_back(Slab_number);};
        
        G4int GetPanelNbr() {return panelNbr;};
        void SetPanelNbr(G4int panel_number) {panelNbr = panel_number;};
        
        
        // Operators -> Utiles ? 
        
        G4double GetEdepTot() {return fEdepTot;};
//        void SetEdepTot(double edepTot) {fEdepTot = edepTot;};
        
        G4double GetHitTimeTot() {return hitTimeTot;};
//        void SetHitTimeTot(double hit_timeTot) {hitTimeTot = hit_timeTot;};
        
        void Add(G4double de, G4double dt);
        
            // methods from base class
        void Draw()  override{}
        void Print() override;

        
        
    
    private:
        
        G4double fEdep = 0; // Energy deposited
        G4double hitTime = 0; 
        vector<G4int> SlabVector;
        G4int panelNbr = 0;
        G4double fEdepTot = 0; // Energy deposition totale
        G4double hitTimeTot = 0;
    
};


using SlabHitsCollection = G4THitsCollection<SlabHit>;


inline void SlabHit::Add(G4double de, G4double dt) {
  fEdepTot += de;
  hitTimeTot += dt;
}
}

#endif
