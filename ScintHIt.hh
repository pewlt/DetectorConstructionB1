#ifndef SLABHIT_HH
#define SLABHIT_HH


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"


class SlabHit
{
    public:
        SlabHit() = default;
        SlabHit(const CalorHit&) = default;
        ~SlabHit() override = default;

        
        G4double GetEdep() {return fEdep;};
        void SetEdep(double edep) {fEdep = edep;};
        
        G4double GetHitTime() {return hitTime;};
        void SetHitTime(double hit_time) {hitTime = hit_time;};
        
        G4int GetChannelNbr() {return channelNbr;};
        void SetChannelNbr(G4int channel_number) {channelNbr = channel_number;};
        
        G4int GetPanelNbr() {return panelNbr;};
        void SetPanelNbr(G4int panel_number) {panelNbr = panel_number;};
    
    private:
        G4double fEdep; // Energy deposition
        G4double hitTime; 
        G4int channelNbr;
        G4int panelNbr;
    
};

using SlabHitsCollection = G4THitsCollection<SlabHit>;

#endif
