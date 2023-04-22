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
        
        //Operators
        ScintHit& operator = (const ScintHit&) = default;
        G4bool operator == (const ScintHit&) const;
        // A quoi vont servir les operators ? 
        
        inline void* operator new(size_t); //What is size_t ?
        inline void operator delete(void*);
        
        //Methods from base class
        void Draw()  override{}
        void Print() override;

        //Methods to handle data
        void Add(G4double de, G4double dl);

        //Get methods 
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
//+++ tout ce qui concerne CalorHitAllocator -> Recherche/comprendre si on en a besoin 
