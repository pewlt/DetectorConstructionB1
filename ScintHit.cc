#include "ScintHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

namespace B1

{

//G4ThreadLocal G4Allocator<CalorHit>* CalorHitAllocator = nullptr;


G4bool ScintHit::operator==(const ScintHit& right) const
{
  return ( this == &right ) ? true : false;
}


void ScintHit::Print()
{
  G4cout
     << "Edep: "
     << std::setw(7) << G4BestUnit(fEdep,"Energy")
     << ”hitTime: ”
     << std::setw(7) << G4BestUnit(hitTime,”hitTime")
     << " track length: "
     << std::setw(7) << G4BestUnit( fTrackLength,"Length")
     << G4endl;
}


}

