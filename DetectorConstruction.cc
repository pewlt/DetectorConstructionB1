//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

namespace B1
{

G4VPhysicalVolume* DetectorConstruction::Construct()
{ 
  //Define Volumes
  return DefineVolumes();
}
  
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{ 
  
    // Geometry parameters
  G4int fNofSlab = 16; //nombre de slab
  G4double SlabLength = 160.*cm; //longueur de la slab
  G4double SlabWidth = 10.*cm; //largeur de la slab
  G4double SlabThickness = 3.*cm; //epaisseur de la slab
  G4double gapThickness = 100.*cm; //ecart entre les deux plaques

  auto TotalThickness = 4 * SlabThickness + gapThickness; //epaisseur totale
  auto TotalWidth = fNofSlab * SlabWidth; //largeur totale
  auto worldSizeX = 1.2 * SlabLength; //longeur du monde primaire
  auto worldSizeY = 1.2 * TotalWidth; //largeur du monde primaire
  auto worldSizeZ = 1.2 * TotalThickness; //epaisseur du monde primaire
  
  
  // Envelope parameters
  //
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool fCheckOverlaps = true;
  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2);

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 world_mat,        // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps
  //
  // Envelope
  //
  auto EnvS = new G4Box("Envelope",                    // its name
    0.5 * SlabLength, 0.5 * TotalWidth, 0.5 * TotalThickness);  // its size

  auto EnvLV = new G4LogicalVolume(EnvS,  // its solid
    env_mat,                                     // its material
    "Envelope");                                 // its name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    EnvLV,                 // its logical volume
    "Envelope",               // its name
    worldLV,               // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);           // overlaps checking

  //
  //Slab
  //

  auto Slab
    = new G4Box("Slab",
                SlabLength/2, SlabWidth/2, SlabThickness/2);

//  auto SlabLV
    auto SlabLV
    = new G4LogicalVolume(
                 Slab,
//                 SLabS,           // its solid
                 slab_mat,  // its material
                 "Slab");         // its name

  //
  //Panel1
  //
  auto Panel1
    = new G4Box("Panel1",
                SlabLength/2, TotalWidth/2, SlabThickness/2);


    auto Panel1LV
    = new G4LogicalVolume(
                 Panel1,
//                 layerS,           // its solid
                 slab_mat,  // its material
                 "Panel1");         // its name

    
  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0 , 0 , -(gapThickness/2 + SlabThickness/2 + SlabThickness)),          // at (0,0,0)
    Panel1LV,                  // its logical volume
    "Panel1",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel1",          // its name
                 SlabLV,          // its logical volume
                 Panel1LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

/// Panel2
  auto Panel2
    = new G4Box("Panel2",
                SlabLength/2, TotalWidth/2, SlabThickness/2);
    
    auto Panel2LV
    = new G4LogicalVolume(
                 Panel2,
                 slab_mat,  // its material
                 "Panel2");         // its name

  
  G4RotationMatrix* r = new G4RotationMatrix;
  r ->rotateZ(M_PI/2.*rad);  
    
  new G4PVPlacement(r,  // no rotation
    G4ThreeVector(0 , 0 , -(gapThickness/2 + SlabThickness/2)),          // at (0,0,0)
    Panel2LV,                  // its logical volume
    "Panel2",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel2",          // its name
                 SlabLV,          // its logical volume
                 Panel2LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

  
  ////Panel3
  auto Panel3
    = new G4Box("Panel3",
                SlabLength/2, TotalWidth/2, SlabThickness/2);


    auto Panel3LV
    = new G4LogicalVolume(
                 Panel3,
                 slab_mat,  // its material
                 "Panel3");         // its name

    
  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0 , 0 , (gapThickness/2 + SlabThickness/2)),          // at (0,0,0)
    Panel3LV,                  // its logical volume
    "Panel3",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  


  new G4PVReplica(
                 "Panel3",          // its name
                 SlabLV,          // its logical volume
                 Panel3LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

/// Panel4
  auto Panel4
    = new G4Box("Panel4",
                SlabLength/2, TotalWidth/2, SlabThickness/2);
    
  auto Panel4LV
    = new G4LogicalVolume(
                 Panel4,
                 slab_mat,  // its material
                 "Panel4");         // its name


    
  new G4PVPlacement(r,  // no rotation
    G4ThreeVector(0 , 0 , (gapThickness/2 + SlabThickness + SlabThickness/2)),          // at (0,0,0)
    Panel4LV,                  // its logical volume
    "Panel4",            // its name
    EnvLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel4",          // its name
                 SlabLV,          // its logical volume
                 Panel4LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica
  
  // A quoi sert le replica au dessus ?? 
 fScoringVolume = Panel4LV; // On a pas besoin du scoring volume donc peut on l'effacer ? 
 
  return worldPV;
}

}
