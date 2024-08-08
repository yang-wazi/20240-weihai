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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4ReplicatedSlice.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  // Define materials
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() {
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Al");
  nistManager->FindOrBuildMaterial("G4_W");
  nistManager->FindOrBuildMaterial("G4_Galactic");

  G4Material* GAGG = new G4Material("GAGG", 6.63*g/cm3, 4);
  GAGG->AddElement(nistManager->FindOrBuildElement("Gd"), 3);
  GAGG->AddElement(nistManager->FindOrBuildElement("Al"), 2);
  GAGG->AddElement(nistManager->FindOrBuildElement("Ga"), 3);
  GAGG->AddElement(nistManager->FindOrBuildElement("O"), 12);
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {
  // Geometry parameters
  G4double moduleXY = 1.5 * cm;
  G4double holeXY = 1.67 * mm;    // Hole: Volumn occupied by a fiber [ (1.67mm)^2 ]
  G4double fiberXY = 1.0 * mm;
  G4double fiber1Z = 4.5 * cm;
  G4double fiber2Z = 10.5 * cm;
  G4double AlThickness = 1.0 * mm; // Modified thickness for clarification
  G4int nofModules = 40;
  G4int nofFibers = 9;

  auto moduleZ = fiber1Z + AlThickness + fiber2Z;
  auto moduleBound = moduleXY - holeXY * nofFibers;
  auto fiberGap = holeXY - fiberXY;
  auto calorXY = nofModules * moduleXY;
  auto calorZ = moduleZ;
  auto worldSizeXY = 1.2 * calorXY;
  auto worldSizeZ  = 3 * m;   //for Particle Gun in the world (2*m)

  // Get materials
  auto air = G4Material::GetMaterial("G4_AIR");
  auto GAGG = G4Material::GetMaterial("GAGG");
  auto Al = G4Material::GetMaterial("G4_Al");
  auto W = G4Material::GetMaterial("G4_W");
  auto vac = G4Material::GetMaterial("G4_Galactic");
  if ( ! air || ! GAGG || ! Al || ! W || ! vac) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  

  //     
  // World
  //
  auto solidWorld = new G4Box("World", worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2);
  auto logicWorld = new G4LogicalVolume(solidWorld, vac, "World");
  auto physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, true);

  //                               
  // Calorimeter
  //  
  auto solidCalor = new G4Box("Calorimeter", calorXY / 2, calorXY / 2, calorZ / 2);
  auto logicCalor = new G4LogicalVolume(solidCalor, vac, "Calorimeter");
  auto physCalor = new G4PVPlacement(nullptr, G4ThreeVector(0,0,-1*m), logicCalor, "Calorimeter", logicWorld, false, 0, true);

  //                               
  // Module Strip
  //  
  auto solidModuleStrip = new G4Box("ModuleStrip", calorXY / 2, moduleXY / 2, calorZ / 2);
  auto logicModuleStrip = new G4LogicalVolume(solidModuleStrip, W, "ModuleStrip");
  new G4PVReplica("ModuleStrip", logicModuleStrip, logicCalor, kYAxis, nofModules, moduleXY);

  //                               
  // Modules
  //  
  auto solidModule = new G4Box("ModuleStrip", moduleXY / 2, moduleXY / 2, calorZ / 2);
  auto logicModule = new G4LogicalVolume(solidModule, W, "Module");
  new G4PVReplica("Module", logicModule, logicModuleStrip, kXAxis, nofModules, moduleXY);

  //                               
  // Three parts of module
  //  
  auto solidSpa1 = new G4Box("Spa1", moduleXY / 2, moduleXY / 2, fiber1Z / 2);
  auto logicSpa1 = new G4LogicalVolume(solidSpa1, W, "Spa1");
  G4double fiber1PosZ = -moduleZ / 2 + fiber1Z / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fiber1PosZ), logicSpa1, "Spa1", logicModule, false, 0, true);

  auto solidMirror = new G4Box("Mirror", moduleXY / 2, moduleXY / 2, AlThickness / 2);
  auto logicMirror = new G4LogicalVolume(solidMirror, Al, "Mirror");
  G4double mirrorPosZ = fiber1PosZ + fiber1Z / 2 + AlThickness / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, mirrorPosZ), logicMirror, "Mirror", logicModule, false, 0, true);

  auto solidSpa2 = new G4Box("Spa2", moduleXY / 2, moduleXY / 2, fiber2Z / 2);
  auto logicSpa2 = new G4LogicalVolume(solidSpa2, W, "Spa2");
  G4double fiber2PosZ = moduleZ / 2 - fiber2Z / 2;
  new G4PVPlacement(nullptr, G4ThreeVector(0, 0, fiber2PosZ), logicSpa2, "Spa2", logicModule, false, 0, true);

  //
  // Hole Strip
  //
  auto solidHoleStrip1 = new G4Box("HoleStrip1", holeXY * nofFibers / 2, holeXY / 2, fiber1Z / 2);
  auto logicHoleStrip1 = new G4LogicalVolume(solidHoleStrip1, W, "HoleStrip1");
  new G4PVDivision("HoleStrip1", logicHoleStrip1, logicSpa1, kYAxis, nofFibers, holeXY, moduleBound);

  auto solidHoleStrip2 = new G4Box("HoleStrip2", holeXY * nofFibers / 2, holeXY / 2, fiber2Z / 2);
  auto logicHoleStrip2 = new G4LogicalVolume(solidHoleStrip2, W, "HoleStrip2");
  new G4PVDivision("HoleStrip2", logicHoleStrip2, logicSpa2, kYAxis, nofFibers, holeXY, moduleBound);

  //
  // Hole
  //
  auto solidHole1 = new G4Box("Hole1", holeXY / 2, holeXY / 2, fiber1Z / 2);
  auto logicHole1 = new G4LogicalVolume(solidHole1, W, "Hole1");
  new G4PVReplica("Hole1", logicHole1, logicHoleStrip1, kXAxis, nofFibers, holeXY);

  auto solidHole2 = new G4Box("Hole2", holeXY / 2, holeXY / 2, fiber2Z / 2);
  auto logicHole2 = new G4LogicalVolume(solidHole2, W, "Hole2");
  new G4PVReplica("Hole1", logicHole2, logicHoleStrip2, kXAxis, nofFibers, holeXY);

  // Fibers
  auto solidFiber1 = new G4Box("Fiber1", fiberXY / 2, fiberXY / 2, fiber1Z / 2);
  auto logicFiber1 = new G4LogicalVolume(solidFiber1, GAGG, "Fiber1");
  new G4PVPlacement(nullptr, G4ThreeVector(), logicFiber1, "Fiber1", logicHole1, false, 0, true);

  auto solidFiber2 = new G4Box("Fiber2", fiberXY / 2, fiberXY / 2, fiber2Z / 2);
  auto logicFiber2 = new G4LogicalVolume(solidFiber2, GAGG, "Fiber2");
  new G4PVPlacement(nullptr, G4ThreeVector(), logicFiber2, "Fiber2", logicHole2, false, 0, true);

  // Visualization attributes
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto moduleVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red for modules
  moduleVisAtt->SetVisibility(true);
  logicModule->SetVisAttributes(moduleVisAtt);

  auto fiberVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // Green for fibers
  fiberVisAtt->SetVisibility(true);
  logicFiber1->SetVisAttributes(fiberVisAtt);
  logicFiber2->SetVisAttributes(fiberVisAtt);

  auto mirrorVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue for Mirror
  mirrorVisAtt->SetVisibility(true);
  logicMirror->SetVisAttributes(mirrorVisAtt);

  return physWorld;
}

void DetectorConstruction::ConstructSDandField() {
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  G4AutoDelete::Register(fMagFieldMessenger);
}

}  // namespace B4

