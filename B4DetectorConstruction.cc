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
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fAbsorberPV(nullptr),
   fGapPV(nullptr),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  G4Element* Al = nistManager->FindOrBuildMaterial("G4_Al");
  G4Element* Ga = nistManager->FindOrBuildMaterial("G4_Ga");
  G4Element* Gd = nistManager->FindOrBuildMaterial("G4_Gd");
  G4Element* O = nistManager->FindOrBuildMaterial("G4_O");
  nistManager->FindOrBuildMaterial("G4_W");  
  // // Liquid argon material
  // G4double a;  // mass of a mole;
  // G4double z;  // z=mean number of protons;  
  // G4double density; 
  // new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  //        // The argon by NIST Manager is a gas with a different density

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);
  
  // GAGG
  G4Material* matGAGG = new G4Material(name = "GAGG",
                                      density = 6.63 * g / cm3,
                                      ncomponents = 4);
  matGAGG->AddElement(Ga, 3);
  matGAGG->AddElement(Gd, 3);
  matGAGG->AddElement(Al, 2);
  matGAGG->AddElement(O, 12);
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
  G4int nofLayers = 10;
  G4double absoThickness = 10.*mm;
  G4double gapThickness =  5.*mm;
  //G4double calorSizeXY  = 10.*cm;
  G4int nofStrips = 40;
  G4int nofCells = 40;
  G4int nofFibers = 9;
  G4double spaThickness1 = 10.5*cm;
  G4double spaThickness2 = 4.5*cm ;
  G4double mirrorThickness = 1*mm ;

  G4double fiberSizeXY = 1*mm;
  G4double fiberGap = 0.67*mm;
  G4double cellMargin = 0.32*mm;
  G4double cellSizeXY = 1.5*cm;

  auto calorThickness = spaThickness1 + mirrorThickness + spaThickness2;
  auto calorSizeXY = nofCells * cellSizeXY;
  auto worldSizeXY = 1.2 * calorSizeXY;
  auto worldSizeZ  = 1.2 * calorThickness;
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_W");
  auto mirrorMaterial = G4Material::GetMaterial("G4_Al");
  auto fiberMaterial = G4Material::GetMaterial("G4_Pb");

  if ( ! defaultMaterial || ! absorberMaterial || ! mirrorMaterial || ! fiberMaterial) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
  
  // colour
  G4VisAttributes* Yellow = new G4VisAttributes(G4Colour(1.,1.,0.,1.0));
  G4VisAttributes* Blue = new G4VisAttributes(G4Colour(0.,0.,1.,0.3)); 
  G4VisAttributes* Red = new G4VisAttributes(G4Colour(1.,0.,0.,1.0));
  G4VisAttributes* Gray = new G4VisAttributes(G4Colour(0.5,0.5,0.5,1.0));
  Yellow->SetForceSolid(true); 
  Blue->SetForceSolid(true);
  Red->SetForceSolid(true);
  Gray->SetForceSolid(true);
  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                               
  // Calorimeter
  //  
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                         
  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 calorLV,          // its logical volume                         
                 "Calorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //                                 
  // spaOne
  //
   auto spaOneS = new G4Box("SpaOne", calorSizeXY/2, calorSizeXY/2, spaThickness1/2);
  auto spaOneLV = new G4LogicalVolume(
                spaOneS,
                absorberMaterial,
                "SpaOne");
  new G4PVPlacement(nullptr,
    G4ThreeVector(0, 0, -0.5*(calorThickness - spaThickness1)),
    spaOneLV,
    "SpaOne",
    calorLV,
    false,
    0,
    fCheckOverlaps);
  spaOneLV->SetVisAttributes(Blue);
  //                               
  // mirror
  //
  auto mirrorS = new G4Box("Mirror", calorSizeXY/2, calorSizeXY/2, mirrorThickness/2);
  auto mirrorLV = new G4LogicalVolume(
                mirrorS,
                mirrorMaterial,
                "Mirror");
  new G4PVPlacement(nullptr,
    G4ThreeVector(0, 0, 0.5*calorThickness - spaThickness2 - 0.5*mirrorThickness),
    mirrorLV,
    "Mirror",
    calorLV,
    false,
    0,
    fCheckOverlaps);
  mirrorLV->SetVisAttributes(Red);
  //                               
  // sapTwo
  //
  auto spaTwoS = new G4Box("SpaTwo", calorSizeXY/2, calorSizeXY/2, spaThickness2/2);
  auto spaTwoLV = new G4LogicalVolume(
                spaTwoS,
                absorberMaterial,
                "SpaTwo");
  new G4PVPlacement(nullptr,
    G4ThreeVector(0, 0, 0.5*(calorThickness - spaThickness2)),
    spaTwoLV,
    "SpaTwo",
    calorLV,
    false,
    0,
    fCheckOverlaps);
  spaTwoLV->SetVisAttributes(Blue);

  G4LogicalVolume** fLogicCellSpaOne = new G4LogicalVolume*[nofCells*nofCells];
  G4LogicalVolume** fLogicCellSpaTwo = new G4LogicalVolume*[nofCells*nofCells];
  auto cellSpaOneS = new G4Box("CellOne", cellSizeXY/2, cellSizeXY/2, spaThickness1/2);
  auto cellSpaTwoS = new G4Box("CellTwo", cellSizeXY/2, cellSizeXY/2, spaThickness2/2);
  auto fiberSpaOneS= new G4Box("FiberOne",fiberSizeXY/2,fiberSizeXY/2,spaThickness1/2);
  auto fiberSpaTwoS= new G4Box("FiberTwo",fiberSizeXY/2,fiberSizeXY/2,spaThickness2/2);
  auto fiberSpaOneLV = new G4LogicalVolume(fiberSpaOneS, fiberMaterial, "FiberOneLV", nullptr,nullptr,nullptr);
  auto fiberSpaTwoLV = new G4LogicalVolume(fiberSpaTwoS, fiberMaterial, "FiberTwoLV", nullptr,nullptr,nullptr);
  fiberSpaOneLV->SetVisAttributes(Yellow);
  fiberSpaTwoLV->SetVisAttributes(Yellow);

  for (G4int cellNi=0; cellNi<nofCells; cellNi++){
    for (G4int cellNj=0; cellNj<nofCells; cellNj++){
      G4double cellposX = (cellNi+0.5)*cellSizeXY-0.5*calorSizeXY;
      G4double cellposY = (cellNj+0.5)*cellSizeXY-0.5*calorSizeXY;
      G4double cellposZ =  -0.5*(calorThickness - spaThickness1); 
      G4int index = cellNi + cellNj * nofCells ;
      fLogicCellSpaOne[index] =
        new G4LogicalVolume(cellSpaOneS, absorberMaterial, "CellOneLV",nullptr,nullptr,nullptr);
      fLogicCellSpaOne[index]->SetVisAttributes(Blue);
      new G4PVPlacement(nullptr,
        G4ThreeVector(cellposX,cellposY,0),
        fLogicCellSpaOne[index], //fLogicCellSpaOne[cellNi][cellNj],
        "CEllOneLV",
        spaOneLV,
        index,
        fCheckOverlaps);
      fLogicCellSpaTwo[index] =
        new G4LogicalVolume(cellSpaTwoS, absorberMaterial, "CellTwoLV",nullptr,nullptr,nullptr);
      fLogicCellSpaTwo[index]->SetVisAttributes(Blue);
      new G4PVPlacement(nullptr,
        G4ThreeVector(cellposX,cellposY,0),
        fLogicCellSpaTwo[index], //fLogicCellSpaTwo[cellNi][cellNj],
        "CEllTwoLV",
        spaTwoLV,
        index,//cellNi + nofCells * cellNj,
        fCheckOverlaps);

  }}
//
   for (G4int cellNi=0; cellNi<nofCells; cellNi++){
   for (G4int cellNj=0; cellNj<nofCells; cellNj++){
        G4int index = cellNi + nofCells * cellNj;
        for (G4int fiberNi=0; fiberNi<nofFibers; fiberNi++){
        for (G4int fiberNj=0; fiberNj<nofFibers; fiberNj++){
          G4double fiberposX = cellMargin+0.5*fiberSizeXY+fiberNi*(fiberGap+fiberSizeXY)-0.5*cellSizeXY;
          G4double fiberposY = cellMargin+0.5*fiberSizeXY+fiberNj*(fiberGap+fiberSizeXY)-0.5*cellSizeXY;
          G4double fiberposZ = 0.5*(calorThickness - spaThickness2);
          // In SpaOne
          //auto fiberSpaOneLV = new G4LogicalVolume(fiberSpaOneS, absorberMaterial, "FiberOneLV", nullptr,nullptr,nullptr);
          new G4PVPlacement(nullptr,
            G4ThreeVector(fiberposX,fiberposY,0),
            fiberSpaOneLV,
            "FiberOneLV",
            fLogicCellSpaOne[index], //fLogicCellSpaOne[cellNi][cellNj],
            //fLogicCellSpaOne[820],
            0,
            fCheckOverlaps);
        }}
    }}
   for (G4int cellNi=0; cellNi<nofCells; cellNi++){
   for (G4int cellNj=0; cellNj<nofCells; cellNj++){
        G4int index = cellNi + nofCells * cellNj;
        for (G4int fiberNi=0; fiberNi<nofFibers; fiberNi++){
        for (G4int fiberNj=0; fiberNj<nofFibers; fiberNj++){
          G4double fiberposX = cellMargin+0.5*fiberSizeXY+fiberNi*(fiberGap+fiberSizeXY)-0.5*cellSizeXY;
          G4double fiberposY = cellMargin+0.5*fiberSizeXY+fiberNj*(fiberGap+fiberSizeXY)-0.5*cellSizeXY;
          G4double fiberposZ = 0.5*(calorThickness - spaThickness2);
          // In SpaTwo
          //auto fiberSpaTwoLV = new G4LogicalVolume(fiberSpaTwoS, absorberMaterial, "FiberTwoLV", nullptr,nullptr,nullptr);
          new G4PVPlacement(nullptr,
            G4ThreeVector(fiberposX,fiberposY,0),
            fiberSpaTwoLV,
            "FiberTwoLV",
            fLogicCellSpaTwo[index], //fLogicCellSpaTwo[cellNi][cellNj],
            //fLogicCellSpaTwo[0][39],
            0,
            fCheckOverlaps);
        }}// end of fiberNj, fiberNi
    }}
  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << calorThickness/mm << "mm of " << absorberMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
