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
/// \file optical/OpNovice2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "DetectorMessenger.hh"
#include "G4RunManager.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction()
  , fDetectorMessenger(nullptr)
{
  fExpHall_x = fExpHall_y = fExpHall_z = 70.0 * cm;
  fTank_x = fTank_y = fTank_z = 1.0 * cm;

  fTank = nullptr;

  fTankMPT    = new G4MaterialPropertiesTable();
  fScintMPT = new G4MaterialPropertiesTable();
  fWorldMPT   = new G4MaterialPropertiesTable();
  fSurfaceMPT = new G4MaterialPropertiesTable();
  fSurfaceMPT2 = new G4MaterialPropertiesTable();

  fSurface2 = new G4OpticalSurface("Surface2");
  fSurface2->SetType(dielectric_dielectric);
  fSurface2->SetFinish(polished);
  fSurface2->SetModel(unified);

  fSurface = new G4OpticalSurface("Surface");
  fSurface->SetType(dielectric_dielectric);
  fSurface->SetFinish(polished);
  fSurface->SetModel(unified);

  const G4int NUM = 6;
  G4double pp[NUM] = {2.0*eV, 2.2*eV, 2.4*eV, 2.6*eV, 2.8*eV, 3.0*eV};
  G4double rindex[NUM] = {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};
  G4double rindex2[NUM] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  G4double rindex3[NUM] = {1.52, 1.52, 1.52, 1.52, 1.52, 1.52};
  G4double reflectivity[NUM] = {0.3, 0.3, 0.3, 0.3, 0.3, 0.3};

  G4double reflectivity2[NUM] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  G4double tran2[NUM] = {0., 0., 0., 0., 0., 0.};

  G4double tran[NUM] = {0.7, 0.7, 0.7, 0.7, 0.7, 0.7};
  G4double absorption[NUM] = {3.448*m, 4.082 * m,  6.329 * m,  9.174 * m,  12.346 * m, 13.889 * m};
  fTankMPT->AddProperty("RINDEX", pp, rindex, NUM);
  //fTankMPT->AddProperty("ABSLENGTH", pp, absorption, NUM);
  //fSurfaceMPT->AddProperty("REFLECTIVITY",pp,reflectivity,NUM);
  //fSurfaceMPT->AddProperty("TRANSMITTANCE",pp,tran,NUM);

  fSurfaceMPT2->AddProperty("REFLECTIVITY", pp, reflectivity2, NUM);
  fSurfaceMPT2->AddProperty("TRANSMITTANCE", pp, tran2, NUM);

  fWorldMPT->AddProperty("RINDEX", pp, rindex2, NUM);
  fScintMPT->AddProperty("RINDEX", pp, rindex3, NUM);


  fSurface2->SetMaterialPropertiesTable(fSurfaceMPT2);

  fSurface->SetMaterialPropertiesTable(fSurfaceMPT);

  fTank_LV  = nullptr;
  fWorld_LV = nullptr;
  rect_mid_LV = nullptr;
  rect_left_LV = nullptr;
  rect_right_LV = nullptr;
  s_left1_LV = nullptr;
  s_left1UP_LV = nullptr;
  s_left2_LV = nullptr;
  s_left2UP_LV = nullptr;
  s_right1_LV = nullptr;
  s_right2_LV = nullptr;
  cone_LV = nullptr;
  rem_cyl_LV = nullptr;
  rem_cyl2_LV = nullptr;
  rem_cyl3_LV = nullptr;
  rem_cyl4_LV = nullptr;
  rec_box_LV = nullptr;
  

  fTankMaterial  = G4NistManager::Instance()->FindOrBuildMaterial("G4_GLASS_PLATE");
  fWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  fScintMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pyrex_Glass");

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete fTankMPT;
  delete fScintMPT;
  delete fWorldMPT;
  delete fSurfaceMPT;
  delete fSurfaceMPT2;
  delete fSurface;
  delete fSurface2;
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  fTankMaterial->SetMaterialPropertiesTable(fTankMPT);
  fTankMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
  fScintMaterial->SetMaterialPropertiesTable(fScintMPT);
  fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);

  // ------------- Volumes --------------
  // The experimental Hall
  G4Box* world_box = new G4Box("World", fExpHall_x, fExpHall_y, fExpHall_z);

  fWorld_LV = new G4LogicalVolume(world_box, fWorldMaterial, "World", 0, 0, 0);

  G4VPhysicalVolume* world_PV =
    new G4PVPlacement(0, G4ThreeVector(), fWorld_LV, "World", 0, false, 0);

  


  G4double length4 = (22.5 - (5 * sqrt(17))) / 2;
  G4double length6 = (19.5 * std::sin(13.5*deg)) - (19.5 * std::sin(13*deg));


  // The rect_mid
  G4Box* rect_mid = new G4Box("Rect_mid", 2.49*fTank_x, 0.25*fTank_y, 15.0*fTank_z + (length6/2)*cm);

  rect_mid_LV = new G4LogicalVolume(rect_mid, fTankMaterial, "Rect_mid", 0, 0, 0);

  rect_mid_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, -10*cm + (length6/2)*cm), rect_mid_LV, "Rect_mid", fWorld_LV, false, 0);
  
  // rect_left (PINK LEFT)
  G4Box* rect_left = new G4Box("Rect_left", 1.99*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  rect_left_LV = new G4LogicalVolume(rect_left, fTankMaterial, "Rect_left", 0, 0, 0);
  rect_left_PV = new G4PVPlacement(0, G4ThreeVector(-4.5*cm, 0, 2.5*cm + length6*cm), rect_left_LV, "Rect_left", fWorld_LV, false, 0);

  // rect_left_2 (GREEN LEFT)
  G4Box* rect_left_2 = new G4Box("Rect_left_2", 1.74*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  G4LogicalVolume* rect_left_2_LV = new G4LogicalVolume(rect_left_2, fTankMaterial, "Rect_left_2", 0, 0, 0);
  G4PVPlacement* rect_left_2_PV = new G4PVPlacement(0, G4ThreeVector(-8.25*cm, 0, 2.5*cm + length6*cm), rect_left_2_LV, "Rect_left_2", fWorld_LV, false, 0);

  // rect_left2 (RED LEFT)
  G4Box* rect_left2 = new G4Box("Rect_left2", 1.24*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  G4LogicalVolume* rect_left2_LV = new G4LogicalVolume(rect_left2, fTankMaterial, "Rect_left2", 0, 0, 0);
  G4PVPlacement* rect_left2_PV = new G4PVPlacement(0, G4ThreeVector(-11.25*cm, 0, 2.5*cm + length6*cm), rect_left2_LV, "Rect_left2", fWorld_LV, false, 0);

  // rect_right (PINK RIGHT)
  G4Box* rect_right = new G4Box("Rect_right", 1.99*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  rect_right_LV = new G4LogicalVolume(rect_right, fTankMaterial, "Rect_right", 0, 0, 0);
  rect_right_PV = new G4PVPlacement(0, G4ThreeVector(4.5*cm, 0, 2.5*cm + length6*cm), rect_right_LV, "Rect_right", fWorld_LV, false, 0);

  // rect_right_2 (GREEN RIGHT)
  G4Box* rect_right_2 = new G4Box("Rect_right_2", 1.74*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  G4LogicalVolume* rect_right_2_LV = new G4LogicalVolume(rect_right_2, fTankMaterial, "Rect_right_2", 0, 0, 0);
  G4PVPlacement* rect_right_2_PV = new G4PVPlacement(0, G4ThreeVector(8.25*cm, 0, 2.5*cm + length6*cm), rect_right_2_LV, "Rect_right_2", fWorld_LV, false, 0);

  // rect_right2 (RED RIGHT)
  G4Box* rect_right2 = new G4Box("Rect_right2", 1.24*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  G4LogicalVolume* rect_right2_LV = new G4LogicalVolume(rect_right2, fTankMaterial, "Rect_right2", 0, 0, 0);
  G4PVPlacement* rect_right2_PV = new G4PVPlacement(0, G4ThreeVector(11.25*cm, 0, 2.5*cm + length6*cm), rect_right2_LV, "Rect_right2", fWorld_LV, false, 0);


  // s_left1
  G4double length1 = (33/8);
  G4Tubs* s_left1 = new G4Tubs("S_left1", (23.5-2)*cm, (23.5)*cm, 10*cm, 0*deg, 27.27*deg);
  s_left1_LV = new G4LogicalVolume(s_left1, fTankMaterial, "S_left1", 0, 0, 0);
  G4RotationMatrix* Rot1 = new G4RotationMatrix();
  Rot1->rotateX(0.*deg);
  //s_left1_PV = new G4PVPlacement(Rot1, G4ThreeVector(0, 0, 0), s_left1_LV, "S_left1", fWorld_LV, false, 0);

  // s_left2
  G4Tubs* s_left2 = new G4Tubs("S_left2", (23.5-2)*cm, (23.5)*cm, 10*cm, 180*deg, 27.27*deg);
  s_left2_LV = new G4LogicalVolume(s_left2, fTankMaterial, "S_left2", 0, 0, 0);
  G4RotationMatrix* Rot2 = new G4RotationMatrix();
  Rot2->rotateX(0.0*deg);
  G4double angle_1 = 30;
  G4double x_1 = (2*23.5-2) * std::cos(27.27*deg);
  G4double y_1 = (2*23.5-2) * std::sin(27.27*deg);
  //s_left2_PV = new G4PVPlacement(0, G4ThreeVector(x_1*cm, y_1*cm, 0), s_left2_LV, "S_left2", fWorld_LV, false, 0);


  // s_right1 (USED FOR PINK UP/DOWN PIECE)
  G4Tubs* s_right1 = new G4Tubs("S_right1", 9.5*cm, 10*cm, 1.99*cm, 0.*deg, 13.5*deg);
  s_right1_LV = new G4LogicalVolume(s_right1, fTankMaterial, "S_right1", 0, 0, 0);
  G4RotationMatrix* Rot3 = new G4RotationMatrix();
  Rot3->rotateY(0.*deg);
  //s_right1_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), s_right1_LV, "S_right1", fWorld_LV, false, 0);

  // s_right2 (USED FOR PINK UP/DOWN PIECE)
  G4Tubs* s_right2 = new G4Tubs("S_right2", 9.5*cm, 10*cm, 1.99*cm, 180.*deg, 13.5*deg);
  s_right2_LV = new G4LogicalVolume(s_right2, fTankMaterial, "S_right2", 0, 0, 0);
  G4RotationMatrix* Rot4 = new G4RotationMatrix();
  Rot4->rotateY(-90.*deg);
  G4double x_2 = 19.5 * std::cos(13.5*deg);
  G4double y_2 = 19.5 * std::sin(13.5*deg);
  //s_right2_PV = new G4PVPlacement(0, G4ThreeVector(x_2*cm, y_2*cm, 0), s_right2_LV, "S_right2", fWorld_LV, false, 0);


  // s_far_right (USED FOR RED UP/DOWN PIECE)
  G4Tubs* s_far_right = new G4Tubs("S_far_right", (3.8529 - 0.5)*cm, (3.8529)*cm, 1.24*cm, 0.*deg, 39.18*deg);
  G4LogicalVolume* s_far_right_LV = new G4LogicalVolume(s_far_right, fTankMaterial, "S_far_right", 0, 0, 0);
  //G4PVPlacement* s_far_right_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 25*cm), s_far_right_LV, "S_far_right", fWorld_LV, false, 0);



  // s_far_left (USED FOR RED UP/DOWN PIECE)
  G4Tubs* s_far_left = new G4Tubs("S_far_left", (3.8529 - 0.5)*cm, (3.8529)*cm, 1.24*cm, 180.*deg, 39.18*deg);
  G4LogicalVolume* s_far_left_LV = new G4LogicalVolume(s_far_left, fTankMaterial, "S_far_left", 0, 0, 0);
  G4double x_3 = (2*3.8529 - 0.5) * std::cos(39.18*deg);
  G4double y_3 = (2*3.8529 - 0.5) * std::sin(39.18*deg);
  //G4PVPlacement* s_far_left_PV = new G4PVPlacement(0, G4ThreeVector(x_3*cm, y_3*cm, 25*cm), s_far_left_LV, "S_far_left", fWorld_LV, false, 0);



  // s_far_right_2_2 (USED FOR GREEN UP/DOWN PIECE)
  G4Tubs* s_far_right_2_2 = new G4Tubs("S_far_right_2_2", 4.8125*cm, 5.3125*cm, 1.74*cm, 0.*deg, 26.72*deg);
  G4LogicalVolume* s_far_right_2_2_LV = new G4LogicalVolume(s_far_right_2_2, fTankMaterial, "S_far_right_2", 0, 0, 0);
  //G4PVPlacement* s_far_right_2_2_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 25*cm), s_far_right_2_2_LV, "S_far_right_2_2", fWorld_LV, false, 0);



  // s_far_left_2_2 (USED FOR GREEN UP/DOWN PIECE)
  G4Tubs* s_far_left_2_2 = new G4Tubs("S_far_left_2_2", 4.8125*cm, 5.3125*cm, 1.74*cm, 180.*deg, 26.72*deg);
  G4LogicalVolume* s_far_left_2_2_LV = new G4LogicalVolume(s_far_left_2_2, fTankMaterial, "S_far_left_2", 0, 0, 0);
  G4double x_3_2 = 10.125 * std::cos(26.72*deg);
  G4double y_3_2 = 10.125 * std::sin(26.72*deg);
  //G4PVPlacement* s_far_left_2_2_PV = new G4PVPlacement(0, G4ThreeVector(x_3*cm, y_3*cm, 25*cm), s_far_left_2_2_LV, "S_far_left_2_2", fWorld_LV, false, 0);


  // s_far_right2 (red)
  G4Tubs* s_far_right2 = new G4Tubs("S_far_right2", (23.125-1.25)*cm, (23.125)*cm, 0.25*cm, 0.*deg, 27.27*deg);
  G4LogicalVolume* s_far_right2_LV = new G4LogicalVolume(s_far_right2, fTankMaterial, "S_far_right2", 0, 0, 0);
  //G4PVPlacement* s_far_right2_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 25*cm), s_far_right2_LV, "S_far_right2", fWorld_LV, false, 0);



  // s_far_left2 (red)
  G4Tubs* s_far_left2 = new G4Tubs("S_far_left2", (23.125-1.25)*cm, (23.125)*cm, 0.25*cm, 180.*deg, 27.27*deg);
  G4LogicalVolume* s_far_left2_LV = new G4LogicalVolume(s_far_left2, fTankMaterial, "S_far_left2", 0, 0, 0);
  G4double x_4 = (2*23.125 -1.25) * std::cos(27.27*deg);
  G4double y_4 = (2*23.125 -1.25) * std::sin(27.27*deg);
  //G4PVPlacement* s_far_left2_PV = new G4PVPlacement(0, G4ThreeVector(x_4*cm, y_4*cm, 25*cm), s_far_left2_LV, "S_far_left2", fWorld_LV, false, 0);


  // s_far_right_2 (green)
  G4Tubs* s_far_right_2 = new G4Tubs("S_far_right_2", 21.625*cm, 23.375*cm, 0.25*cm, 0.*deg, 27.27*deg);
  G4LogicalVolume* s_far_right_2_LV = new G4LogicalVolume(s_far_right_2, fTankMaterial, "S_far_right_2", 0, 0, 0);
  //G4PVPlacement* s_far_right_2_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 25*cm), s_far_right_2_LV, "S_far_right_2", fWorld_LV, false, 0);



  // s_far_left_2 (green)
  G4Tubs* s_far_left_2 = new G4Tubs("S_far_left_2", 21.625*cm, 23.375*cm, 0.25*cm, 180.*deg, 27.27*deg);
  G4LogicalVolume* s_far_left_2_LV = new G4LogicalVolume(s_far_left_2, fTankMaterial, "S_far_left_2", 0, 0, 0);
  G4double x_4_2 = (2*23.375 -1.75) * std::cos(27.27*deg);
  G4double y_4_2 = (2*23.375 -1.75) * std::sin(27.27*deg);
  //G4PVPlacement* s_far_left_2_PV = new G4PVPlacement(0, G4ThreeVector(x_4*cm, y_4*cm, 25*cm), s_far_left_2_LV, "S_far_left_2", fWorld_LV, false, 0);


  // RED FAR PIECE
  G4Tubs* red_piece1 = new G4Tubs("Red_piece1", (13.505-2.5)*cm, (13.505)*cm, 0.25*cm, 0.*deg, 57.25*deg);
  G4LogicalVolume* red_piece1_LV = new G4LogicalVolume(red_piece1, fTankMaterial, "Red_piece1", 0, 0, 0);
  G4double x_red = (2*13.505 -2.5) * std::cos(57.25*deg);
  G4double y_red = (2*13.505 -2.5) * std::sin(57.25*deg);

  G4Tubs* red_piece2 = new G4Tubs("Red_piece2", (13.505-2.5)*cm, (13.505)*cm, 0.25*cm, 180.*deg, 57.25*deg);
  G4LogicalVolume* red_piece2_LV = new G4LogicalVolume(red_piece2, fTankMaterial, "Red_piece2", 0, 0, 0);

  G4VSolid* red_piece = new G4UnionSolid("Red_piece", red_piece1, red_piece2, 0, G4ThreeVector(x_red*cm, y_red*cm, 0));
  G4LogicalVolume* red_piece_LV = new G4LogicalVolume(red_piece, fTankMaterial, "Red_piece");
  G4RotationMatrix* Rot_red = new G4RotationMatrix();
  Rot_red->rotateX(90.*deg);
  G4RotationMatrix* Rot_red2 = new G4RotationMatrix();
  Rot_red2->rotateX(-90.*deg);
  Rot_red2->rotateZ(180.*deg);
  G4PVPlacement* red_piece_PV = new G4PVPlacement(Rot_red, G4ThreeVector(-1*cm,1.62*cm,-2.5*cm - 2*length4*cm), red_piece_LV, "Red_piece", fWorld_LV, false, 0);
  G4PVPlacement* red_piece_PV2 = new G4PVPlacement(Rot_red2, G4ThreeVector(1*cm,-1.62*cm,-2.5*cm - 2*length4*cm), red_piece_LV, "Red_piece2", fWorld_LV, false, 0);


  // GREEN FAR PIECE
  G4Tubs* green_piece1 = new G4Tubs("Green_piece1", (16.6887-3.5)*cm, (16.6887)*cm, 0.25*cm, 0.*deg, 43.62*deg);
  G4LogicalVolume* green_piece1_LV = new G4LogicalVolume(green_piece1, fTankMaterial, "Green_piece1", 0, 0, 0);
  G4double x_green = (2*16.6887 -3.5) * std::cos(43.62*deg);
  G4double y_green = (2*16.6887 -3.5) * std::sin(43.62*deg);

  G4Tubs* green_piece2 = new G4Tubs("Green_piece2", (16.6887-3.5)*cm, (16.6887)*cm, 0.25*cm, 180.*deg, 43.62*deg);
  G4LogicalVolume* green_piece2_LV = new G4LogicalVolume(green_piece2, fTankMaterial, "Green_piece2", 0, 0, 0);

  G4VSolid* green_piece = new G4UnionSolid("Green_piece", green_piece1, green_piece2, 0, G4ThreeVector(x_green*cm, y_green*cm, 0));
  G4LogicalVolume* green_piece_LV = new G4LogicalVolume(green_piece, fTankMaterial, "Green_piece");
  G4RotationMatrix* Rot_green = new G4RotationMatrix();
  Rot_green->rotateX(90.*deg);
  G4RotationMatrix* Rot_green2 = new G4RotationMatrix();
  Rot_green2->rotateX(-90.*deg);
  Rot_green2->rotateZ(180.*deg);
  G4PVPlacement* green_piece_PV = new G4PVPlacement(Rot_green, G4ThreeVector(-6.68*cm, 1.08*cm, -2.5*cm - 2*length4*cm), green_piece_LV, "Green_piece", fWorld_LV, false, 0);
  G4PVPlacement* green_piece_PV2 = new G4PVPlacement(Rot_green2, G4ThreeVector(6.68*cm, -1.08*cm, -2.5*cm - 2*length4*cm), green_piece_LV, "Green_piece2", fWorld_LV, false, 0);


  G4double length7 = (19.5 * std::cos(13*deg)) - (19.5 * std::cos(13.5*deg));

  // PINK INNER PIECE
  G4Tubs* pink_piece1 = new G4Tubs("Pink_piece1", (26.7314-4)*cm, (26.7314)*cm, 0.25*cm, 0.*deg, 24.63*deg);
  G4LogicalVolume* pink_piece1_LV = new G4LogicalVolume(pink_piece1, fTankMaterial, "Pink_piece1", 0, 0, 0);
  G4double x_pink = (2*26.7314 -4) * std::cos(24.63*deg);
  G4double y_pink = (2*26.7314 -4) * std::sin(24.63*deg);

  G4Tubs* pink_piece2 = new G4Tubs("Pink_piece2", (26.7314-4)*cm, (26.7314)*cm, 0.25*cm, 180.*deg, 24.63*deg);
  G4LogicalVolume* pink_piece2_LV = new G4LogicalVolume(pink_piece2, fTankMaterial, "Pink_piece2", 0, 0, 0);

  G4VSolid* pink_piece = new G4UnionSolid("Pink_piece", pink_piece1, pink_piece2, 0, G4ThreeVector(x_pink*cm, y_pink*cm, 0));
  G4LogicalVolume* pink_piece_LV = new G4LogicalVolume(pink_piece, fTankMaterial, "Pink_piece");
  G4RotationMatrix* Rot_pink = new G4RotationMatrix();
  Rot_pink->rotateX(90.*deg);
  G4RotationMatrix* Rot_pink2 = new G4RotationMatrix();
  Rot_pink2->rotateX(-90.*deg);
  Rot_pink2->rotateZ(180.*deg);
  G4PVPlacement* pink_piece_PV = new G4PVPlacement(Rot_pink, G4ThreeVector(-20.22*cm, 0.54*cm, -2.5*cm - 2*length4*cm), pink_piece_LV, "Pink_piece", fWorld_LV, false, 0);
  G4PVPlacement* pink_piece_PV2 = new G4PVPlacement(Rot_pink2, G4ThreeVector(20.22*cm, -0.54*cm, -2.5*cm - 2*length4*cm), pink_piece_LV, "Pink_piece2", fWorld_LV, false, 0);

  // Box
  G4double length2 = ((5*sqrt(17)) - ((sqrt(77) * 19.5) / 39)) / 2;
  G4double length3 = ((sqrt(77) * 19.5) / 39);
  G4double L = 25 - 19.5*std::sin(13*deg);
  G4Box* rec_box = new G4Box("Rec_box", 0.25*cm, 22.5*cm, 10*cm);
  rec_box_LV = new G4LogicalVolume(rec_box, fTankMaterial, "Rec_box", 0, 0, 0);
  //rec_box_PV = new G4PVPlacement(0, G4ThreeVector(17.5*cm, length3*cm + length2*cm, 5*cm), rec_box_LV, "Rec_box", fWorld_LV, false, 0);



  // union left 1 and 2
  G4VSolid* union1 = new G4UnionSolid("Union1", s_left1, s_left2, 0, G4ThreeVector(x_1*cm, y_1*cm, 0));
  G4LogicalVolume* union1_LV = new G4LogicalVolume(union1, fTankMaterial, "Union1");
  //G4PVPlacement* union1_PV = new G4PVPlacement(0, G4ThreeVector(0,0,0), union1_LV, "Union1", fWorld_LV, false, 0);


  // union right 1 and 2 (PINK UP/DOWN PIECE)
  G4VSolid* union2 = new G4UnionSolid("Union2", s_right1, s_right2, 0, G4ThreeVector(x_2*cm, y_2*cm, 0));
  G4LogicalVolume* union2_LV = new G4LogicalVolume(union2, fTankMaterial, "Pink_piece");
  G4RotationMatrix* Rot_new = new G4RotationMatrix();
  Rot_new->rotateZ(90.*deg);
  Rot_new->rotateX(90.*deg);
  G4RotationMatrix* Rot_new2 = new G4RotationMatrix();
  Rot_new2->rotateZ(270.*deg);
  Rot_new2->rotateX(90.*deg);


  G4PVPlacement* union_2_PV = new G4PVPlacement(Rot_new, G4ThreeVector(4.5*cm, 9.75*cm, length6*cm), union2_LV, "Union_2", fWorld_LV, false, 0);
  G4PVPlacement* union_2_PV2 = new G4PVPlacement(Rot_new2, G4ThreeVector(-4.5*cm, -9.75*cm, length6*cm), union2_LV, "Union_2_2", fWorld_LV, false, 0);


  // union 3 (RED UP/DOWN PIECE)
  G4VSolid* union3 = new G4UnionSolid("Union3", s_far_right, s_far_left, 0, G4ThreeVector(x_3*cm, y_3*cm, 0));
  G4LogicalVolume* union3_LV = new G4LogicalVolume(union3, fTankMaterial, "Union3");
  G4PVPlacement* union3_PV = new G4PVPlacement(Rot_new, G4ThreeVector(11.25*cm, (3.8529 - 0.25)*cm, length6*cm), union3_LV, "Union3", fWorld_LV, false, 0);
  G4PVPlacement* union3_PV2 = new G4PVPlacement(Rot_new2, G4ThreeVector(-11.25*cm, -(3.8529 - 0.25)*cm, length6*cm), union3_LV, "Union3_2", fWorld_LV, false, 0);


  // union3 2 (GREEN UP/DOWN BEND PIECE)
  G4VSolid* union3_2 = new G4UnionSolid("Union3_2", s_far_right_2_2, s_far_left_2_2, 0, G4ThreeVector(x_3_2*cm, y_3_2*cm, 0));
  G4LogicalVolume* union3_2_LV = new G4LogicalVolume(union3_2, fTankMaterial, "Union2");
  G4PVPlacement* union3_2_PV = new G4PVPlacement(Rot_new2, G4ThreeVector(-8.25*cm, -5.0625*cm, length6*cm), union3_2_LV, "Union32_2", fWorld_LV, false, 0);
  G4PVPlacement* union3_2_PV2 = new G4PVPlacement(Rot_new, G4ThreeVector(8.25*cm, 5.0625*cm, length6*cm), union3_2_LV, "Union3_2_2", fWorld_LV, false, 0);



  // intersect rect
  G4double length5 = -2.5 - (5 * sqrt(17));


  G4VSolid* intersect_rect = new G4IntersectionSolid("Intersect_rect", union1, rec_box, Rot4, G4ThreeVector(17.5*cm, length3*cm + length2*cm, 5*cm));
  G4LogicalVolume* intersect_rect_LV = new G4LogicalVolume(intersect_rect, fTankMaterial, "Intersect_rect");
  G4RotationMatrix* Rot_6 = new G4RotationMatrix();
  Rot_6->rotateX(-90.*deg);
  Rot_6->rotateZ(180.*deg);
  G4RotationMatrix* Rot_7 = new G4RotationMatrix();
  Rot_7->rotateX(90.*deg);
  //G4PVPlacement* intersect_rect_PV = new G4PVPlacement(Rot_7, G4ThreeVector(-17.75*cm, -4.5*cm + length7*cm, -2.5*cm - 2*length4*cm), intersect_rect_LV, "Intersect_rect", fWorld_LV, false, 0);
  //G4PVPlacement* intersect_rect2_PV = new G4PVPlacement(Rot_6, G4ThreeVector(19*cm, 4.5*cm - length7*cm, -2.5*cm - 2*length4*cm), intersect_rect_LV, "Intersect_rect2", fWorld_LV, false, 0);
  //G4PVPlacement* intersect_rect_2_PV = new G4PVPlacement(Rot_7, G4ThreeVector(-19*cm, -4.5*cm + length7*cm, -2.5*cm - 2*length4*cm), intersect_rect_LV, "Intersect_rect_2", fWorld_LV, false, 0);
  //G4PVPlacement* intersect_rect_2_PV2 = new G4PVPlacement(Rot_6, G4ThreeVector(17.75*cm, 4.5*cm - length7*cm, -2.5*cm - 2*length4*cm), intersect_rect_LV, "Intersect_rect_2_2", fWorld_LV, false, 0);


  // union4 (farther) (long piece)

  G4VSolid* union4 = new G4UnionSolid("Union4", s_far_right2, s_far_left2, 0, G4ThreeVector(x_4*cm, y_4*cm, 0));
  G4LogicalVolume* union4_LV = new G4LogicalVolume(union4, fTankMaterial, "Union4");
  //G4PVPlacement* union4_PV = new G4PVPlacement(Rot_7, G4ThreeVector(-15.625*cm, 1.62*cm, -2.5*cm - 2*length4*cm), union4_LV, "Union4", fWorld_LV, false, 0);
  //G4PVPlacement* union4_PV2 = new G4PVPlacement(Rot_6, G4ThreeVector(15.625*cm, -1.62*cm, -2.5*cm - 2*length4*cm), union4_LV, "Union4_2", fWorld_LV, false, 0);


  //union4 2 (closer) (long piece)
  G4VSolid* union4_2 = new G4UnionSolid("Union4_2", s_far_right_2, s_far_left_2, 0, G4ThreeVector(x_4_2*cm, y_4_2*cm, 0));
  G4LogicalVolume* union4_2_LV = new G4LogicalVolume(union4_2, fTankMaterial, "Union2");
  //G4PVPlacement* union4_2_PV = new G4PVPlacement(Rot_7, G4ThreeVector(-17.125*cm, 1.08*cm, -2.5*cm - 2*length4*cm), union4_2_LV, "Union4_2", fWorld_LV, false, 0);
  //G4PVPlacement* union4_2_PV2 = new G4PVPlacement(Rot_6, G4ThreeVector(17.125*cm, -1.08*cm, -2.5*cm - 2*length4*cm), union4_2_LV, "Union4_2_2", fWorld_LV, false, 0);

  // scintillator
  G4Box* scint = new G4Box("Scint", 12.5*cm, 0.25*cm, 7.5*cm);
  G4LogicalVolume* scint_LV = new G4LogicalVolume(scint, fTankMaterial, "Scint", 0, 0, 0);
  scint_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 12.5*cm + length6*cm), scint_LV, "Scint", fWorld_LV, false, 0);



  //box for PMT view picture
  //G4Box* photo = new G4Box("Photo", 10*cm, 10*cm, 5*cm);
  //G4LogicalVolume* photo_LV = new G4LogicalVolume(photo, fTankMaterial, "Photo", 0, 0, 0);
  //G4PVPlacement* photo_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, -19.9*cm), photo_LV, "Photo", fWorld_LV, false, 0);




  // cylinder
  G4Tubs* cone = new G4Tubs("Cone", 0., 2.75*cm, 5.0*cm, 0.*deg, 360.0*deg);
  cone_LV = new G4LogicalVolume(cone, fWorldMaterial, "Cone", 0, 0, 0);
  cone_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, -30*cm), cone_LV, "Cone", fWorld_LV, false, 0);

  //detVolume = cone_LV;




  // ------------- Surface --------------

  //G4LogicalBorderSurface* surface1 =
  //  new G4LogicalBorderSurface("Surface1", rect_mid_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface2 =
  //  new G4LogicalBorderSurface("Surface2", rect_left_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface3 =
  //  new G4LogicalBorderSurface("Surface3", rect_right_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface4 =
  //  new G4LogicalBorderSurface("Surface4", union2_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface5 =
  //  new G4LogicalBorderSurface("Surface5", union2_PV2, world_PV, fSurface);

  //G4LogicalBorderSurface* surface6 =
  //  new G4LogicalBorderSurface("Surface6", intersect_rect_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface7 =
  //  new G4LogicalBorderSurface("Surface7", intersect_rect2_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface8 =
  //  new G4LogicalBorderSurface("Surface8", rect_left2_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface9 =
  //  new G4LogicalBorderSurface("Surface9", rect_right2_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface10 =
  //  new G4LogicalBorderSurface("Surface10", cone_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface11 =
  //  new G4LogicalBorderSurface("Surface11", union3_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface12 =
  //  new G4LogicalBorderSurface("Surface12", union3_PV2, world_PV, fSurface);

  //G4LogicalBorderSurface* surface13 =
  //  new G4LogicalBorderSurface("Surface13", union4_PV, world_PV, fSurface);

  //G4LogicalBorderSurface* surface14 =
  //  new G4LogicalBorderSurface("Surface14", union4_PV2, world_PV, fSurface);

  //G4LogicalBorderSurface* surface15 =
  //  new G4LogicalBorderSurface("Surface15", scint_PV, world_PV, fSurface);

  //G4OpticalSurface* opticalSurface = dynamic_cast<G4OpticalSurface*>(
  //  surface->GetSurface(fTank, world_PV)->GetSurfaceProperty());
  //G4cout << "******  opticalSurface->DumpInfo:" << G4endl;
  //if(opticalSurface)
  //{
  //  opticalSurface->DumpInfo();
  //}
  //G4cout << "******  end of opticalSurface->DumpInfo" << G4endl;

  return world_PV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfaceSigmaAlpha(G4double v)
{
  fSurface->SetSigmaAlpha(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface sigma alpha set to: " << fSurface->GetSigmaAlpha()
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSurfacePolish(G4double v)
{
  fSurface->SetPolish(v);
  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  G4cout << "Surface polish set to: " << fSurface->GetPolish() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddTankMPV(const G4String& prop,
                                      G4MaterialPropertyVector* mpv)
{
  fTankMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the box is now: " << G4endl;
  fTankMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPV(const G4String& prop,
                                       G4MaterialPropertyVector* mpv)
{
  fWorldMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPV(const G4String& prop,
                                         G4MaterialPropertyVector* mpv)
{
  fSurfaceMPT->AddProperty(prop, mpv);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddTankMPC(const G4String& prop, G4double v)
{
  fTankMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the box is now: " << G4endl;
  fTankMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddWorldMPC(const G4String& prop, G4double v)
{
  fWorldMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the world is now: " << G4endl;
  fWorldMPT->DumpTable();
  G4cout << "............." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddSurfaceMPC(const G4String& prop, G4double v)
{
  fSurfaceMPT->AddConstProperty(prop, v);
  G4cout << "The MPT for the surface is now: " << G4endl;
  fSurfaceMPT->DumpTable();
  G4cout << "............." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && fWorldMaterial != pmat)
  {
    fWorldMaterial = pmat;
    if(fWorld_LV)
    {
      fWorld_LV->SetMaterial(fWorldMaterial);
      fWorldMaterial->SetMaterialPropertiesTable(fWorldMPT);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "World material set to " << fWorldMaterial->GetName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetTankMaterial(const G4String& mat)
{
  G4Material* pmat = G4NistManager::Instance()->FindOrBuildMaterial(mat);
  if(pmat && fTankMaterial != pmat)
  {
    fTankMaterial = pmat;
    if(fTank_LV)
    {
      fTank_LV->SetMaterial(fTankMaterial);
      fTankMaterial->SetMaterialPropertiesTable(fTankMPT);
      fTankMaterial->GetIonisation()->SetBirksConstant(0.126 * mm / MeV);
    }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    G4cout << "Tank material set to " << fTankMaterial->GetName() << G4endl;
  }
}

