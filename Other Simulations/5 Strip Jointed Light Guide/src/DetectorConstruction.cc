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
  
  // rect_left
  G4Box* rect_left = new G4Box("Rect_left", 1.5*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  rect_left_LV = new G4LogicalVolume(rect_left, fTankMaterial, "Rect_left", 0, 0, 0);
  rect_left_PV = new G4PVPlacement(0, G4ThreeVector(-4*cm, 0, 2.5*cm + length6*cm), rect_left_LV, "Rect_left", fWorld_LV, false, 0);

  // rect_left2
  G4Box* rect_left2 = new G4Box("Rect_left2", 0.99*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  G4LogicalVolume* rect_left2_LV = new G4LogicalVolume(rect_left2, fTankMaterial, "Rect_left2", 0, 0, 0);
  G4PVPlacement* rect_left2_PV = new G4PVPlacement(0, G4ThreeVector(-6.5*cm, 0, 2.5*cm + length6*cm), rect_left2_LV, "Rect_left2", fWorld_LV, false, 0);

  // rect_right
  G4Box* rect_right = new G4Box("Rect_right", 1.5*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  rect_right_LV = new G4LogicalVolume(rect_right, fTankMaterial, "Rect_right", 0, 0, 0);
  rect_right_PV = new G4PVPlacement(0, G4ThreeVector(4*cm, 0, 2.5*cm + length6*cm), rect_right_LV, "Rect_right", fWorld_LV, false, 0);

  // rect_right2
  G4Box* rect_right2 = new G4Box("Rect_right2", 0.99*fTank_x, 0.25*fTank_y, 2.5*fTank_z);
  G4LogicalVolume* rect_right2_LV = new G4LogicalVolume(rect_right2, fTankMaterial, "Rect_right2", 0, 0, 0);
  G4PVPlacement* rect_right2_PV = new G4PVPlacement(0, G4ThreeVector(6.5*cm, 0, 2.5*cm + length6*cm), rect_right2_LV, "Rect_right2", fWorld_LV, false, 0);


  // s_left1
  G4double length1 = (33/8);
  G4Tubs* s_left1 = new G4Tubs("S_left1", 21*cm, 24*cm, 10*cm, 0*deg, 27.27*deg);
  s_left1_LV = new G4LogicalVolume(s_left1, fTankMaterial, "S_left1", 0, 0, 0);
  G4RotationMatrix* Rot1 = new G4RotationMatrix();
  Rot1->rotateX(0.*deg);
  //s_left1_PV = new G4PVPlacement(Rot1, G4ThreeVector(0, 0, 0), s_left1_LV, "S_left1", fWorld_LV, false, 0);

  // s_left2
  G4Tubs* s_left2 = new G4Tubs("S_left2", 21*cm, 24*cm, 10*cm, 180*deg, 27.27*deg);
  s_left2_LV = new G4LogicalVolume(s_left2, fTankMaterial, "S_left2", 0, 0, 0);
  G4RotationMatrix* Rot2 = new G4RotationMatrix();
  Rot2->rotateX(0.0*deg);
  G4double angle_1 = 30;
  G4double x_1 = 45 * std::cos(27.27*deg);
  G4double y_1 = 45 * std::sin(27.27*deg);
  //s_left2_PV = new G4PVPlacement(0, G4ThreeVector(x_1*cm, y_1*cm, 0), s_left2_LV, "S_left2", fWorld_LV, false, 0);


  // s_right1
  G4Tubs* s_right1 = new G4Tubs("S_right1", 9.5*cm, 10*cm, 1.5*cm, 0.*deg, 13.5*deg);
  s_right1_LV = new G4LogicalVolume(s_right1, fTankMaterial, "S_right1", 0, 0, 0);
  G4RotationMatrix* Rot3 = new G4RotationMatrix();
  Rot3->rotateY(0.*deg);
  //s_right1_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), s_right1_LV, "S_right1", fWorld_LV, false, 0);

  // s_right2
  G4Tubs* s_right2 = new G4Tubs("S_right2", 9.5*cm, 10*cm, 1.5*cm, 180.*deg, 13.5*deg);
  s_right2_LV = new G4LogicalVolume(s_right2, fTankMaterial, "S_right2", 0, 0, 0);
  G4RotationMatrix* Rot4 = new G4RotationMatrix();
  Rot4->rotateY(-90.*deg);
  G4double x_2 = 19.5 * std::cos(13.5*deg);
  G4double y_2 = 19.5 * std::sin(13.5*deg);
  //s_right2_PV = new G4PVPlacement(0, G4ThreeVector(x_2*cm, y_2*cm, 0), s_right2_LV, "S_right2", fWorld_LV, false, 0);


  // s_far_right
  G4Tubs* s_far_right = new G4Tubs("S_far_right", 4.8125*cm, 5.3125*cm, 0.99*cm, 0.*deg, 26.72*deg);
  G4LogicalVolume* s_far_right_LV = new G4LogicalVolume(s_far_right, fTankMaterial, "S_far_right", 0, 0, 0);
  //G4PVPlacement* s_far_right_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 25*cm), s_far_right_LV, "S_far_right", fWorld_LV, false, 0);



  // s_far_left
  G4Tubs* s_far_left = new G4Tubs("S_far_left", 4.8125*cm, 5.3125*cm, 0.99*cm, 180.*deg, 26.72*deg);
  G4LogicalVolume* s_far_left_LV = new G4LogicalVolume(s_far_left, fTankMaterial, "S_far_left", 0, 0, 0);
  G4double x_3 = 10.125 * std::cos(26.72*deg);
  G4double y_3 = 10.125 * std::sin(26.72*deg);
  //G4PVPlacement* s_far_left_PV = new G4PVPlacement(0, G4ThreeVector(x_3*cm, y_3*cm, 25*cm), s_far_left_LV, "S_far_left", fWorld_LV, false, 0);




  // s_far_right2
  G4Tubs* s_far_right2 = new G4Tubs("S_far_right2", 21.51*cm, 23.49*cm, 0.25*cm, 0.*deg, 27.27*deg);
  G4LogicalVolume* s_far_right2_LV = new G4LogicalVolume(s_far_right2, fTankMaterial, "S_far_right2", 0, 0, 0);
  //G4PVPlacement* s_far_right2_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 25*cm), s_far_right2_LV, "S_far_right2", fWorld_LV, false, 0);



  // s_far_left2
  G4Tubs* s_far_left2 = new G4Tubs("S_far_left2", 21.51*cm, 23.49*cm, 0.25*cm, 180.*deg, 27.27*deg);
  G4LogicalVolume* s_far_left2_LV = new G4LogicalVolume(s_far_left2, fTankMaterial, "S_far_left2", 0, 0, 0);
  G4double x_4 = 45 * std::cos(27.27*deg);
  G4double y_4 = 45 * std::sin(27.27*deg);
  //G4PVPlacement* s_far_left2_PV = new G4PVPlacement(0, G4ThreeVector(x_4*cm, y_4*cm, 25*cm), s_far_left2_LV, "S_far_left2", fWorld_LV, false, 0);





  // Box
  G4double length2 = ((5*sqrt(17)) - ((sqrt(77) * 19.5) / 39)) / 2;
  G4double length3 = ((sqrt(77) * 19.5) / 39);
  G4Box* rec_box = new G4Box("Rec_box", 0.25*cm, 22.5*cm, 10*cm);
  rec_box_LV = new G4LogicalVolume(rec_box, fTankMaterial, "Rec_box", 0, 0, 0);
  //rec_box_PV = new G4PVPlacement(0, G4ThreeVector(17.5*cm, length3*cm + length2*cm, 5*cm), rec_box_LV, "Rec_box", fWorld_LV, false, 0);



  // union left 1 and 2
  G4VSolid* union1 = new G4UnionSolid("Union1", s_left1, s_left2, 0, G4ThreeVector(x_1*cm, y_1*cm, 0));
  G4LogicalVolume* union1_LV = new G4LogicalVolume(union1, fTankMaterial, "Union1");
  //G4PVPlacement* union1_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 0), union1_LV, "Union1", fWorld_LV, false, 0);


  // union right 1 and 2
  G4VSolid* union2 = new G4UnionSolid("Union2", s_right1, s_right2, 0, G4ThreeVector(x_2*cm, y_2*cm, 0));
  G4LogicalVolume* union2_LV = new G4LogicalVolume(union2, fTankMaterial, "Union2");
  G4RotationMatrix* Rot_new = new G4RotationMatrix();
  Rot_new->rotateZ(90.*deg);
  Rot_new->rotateX(90.*deg);
  G4RotationMatrix* Rot_new2 = new G4RotationMatrix();
  Rot_new2->rotateZ(270.*deg);
  Rot_new2->rotateX(90.*deg);

  G4PVPlacement* union2_PV = new G4PVPlacement(Rot_new, G4ThreeVector(4*cm, 9.75*cm, length6*cm), union2_LV, "Union2", fWorld_LV, false, 0);
  G4PVPlacement* union2_PV2 = new G4PVPlacement(Rot_new2, G4ThreeVector(-4*cm, -9.75*cm, length6*cm), union2_LV, "Union2_2", fWorld_LV, false, 0);



  // union 3
  G4VSolid* union3 = new G4UnionSolid("Union3", s_far_right, s_far_left, 0, G4ThreeVector(x_3*cm, y_3*cm, 0));
  G4LogicalVolume* union3_LV = new G4LogicalVolume(union3, fTankMaterial, "Union3");
  G4PVPlacement* union3_PV = new G4PVPlacement(Rot_new, G4ThreeVector(6.5*cm, 5.0625*cm, length6*cm), union3_LV, "Union3", fWorld_LV, false, 0);
  G4PVPlacement* union3_PV2 = new G4PVPlacement(Rot_new2, G4ThreeVector(-6.5*cm, -5.0625*cm, length6*cm), union3_LV, "Union3_2", fWorld_LV, false, 0);



  // intersect rect
  G4double length5 = -2.5 - (5 * sqrt(17));
  G4double length7 = (19.5 * std::cos(13*deg)) - (19.5 * std::cos(13.5*deg));


  G4VSolid* intersect_rect = new G4IntersectionSolid("Intersect_rect", union1, rec_box, Rot4, G4ThreeVector(17.5*cm, length3*cm + length2*cm, 5*cm));
  G4LogicalVolume* intersect_rect_LV = new G4LogicalVolume(intersect_rect, fTankMaterial, "Intersect_rect");
  G4RotationMatrix* Rot_6 = new G4RotationMatrix();
  Rot_6->rotateX(-90.*deg);
  Rot_6->rotateZ(180.*deg);
  G4RotationMatrix* Rot_7 = new G4RotationMatrix();
  Rot_7->rotateX(90.*deg);
  //G4PVPlacement* intersect_rect_PV = new G4PVPlacement(Rot_7, G4ThreeVector(-18.5*cm, -4.5*cm + length7*cm, -2.5*cm - 2*length4*cm), intersect_rect_LV, "Intersect_rect", fWorld_LV, false, 0);
  //G4PVPlacement* intersect_rect2_PV = new G4PVPlacement(Rot_6, G4ThreeVector(18.5*cm, 4.5*cm - length7*cm, -2.5*cm - 2*length4*cm), intersect_rect_LV, "Intersect_rect2", fWorld_LV, false, 0);



  // union4 

  G4VSolid* union4 = new G4UnionSolid("Union4", s_far_right2, s_far_left2, 0, G4ThreeVector(x_4*cm, y_4*cm, 0));
  G4LogicalVolume* union4_LV = new G4LogicalVolume(union4, fTankMaterial, "Union4");
  //G4PVPlacement* union4_PV = new G4PVPlacement(Rot_7, G4ThreeVector(-16*cm, 1.08*cm, -2.5*cm - 2*length4*cm), union4_LV, "Union4", fWorld_LV, false, 0);
  //G4PVPlacement* union4_PV2 = new G4PVPlacement(Rot_6, G4ThreeVector(16*cm, -1.08*cm, -2.5*cm - 2*length4*cm), union4_LV, "Union4_2", fWorld_LV, false, 0);


  // cylinder
  G4Tubs* cone = new G4Tubs("Cone", 0., 3.0*cm, 20.0*cm, 0.*deg, 360.0*deg);
  cone_LV = new G4LogicalVolume(cone, fWorldMaterial, "Cone", 0, 0, 0);
  cone_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, -45*cm), cone_LV, "Cone", fWorld_LV, false, 0);




  //sub_box
  G4Box* sub_box = new G4Box("Sub_box", 10*cm, 12.5*cm, 10*cm);
  G4LogicalVolume* sub_box_LV = new G4LogicalVolume(sub_box, fTankMaterial, "Sub_box", 0, 0, 0);
  //G4PVPlacement* sub_box_PV = new G4PVPlacement(0, G4ThreeVector(40*cm, 10*cm, 0), sub_box_LV, "Sub_box", fWorld_LV, false, 0);

  //curves
  G4Tubs* curve1 = new G4Tubs("Curve1", 45*cm, 48*cm, 0.25*cm, 0.*deg, 30*deg);
  G4LogicalVolume* curve1_LV = new G4LogicalVolume(curve1, fTankMaterial, "Curve1", 0, 0, 0);

  G4Tubs* curve2 = new G4Tubs("Curve2", 46.02*cm, 48*cm, 0.25*cm, 0.*deg, 30*deg);
  G4LogicalVolume* curve2_LV = new G4LogicalVolume(curve2, fTankMaterial, "Curve2", 0, 0, 0);

  G4RotationMatrix* Rot_right = new G4RotationMatrix();
  Rot_right->rotateX(90.*deg);
  G4RotationMatrix* Rot_left = new G4RotationMatrix();
  Rot_left->rotateX(-90.*deg);
  Rot_left->rotateZ(180.*deg);
  //G4PVPlacement* curve1_PV = new G4PVPlacement(Rot_right, G4ThreeVector(-42.5*cm, 0.5*cm + length7*cm, -2.5*cm - 2*length4*cm), curve1_LV, "Curve1", fWorld_LV, false, 0);  

  G4VSolid* curve1_sub = new G4IntersectionSolid("Curve1_sub", curve1, sub_box, 0, G4ThreeVector(40*cm, 8.12*cm, 0));
  G4LogicalVolume* curve1_sub_LV = new G4LogicalVolume(curve1_sub, fTankMaterial, "Curve1_sub");
  G4PVPlacement* curve1_sub_PV = new G4PVPlacement(Rot_right, G4ThreeVector(-42.5*cm, 0.5*cm + length7*cm, -2.5*cm - 2*length4*cm), curve1_sub_LV, "Curve1_sub", fWorld_LV, false, 0);
  G4PVPlacement* curve1_sub_PV2 = new G4PVPlacement(Rot_left, G4ThreeVector(42.5*cm, -0.5*cm - length7*cm, -2.5*cm - 2*length4*cm), curve1_sub_LV, "Curve1_sub2", fWorld_LV, false, 0);

  G4VSolid* curve2_sub = new G4IntersectionSolid("Curve2_sub", curve2, sub_box, 0, G4ThreeVector(40*cm, 8.12*cm, 0));
  G4LogicalVolume* curve2_sub_LV = new G4LogicalVolume(curve2_sub, fTankMaterial, "Curve2_sub");
  G4PVPlacement* curve2_sub_PV = new G4PVPlacement(Rot_right, G4ThreeVector(-40.51*cm, 1.08*cm, -2.5*cm - 2*length4*cm), curve2_sub_LV, "Curve2_sub", fWorld_LV, false, 0);
  G4PVPlacement* curve2_sub_PV2 = new G4PVPlacement(Rot_left, G4ThreeVector(40.51*cm, -1.08*cm, -2.5*cm - 2*length4*cm), curve2_sub_LV, "Curve2_sub2", fWorld_LV, false, 0);



  // scintillator
  G4Box* scint = new G4Box("Scint", 7.5*cm, 0.25*cm, 7.5*cm);
  G4LogicalVolume* scint_LV = new G4LogicalVolume(scint, fTankMaterial, "Scint", 0, 0, 0);
  scint_PV = new G4PVPlacement(0, G4ThreeVector(0, 0, 12.5*cm + length6*cm), scint_LV, "Scint", fWorld_LV, false, 0);






  // ------------- Surface --------------

  G4LogicalBorderSurface* surface1 =
    new G4LogicalBorderSurface("Surface1", rect_mid_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface2 =
    new G4LogicalBorderSurface("Surface2", rect_left_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface3 =
    new G4LogicalBorderSurface("Surface3", rect_right_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface4 =
    new G4LogicalBorderSurface("Surface4", union2_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface5 =
    new G4LogicalBorderSurface("Surface5", union2_PV2, world_PV, fSurface);

  G4LogicalBorderSurface* surface6 =
    new G4LogicalBorderSurface("Surface6", curve1_sub_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface7 =
    new G4LogicalBorderSurface("Surface7", curve1_sub_PV2, world_PV, fSurface);

  G4LogicalBorderSurface* surface8 =
    new G4LogicalBorderSurface("Surface8", rect_left2_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface9 =
    new G4LogicalBorderSurface("Surface9", rect_right2_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface10 =
    new G4LogicalBorderSurface("Surface10", cone_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface11 =
    new G4LogicalBorderSurface("Surface11", union3_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface12 =
    new G4LogicalBorderSurface("Surface12", union3_PV2, world_PV, fSurface);

  G4LogicalBorderSurface* surface13 =
    new G4LogicalBorderSurface("Surface13", curve2_sub_PV, world_PV, fSurface);

  G4LogicalBorderSurface* surface14 =
    new G4LogicalBorderSurface("Surface14", curve2_sub_PV2, world_PV, fSurface);

  G4LogicalBorderSurface* surface15 =
    new G4LogicalBorderSurface("Surface15", scint_PV, world_PV, fSurface);

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

