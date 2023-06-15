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
// The code was written by :
//	^Claudio Andenna  claudio.andenna@ispesl.it, claudio.andenna@iss.infn.it
//      *Barbara Caccia barbara.caccia@iss.it
//      with the support of Pablo Cirrone (LNS, INFN Catania Italy)
//	with the contribute of Alessandro Occhigrossi*
//
// ^INAIL DIPIA - ex ISPESL and INFN Roma, gruppo collegato Sanità, Italy
// *Istituto Superiore di Sanità and INFN Roma, gruppo collegato Sanità, Italy
//  Viale Regina Elena 299, 00161 Roma (Italy)
//  tel (39) 06 49902246
//  fax (39) 06 49387075
//
// more information:
// http://g4advancedexamples.lngs.infn.it/Examples/medical-linac
//
//*******************************************************//


#include "ML2WorldConstruction.hh"
#include "G4SystemOfUnits.hh"

  G4double density,      // density
    a,                   // atomic mass
    z;                   // atomic number
  G4String name,         // name
    symbol;              // symbol
  G4int ncomponents,     // n components
    iz,                  // number of protons
    in;                  // number of nuceons
  G4double abundance,    // abundance
    temperature,         // temperature
    pressure;            // pressure




CML2WorldConstruction::CML2WorldConstruction():acceleratorEnv(0),phantomEnv(0),PVWorld(0),phaseSpace(0),backScatteredPlane(0)
{
    phantomEnv = CML2PhantomConstruction::GetInstance();
    acceleratorEnv = CML2AcceleratorConstruction::GetInstance();
//    wallsEnv=CML2Walls::GetInstance();
    bWorldCreated = false;
    bOnlyVisio = 0;
}

CML2WorldConstruction::~CML2WorldConstruction(void)
{
    delete phaseSpace;
    delete backScatteredPlane;
}

CML2WorldConstruction* CML2WorldConstruction::instance = 0;

CML2WorldConstruction* CML2WorldConstruction::GetInstance()
{
    if (instance == 0)
    {
        instance = new CML2WorldConstruction();
    }
    return instance;
}

G4VPhysicalVolume* CML2WorldConstruction::Construct()
{
    return PVWorld;
}



bool CML2WorldConstruction::create(SInputData *inputData, bool bOV)
{
    // create the world box
    bOnlyVisio = bOV;
    G4double halfSize = 6000.*mm;
    G4Material *Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Box *worldB = new G4Box("worldG", halfSize, halfSize, halfSize);
    G4LogicalVolume *worldLV = new G4LogicalVolume(worldB, Vacuum, "worldL", 0, 0, 0);
    G4VisAttributes* simpleWorldVisAtt = new G4VisAttributes(G4Colour::Red());
    simpleWorldVisAtt -> SetVisibility(true);
    worldLV -> SetVisAttributes(simpleWorldVisAtt);
    PVWorld = new G4PVPlacement(0,  G4ThreeVector(0.,0.,0.), "worldPV", worldLV, 0, false, 0);
    
    
    // create the phantom-world box
    //bOnlyVisio = bOV;
    G4double halfSizePHANX = 150.*mm;
    G4double halfSizePHANY = 3000.*mm;
    G4double halfSizePHANZ = 3000.*mm;
    //G4Material *Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4Box *worldB1 = new G4Box("worldG", halfSizePHANX, halfSizePHANY, halfSizePHANZ);
    G4LogicalVolume *worldLV1 = new G4LogicalVolume(worldB1, Vacuum, "worldL1", 0, 0, 0);
    G4VisAttributes* simpleWorld1VisAtt = new G4VisAttributes(G4Colour::Blue());
    simpleWorldVisAtt -> SetVisibility(true);
    worldLV1 -> SetVisAttributes(simpleWorld1VisAtt);
    PVWorld1 = new G4PVPlacement(0,  G4ThreeVector(-3200.,0.,0.), "world1PV", worldLV1, PVWorld, false, 0);/////////////sdfsdfsddf,jkjdfbgkdjfbg
    
    
    
    
    
    
    //concrete
  //G4Element* H = new G4Element
  //  (name="Hydrogen",symbol="H" , z= 1., a=1.00794*g/mole);
  //G4Element* Ca = new G4Element
  //  (name="Calcium",symbol="Ca" , z= 20., a=40.078*g/mole);
  G4Material* concrete = new G4Material
    (name="Concrete", density=2.3*g/cm3, ncomponents=6);
  G4Material *Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");
  G4Material *O = G4NistManager::Instance()->FindOrBuildMaterial("G4_O");
  G4Material *H = G4NistManager::Instance()->FindOrBuildMaterial("G4_H");
  G4Material *Ca = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ca");
  G4Material *Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  G4Material *Fe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe");
  concrete->AddMaterial(Si, 0.227915);
  concrete->AddMaterial(O, 0.60541);
  concrete->AddMaterial(H, 0.09972);
  concrete->AddMaterial(Ca, 0.04986);
  concrete->AddMaterial(Al, 0.014245);
  concrete->AddMaterial(Fe, 0.00285);


    //Creamos las paredes las cuales debemos antes ubicar en el
    //Plano para asi conocer y verificar la geometria y overlaps.

    //Pared 1
    G4double halfSizeX1 = 100.*mm;
    G4double halfSizeY1 = 3000.*mm;
    G4double halfSizeZ1 = 3000.*mm;
    G4Box *Pared1B = new G4Box("Pared1G", halfSizeX1, halfSizeY1, halfSizeZ1);
    G4LogicalVolume *Pared1LV = new G4LogicalVolume(Pared1B, concrete, "Pared1L", 0, 0, 0);
    G4VisAttributes* simplePared1VisAtt = new G4VisAttributes(G4Colour::Blue());
    simplePared1VisAtt -> SetVisibility(true);
    Pared1LV -> SetVisAttributes(simplePared1VisAtt);
    PVPared1 = new G4PVPlacement(0,  G4ThreeVector(2900.,0.,0.), "Pared1PV", Pared1LV, PVWorld, false, 0);

    //Pared Pared 2
    G4double halfSizeX2 = 2800.*mm;
    G4double halfSizeY2 = 3000.*mm;
    G4double halfSizeZ2 = 100.*mm;
    G4Box *Pared2B = new G4Box("Pared2G", halfSizeX2, halfSizeY2, halfSizeZ2);
    G4LogicalVolume *Pared2LV = new G4LogicalVolume(Pared2B, concrete, "Pared2L", 0, 0, 0);
    G4VisAttributes* simplePared2VisAtt = new G4VisAttributes(G4Colour::Blue());
    simplePared2VisAtt -> SetVisibility(true);
    Pared2LV -> SetVisAttributes(simplePared2VisAtt);
    PVPared2 = new G4PVPlacement(0,  G4ThreeVector(0.,0.,2900.), "Pared2PV", Pared2LV, PVWorld, false, 0);
    
    //Pared 3
    G4double halfSizeX3 = 100.*mm;
    G4double halfSizeY3 = 3000.*mm;
    G4double halfSizeZ3 = 3000.*mm;
    G4Box *Pared3B = new G4Box("Pared3G", halfSizeX3, halfSizeY3, halfSizeZ3);
    G4LogicalVolume *Pared3LV = new G4LogicalVolume(Pared3B, concrete, "Pared3L", 0, 0, 0);
    G4VisAttributes* simplePared3VisAtt = new G4VisAttributes(G4Colour::Blue());
    simplePared3VisAtt -> SetVisibility(true);
    Pared3LV -> SetVisAttributes(simplePared3VisAtt);
   PVPared3 = new G4PVPlacement(0,  G4ThreeVector(-2900.,0.,0.), "Pared3PV", Pared3LV, PVWorld, false, 0);
    
    //Pared 4
    G4double halfSizeX4 = 2800.*mm;
    G4double halfSizeY4 = 3000.*mm;
    G4double halfSizeZ4 = 100.*mm;
    G4Box *Pared4B = new G4Box("Pared4G", halfSizeX4, halfSizeY4, halfSizeZ4);
    G4LogicalVolume *Pared4LV = new G4LogicalVolume(Pared4B, concrete, "Pared4L", 0, 0, 0);
    G4VisAttributes* simplePared4VisAtt = new G4VisAttributes(G4Colour::Blue());
    simplePared4VisAtt -> SetVisibility(true);
    Pared4LV -> SetVisAttributes(simplePared4VisAtt);
   PVPared4 = new G4PVPlacement(0,  G4ThreeVector(0.,0.,-2900.), "Pared4PV", Pared4LV, PVWorld, false, 0);
    

   
    // create the accelerator-world box
    if (!acceleratorEnv -> Construct(PVWorld, bOV))
    {
        G4cout << "\n\n The macro file '" << inputData->generalData.StartFileInputData <<
        		"' refers to a not defined accelerator.\n" << acceleratorEnv->getAcceleratorName() <<
				"\n\nSTOP\n\n" << G4endl;
        return false;
    }

    // create the phantom-world box ///////////////AQUIUI!!!!!!!!!!111
    if ( !phantomEnv->Construct(PVWorld1,
                        inputData->voxelSegmentation.nX,
                        inputData->voxelSegmentation.nY,
                        inputData->voxelSegmentation.nZ,
			bOV) )
    {
        G4cout << "\n\n The macro file '" << inputData->generalData.StartFileInputData <<
        		"' refers to a not defined phantom.\n" << phantomEnv->getPhantomName() <<
				"\n\nSTOP\n\n" << G4endl;
        return false;
    }

    // if the bSavePhaseSpace flag is true create a phase plane
/*   
 if (inputData -> generalData.bSavePhaseSpace)
    {
        phaseSpace = new CML2PhaseSpaces();
        if (inputData -> generalData.bForcePhaseSpaceBeforeJaws)
        {
        	inputData -> generalData.centrePhaseSpace.setZ(acceleratorEnv->getZ_Value_PhaseSpaceBeforeJaws());
        }

        phaseSpace -> createPlane(idSD_PhaseSpace,
        		inputData->generalData.max_N_particles_in_PhSp_File,
        		inputData->generalData.seed,
				inputData->generalData.nMaxParticlesInRamPlanePhaseSpace,
				acceleratorEnv->getPhysicalVolume(), "PhSp",
				inputData->generalData.PhaseSpaceOutFile,
				inputData->generalData.bSavePhaseSpace,
				inputData->generalData.bStopAtPhaseSpace,
				inputData->generalData.centrePhaseSpace,
				inputData->generalData.halfSizePhaseSpace,
				&inputData->primaryParticleData,
				acceleratorEnv->getAcceleratorIsoCentre()); // phase space plane, yellow
    }

    // create a killer plane to destroy the particles back scattered from the target
    backScatteredPlane = new CML2PhaseSpaces();
    backScatteredPlane -> createPlane(acceleratorEnv->getPhysicalVolume(),
    		"killerPlane", G4ThreeVector(0, 0, -50*mm), G4ThreeVector(200*mm, 200*mm, 1*mm)); // killer plane, cyan
*/
    bWorldCreated = true;
    return bWorldCreated;
}
void CML2WorldConstruction::checkVolumeOverlap()
{
    // loop inside all the daughters volumes
	G4cout<< G4endl;
    //        bool bCheckOverlap;
    //        bCheckOverlap=false;

    int nSubWorlds, nSubWorlds2;
    for (int i=0; i<(int) PVWorld->GetLogicalVolume()->GetNoDaughters(); i++)
    {
        PVWorld->GetLogicalVolume()->GetDaughter(i)->CheckOverlaps();
        nSubWorlds=(int) PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetNoDaughters();
        for (int j=0; j<nSubWorlds; j++)
        {
            PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->CheckOverlaps();
            nSubWorlds2=(int) PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->GetLogicalVolume()->GetNoDaughters();
            for (int k=0; k<nSubWorlds2; k++)
            {
                PVWorld->GetLogicalVolume()->GetDaughter(i)->GetLogicalVolume()->GetDaughter(j)->GetLogicalVolume()->GetDaughter(k)->CheckOverlaps();
            }
        }
    }
    G4cout<< G4endl;
}
bool CML2WorldConstruction::newGeometry()
{

    G4bool bNewRotation = false;
    G4bool bNewCentre = false;
    G4bool bNewGeometry = false;
    bNewCentre = phantomEnv -> applyNewCentre();
    G4RotationMatrix *rmInv = acceleratorEnv -> rotateAccelerator();
    if (rmInv!=0)
    {
        CML2PrimaryGenerationAction::GetInstance()->setRotation(rmInv);
        bNewRotation = true;
    }
    if (bNewRotation || bNewCentre)
    {
    	bNewGeometry = true;
    }

    return bNewGeometry;
}

