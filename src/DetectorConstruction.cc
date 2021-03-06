#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <boost/algorithm/string.hpp>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnionSolid.hh"

using namespace std;

//

DetectorConstruction::DetectorConstruction(G4int ver) :
                version_(ver), hcalThickness(0){

//	if (version_ == 2)
//		buildTracker();

	buildECal();
//	buildHCal();
      
	DefineMaterials();
	m_detectorMessenger = new DetectorMessenger(this);
	UpdateCalorSize();
        SetMagField(0);

}

void DetectorConstruction::buildTracker(){
	initLayer(1);

	std::vector<std::pair <G4double,std::string> > iEleL;
        iEleL.push_back(make_pair(150*cm,"G4_Galactic"));

	for (int i = 0; i < 6; i ++){
		iEleL.push_back(make_pair(.7*mm,"Si"));
		iEleL.push_back(make_pair(99.3*mm,"G4_Galactic"));
	}


	iEleL.push_back(make_pair(.7*mm,"Si"));
	iEleL.push_back(make_pair(6.8*mm,"G4_Galactic"));

	iEleL.push_back(make_pair(.3504*mm,"W"));
	iEleL.push_back(make_pair(6.8*mm,"G4_Galactic"));


	for (int i = 0; i < 3; i ++){
		iEleL.push_back(make_pair(.7*mm,"Si"));
		if (i < 3){
			iEleL.push_back(make_pair(15*mm,"G4_Galactic"));

		}
	}

	iEleL.push_back(make_pair(.7*mm,"Si"));
	iEleL.push_back(make_pair(44.8*mm,"G4_Galactic"));
	iEleL.push_back(make_pair(.7*mm,"Si"));
	iEleL.push_back(make_pair(84.3*mm,"G4_Galactic"));
	iEleL.push_back(make_pair(.7*mm,"Si"));
	iEleL.push_back(make_pair(14.3*mm,"G4_Galactic"));


	m_caloStruct.push_back( SamplingSection(iEleL) );
}

void DetectorConstruction::buildECal(){

        std::vector<std::pair <G4double,std::string> > iEleL;

		double wThick = 2.25, cuThick = 5;//5;

		iEleL.clear();
		iEleL.push_back(make_pair(2*mm,"Air"));
		iEleL.push_back(make_pair(2*mm,"PCB"));
		iEleL.push_back(make_pair(0.5*mm,"CFMix"));
		iEleL.push_back(make_pair(0.5*mm,"Si"));
		iEleL.push_back(make_pair(0.5*mm,"CFMix"));
		iEleL.push_back(make_pair(cuThick*mm,"Cu"));
		iEleL.push_back(make_pair(wThick,"W"));

		int Nmodule=40;
		for(unsigned i=0; i<Nmodule; i++) {
				m_caloStruct.push_back( SamplingSection(iEleL) );
		}
}
/*
void DetectorConstruction::buildECal(){

	std::vector<std::pair <G4double,std::string> > initBlock,iEleL,iEleR;

	//atleast 6 eles in first layer of ECAL
	initBlock.push_back(make_pair(1*mm, "Al"));
	initBlock.push_back(make_pair(1*mm, "Al"));
	initBlock.push_back(make_pair(1.5*mm,"PCB"));
	initBlock.push_back(make_pair(0.05*mm,"CFMix"));
	initBlock.push_back(make_pair(0.5*mm,"Si"));
	initBlock.push_back(make_pair(0.2*mm,"CFMix"));


	///Single element in block
	iEleL.push_back(make_pair(0.0*mm,"W"));

	iEleL.push_back(make_pair(1*mm,"Air"));
	iEleL.push_back(make_pair(1.5*mm,"PCB"));
	iEleL.push_back(make_pair(0.05*mm,"CFMix"));
	iEleL.push_back(make_pair(0.5*mm,"Si"));
	iEleL.push_back(make_pair(0.2*mm,"CFMix"));

	iEleR.push_back(make_pair(0.0*mm,"W"));

	iEleR.push_back(make_pair(0.35*mm,"C"));
	iEleR.push_back(make_pair(1.3*mm,"Air"));
	iEleR.push_back(make_pair(0.35*mm,"C"));

	iEleR.push_back(make_pair(0.0*mm,"W"));

	iEleR.push_back(make_pair(0.2*mm,"CFMix"));
	iEleR.push_back(make_pair(0.5*mm,"Si"));
	iEleR.push_back(make_pair(0.05*mm,"CFMix"));
	iEleR.push_back(make_pair(1.5*mm,"PCB"));
	iEleR.push_back(make_pair(1*mm,"Air"));

	unsigned Nmodule=3;
	for(unsigned i=0; i<Nmodule; i++) {
		switch(i){
			case 0: m_caloStruct.push_back( SamplingSection(initBlock) ); break;
			case 1: iEleL[0].first = 0.75 *mm;iEleR[0].first = 0.75 *mm;iEleR[4].first = 0.75 *mm; break;
			case 2: iEleL[0].first = 2.25 *mm;iEleR[0].first = 1.5 *mm;iEleR[4].first = 1.5 *mm; break;
		}
		if (i > 0)
			m_caloStruct.push_back( SamplingSection(iEleL) );
			m_caloStruct.push_back( SamplingSection(iEleR) );

	}


	iEleL[0].first = 3.5 *mm;iEleR[0].first = 1.75 *mm;iEleR[4].first = 1.75 *mm;
	Nmodule=9;
	for(unsigned i=0; i<Nmodule; i++) {
		m_caloStruct.push_back( SamplingSection(iEleL) );
		m_caloStruct.push_back( SamplingSection(iEleR) );

	}

	iEleL[0].first = 7 *mm;iEleR[0].first = 3.5 *mm;iEleR[4].first = 3.5 *mm;
	Nmodule=5;
	for(unsigned i=0; i<Nmodule; i++) {
		m_caloStruct.push_back( SamplingSection(iEleL) );
		m_caloStruct.push_back( SamplingSection(iEleR) );

	}


}
*/
void DetectorConstruction::buildHCal(){
	std::vector<std::pair <G4double,std::string> > iEleL;

	iEleL.push_back(make_pair(50*mm,"Steel"));
	iEleL.push_back(make_pair(2*mm,"Air"));
	iEleL.push_back(make_pair(6*mm,"Scintillator"));
	iEleL.push_back(make_pair(2*mm,"Air"));
	unsigned Nmodule=50;

	for(unsigned i=0; i<Nmodule; i++) {
		m_caloStruct.push_back( SamplingSection(iEleL) );
	}

}

DetectorConstruction::~DetectorConstruction() {
	delete m_detectorMessenger;
}

void DetectorConstruction::DefineMaterials() {
	G4NistManager* nistManager = G4NistManager::Instance();

	m_materials["Al"] = nistManager->FindOrBuildMaterial("G4_Al", false);
	m_dEdx["Al"] = 0.4358;
	m_materials["W"] = nistManager->FindOrBuildMaterial("G4_W", false);
	m_dEdx["W"] = 2.210;
	m_materials["Pb"] = nistManager->FindOrBuildMaterial("G4_Pb", false);
	m_dEdx["Pb"] = 1.274;
	m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu", false);
	m_dEdx["Cu"] = 1.257;
	m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si", false);
	m_dEdx["Si"] = 0.3876;
	m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn", false);
	m_dEdx["Zn"] = 1.007;
	m_materials["Air"] = nistManager->FindOrBuildMaterial("G4_AIR", false);
	m_dEdx["Air"] = 0;
	m_materials["G4_Galactic"] = nistManager->FindOrBuildMaterial("G4_Galactic",false);
	m_dEdx["G4_Galactic"] = 0;
	m_materials["Fe"] = nistManager->FindOrBuildMaterial("G4_Fe", false);
	m_dEdx["Fe"] = 1.143;
	m_materials["Mn"] = nistManager->FindOrBuildMaterial("G4_Mn", false);
	m_dEdx["Mn"] = 1.062;
	m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C", false);
	m_dEdx["C"] = 0.3952;
	m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H", false);
	m_dEdx["H"] = 0;
	m_materials["Cl"] = nistManager->FindOrBuildMaterial("G4_Cl", false);
	m_dEdx["Cl"] = 0;
	m_materials["Cr"] = nistManager->FindOrBuildMaterial("G4_Cr", false);
	m_dEdx["Cr"] = 1.046;
	m_materials["Ni"] = nistManager->FindOrBuildMaterial("G4_Ni", false);
	m_dEdx["Ni"] = 1.307;
	m_materials["O"] = nistManager->FindOrBuildMaterial("G4_O", false);
	m_materials["Br"] = nistManager->FindOrBuildMaterial("G4_Br", false);
	m_materials["glass"] = nistManager->FindOrBuildMaterial("G4_GLASS_PLATE", false);

	m_materials["PCB"] = new G4Material("FR4", 5.68 * g / cm3, 2);
	m_materials["PCB"]->AddMaterial(m_materials["Cu"], 0.5);
	m_materials["PCB"]->AddMaterial(m_materials["glass"], 0.5);
	m_dEdx["PCB"] = 0;

	m_materials["Brass"] = new G4Material("Brass", 8.53 * g / cm3, 2);
	m_materials["Brass"]->AddMaterial(m_materials["Cu"], 70 * perCent);
	m_materials["Brass"]->AddMaterial(m_materials["Zn"], 30 * perCent);
	m_dEdx["Brass"] = 0.7 * m_dEdx["Cu"] + 0.3 * m_dEdx["Zn"];
	m_materials["Steel"] = new G4Material("Steel", 7.87 * g / cm3, 3);
	m_materials["Steel"]->AddMaterial(m_materials["Fe"], 0.9843);
	m_materials["Steel"]->AddMaterial(m_materials["Mn"], 0.014);
	m_materials["Steel"]->AddMaterial(m_materials["C"], 0.0017);
	m_dEdx["Steel"] = 0.9843 * m_dEdx["Fe"] + 0.014 * m_dEdx["Mn"]
			+ 0.0017 * m_dEdx["C"];
	m_materials["SSteel"] = new G4Material("SSteel", 8.02 * g / cm3, 4);
	m_materials["SSteel"]->AddMaterial(m_materials["Fe"], 0.70);
	m_materials["SSteel"]->AddMaterial(m_materials["Mn"], 0.01);
	m_materials["SSteel"]->AddMaterial(m_materials["Cr"], 0.19);
	m_materials["SSteel"]->AddMaterial(m_materials["Ni"], 0.10);
	m_dEdx["SSteel"] = 0.7 * m_dEdx["Fe"] + 0.01 * m_dEdx["Mn"]
			+ 0.19 * m_dEdx["Cr"] + 0.1 * m_dEdx["Ni"];

	m_materials["Scintillator"] = new G4Material("Scintillator",
			1.032 * g / cm3, 2);
	m_materials["Scintillator"]->AddMaterial(m_materials["C"],
			91.512109 * perCent);
	m_materials["Scintillator"]->AddMaterial(m_materials["H"],
			8.4878906 * perCent);
	m_dEdx["Scintillator"] = m_dEdx["C"];

	G4cout << m_materials["Scintillator"] << G4endl;
	m_materials["Polystyrole"] = new G4Material("Polystyrole", 1.065 * g / cm3,
			2);
	m_materials["Polystyrole"]->AddMaterial(m_materials["H"], 50 * perCent);
	m_materials["Polystyrole"]->AddMaterial(m_materials["C"], 50 * perCent);
	m_dEdx["Polystyrole"] = 0.5 * m_dEdx["C"];

	m_materials["PVC"] = new G4Material("PVC", 1.350 * g / cm3, 3);
	m_materials["PVC"]->AddMaterial(m_materials["H"], 50 * perCent);
	m_materials["PVC"]->AddMaterial(m_materials["C"], 33.33 * perCent);
	m_materials["PVC"]->AddMaterial(m_materials["Cl"], 16.67 * perCent);
	m_dEdx["PVC"] = 0.33 * m_dEdx["C"];

	m_materials["CFMix"] = new G4Material("CFMix", 0.845 * g / cm3, 3);
	m_materials["CFMix"]->AddMaterial(m_materials["C"], 0.84491305);
	m_materials["CFMix"]->AddMaterial(m_materials["H"], 0.042542086);
	m_materials["CFMix"]->AddMaterial(m_materials["O"], 0.11254487);
	m_dEdx["CFMix"] = 0;

	m_materials["Foam"] = new G4Material("Foam", 0.0999 * g / cm3, 2);
	m_materials["Foam"]->AddMaterial(m_materials["C"], 0.856);
	m_materials["Foam"]->AddMaterial(m_materials["H"], 0.144);
	m_dEdx["Foam"] = 1.749 * 0.856 * 0.0999 / 10.;

	m_materials["WCu"] = new G4Material("WCu", 14.979 * g / cm3, 2);
	m_materials["WCu"]->AddMaterial(m_materials["W"], 75 * perCent);
	m_materials["WCu"]->AddMaterial(m_materials["Cu"], 25 * perCent);
	m_dEdx["WCu"] = 0.75 * m_dEdx["W"] + 0.25 * m_dEdx["Cu"];

	m_materials["NeutMod"] = new G4Material("NeutMod", 0.950 * g / cm3, 2);
	m_materials["NeutMod"]->AddMaterial(m_materials["C"], 0.85628);
	m_materials["NeutMod"]->AddMaterial(m_materials["H"], 0.14372);
	m_dEdx["NeutMod"] = 1.749 * 0.86 * 0.950 / 10.;

}

void DetectorConstruction::UpdateCalorSize() {
	m_CalorSizeZ = 0;
	m_ECALSizeZ = 0;
	for (size_t i = 0; i < m_caloStruct.size(); i++) {
		m_CalorSizeZ = m_CalorSizeZ + m_caloStruct[i].Total_thick;
		if (i >= initLayer_ && m_caloStruct[i].n_elements >= 6){
			m_ECALSizeZ += m_caloStruct[i].Total_thick;
		}
	//G4cout << "m_CalorSizeZ = " << m_CalorSizeZ << G4endl;
	//G4cout << "m_ECALSizeZ = " << m_ECALSizeZ << G4endl;
	}
	m_nSectors = 1;



	m_ECALSizeXY = 6*85*1.25;//700; //use full length for making hexagon map
	m_HCALSizeXY = 3100; //use full length for making square map

	std::cout << "The calorimeter size is " << m_ECALSizeXY << "the width is " << m_sectorWidth << std::endl;
	m_sectorWidth = m_ECALSizeXY;

	for (size_t i = 0; i < m_caloStruct.size(); i++)
		m_caloStruct[i].setNumberOfSectors(m_nSectors);

	m_WorldSizeZ = m_CalorSizeZ + 2.5 * cm;
        m_WorldSizeXY = m_ECALSizeXY * 1.1;
//Uncomment for HCAL        m_WorldSizeXY = m_HCALSizeXY * 1.1;


	G4cout << "[DetectorConstruction][UpdateECALSize] Z x XY = "
			<< m_CalorSizeZ << " x " << m_ECALSizeXY << " mm "
			<< ", nsectors = " << m_nSectors << G4endl;
	G4cout << "[DetectorConstruction][UpdateHCALSize] Z x XY = "
			<< m_CalorSizeZ << " x " << m_HCALSizeXY << " mm "
			<< ", nsectors = " << m_nSectors << G4endl;

}

G4VPhysicalVolume* DetectorConstruction::Construct() {
	//clean old geometry
	G4GeometryManager::GetInstance()->OpenGeometry();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();

	//world
	G4double expHall_z = 20 * m;
	G4double expHall_x = 10* m;
	G4double expHall_y = 10 * m;

	G4Box* experimentalHall_box = new G4Box("expHall_box", expHall_x, expHall_y,
			expHall_z);

	G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(
			experimentalHall_box, m_materials["Air"], "expHall_log");
	G4VPhysicalVolume* experimentalHall_phys = new G4PVPlacement(0, // no rotation
			G4ThreeVector(0., 0., 0.), // translation position
			experimentalHall_log,    // its logical volume
			"expHall",               // its name
			0,                       // its mother volume
			false,                   // no boolean operations
			0);                      // its copy number

	//detector's World
	G4double pos_x = 0.;
	G4double pos_y = 0.;
	G4double pos_z = 0.;

	m_solidWorld = new G4Box("Wbox", m_WorldSizeXY / 2, m_WorldSizeXY / 2,
			m_WorldSizeZ / 2);

	cout << "m_WorldSizeXY = " << m_WorldSizeXY << ", m_WorldSizeZ = "
		<< m_WorldSizeZ << endl;

	m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"],
			"Wlog");
	m_physWorld = new G4PVPlacement(0, G4ThreeVector(pos_x, pos_y, pos_z),
			m_logicWorld, "Wphys", experimentalHall_log, false, 0);

	for (unsigned iS(0); iS < m_nSectors; ++iS) {
		G4double minL = m_sectorWidth * iS;
		buildSectorStack(iS, minL, m_sectorWidth);
		if (m_nSectors > 1)
			fillInterSectorSpace(iS, minL + m_sectorWidth, m_interSectorWidth);
	}
	// Visualization attributes
	//
	m_logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

	return experimentalHall_phys;
}

void DetectorConstruction::buildSectorStack(const unsigned sectorNum,
		const G4double & minL, const G4double & width) {

	//build the stack
	G4double zOffset(-m_CalorSizeZ / 2), zOverburden(0.);
	char nameBuf[10];
	G4VSolid *solid;

	G4double totalLengthX0 = 0;
	G4double totalLengthL0 = 0;

	for (size_t i = 0; i < m_caloStruct.size(); i++) {
		const unsigned nEle = m_caloStruct[i].n_elements;
		//index for counting Si sensitive layers
		unsigned idx = 0;

		for (unsigned ie(0); ie < nEle; ++ie) {
			std::string eleName = m_caloStruct[i].ele_name[ie];

			if (m_nSectors == 1)
				sprintf(nameBuf, "%s%d", eleName.c_str(), int(i + 1));

			else
				sprintf(nameBuf, "%s%d_%d", eleName.c_str(), int(sectorNum),
						int(i + 1));

			if (eleName == "Si") {
				if (m_nSectors == 1)
					sprintf(nameBuf, "Si%d_%d", int(i + 1), idx);

				else
					sprintf(nameBuf, "Si%d_%d_%d", int(sectorNum), int(i + 1),
							idx);
				idx++;
			}
			std::string baseName(nameBuf);
			G4double thick = m_caloStruct[i].sublayer_thick[ie];
			//

			if (thick > 0) {
#if 0
				cout << "solid = constructSolid("<<baseName
					<<",thick="<<thick
					<<",zOffset+zOverburden="<<zOffset+zOverburden
					<<",width="<<width<<");"<<endl;
#endif
				solid = constructSolid(baseName, thick, zOffset + zOverburden, minL, width, i, nEle);
				G4LogicalVolume *logi = new G4LogicalVolume(solid, m_materials[eleName], baseName + "log");
				m_caloStruct[i].sublayer_X0[ie] = m_materials[eleName]->GetRadlen();
				m_caloStruct[i].sublayer_dEdx[ie] = m_dEdx[eleName];
				m_caloStruct[i].sublayer_L0[ie] = m_materials[eleName]->GetNuclearInterLength();
				if (sectorNum == 0 || sectorNum == m_nSectors - 1) {
					G4cout << "************ " << eleName;
					if (m_nSectors > 1)
						G4cout << " sector " << sectorNum;
					G4cout << " layer " << i << " dEdx="
							<< m_caloStruct[i].sublayer_dEdx[ie] << " X0="
							<< m_caloStruct[i].sublayer_X0[ie] << " L0="
							<< m_caloStruct[i].sublayer_L0[ie] << " zpos="
							<< zOffset + zOverburden << "mm w="
							<< m_caloStruct[i].sublayer_thick[ie] << "mm";

					totalLengthX0 += m_caloStruct[i].sublayer_thick[ie]
							/ m_caloStruct[i].sublayer_X0[ie];
					G4cout << " TotX0=" << totalLengthX0;
					totalLengthL0 += m_caloStruct[i].sublayer_thick[ie]
							/ m_caloStruct[i].sublayer_L0[ie];
					G4cout << " TotLambda=" << totalLengthL0;
					G4cout << " Number of elements=" << nEle << G4endl;
				}

				if (m_caloStruct[i].isSensitiveElement(ie))
					m_logicSi.push_back(logi);

				G4double xpvpos = -m_ECALSizeXY / 2. + minL + width / 2;
#if 0
				cout << "m_caloStruct[i].ele_vol[nEle*sectorNum+ie]=new G4PVPlacement(0, G4ThreeVector(xpvpos="<<xpvpos
				<< ",0.,zOffset+zOverburden+thick/2="<<zOffset+zOverburden+thick/2
				<< "), logi,"
				<< baseName+"phys, m_logicWorld, false, 0);" << endl;
#endif

 				G4RotationMatrix* rot = new G4RotationMatrix(0,0,1);


				if ( (i == initLayer_ - 1) || nEle < 6){

				G4double extraOffset = 0 * mm;
				if (i != initLayer_ - 1) extraOffset = 20 *mm;
				m_caloStruct[i].sublayer_vol[nEle * sectorNum + ie] =
								new G4PVPlacement(0,
												G4ThreeVector(xpvpos, 0.,
																zOffset + zOverburden + thick / 2 + extraOffset),
												logi, baseName + "phys", m_logicWorld, false,
												0);
				}
				else{
			 	G4double extraOffset = 10 * mm;
				m_caloStruct[i].sublayer_vol[nEle * sectorNum + ie] =
								new G4PVPlacement(G4Transform3D(*rot,
												G4ThreeVector(xpvpos, 0.,
																zOffset + zOverburden + extraOffset )),
												logi, baseName + "phys", m_logicWorld, false,
												0);
				}


				G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(
						m_caloStruct[i].g4Colour(ie));
				simpleBoxVisAtt->SetVisibility(true);
				logi->SetVisAttributes(simpleBoxVisAtt);
				zOverburden = zOverburden + thick;
				//for sensitive volumes
				//add region to be able to set specific cuts for it
				//just for Si
				if (eleName == "Si") {
					unsigned nlogicsi = m_logicSi.size();
					G4Region* aRegion = new G4Region(baseName + "Reg");
					m_logicSi[nlogicsi - 1]->SetRegion(aRegion);
					aRegion->AddRootLogicalVolume(m_logicSi[nlogicsi - 1]);
				}
			}

		}	  //loop on elements
	}	  //loop on layers

}	  //buildstack

void DetectorConstruction::fillInterSectorSpace(const unsigned sectorNum,
		const G4double & minL, const G4double & width) {
	//build the stack
	G4double zOffset(-m_CalorSizeZ / 2), zOverburden(0.);
	char nameBuf[10];
	G4VSolid *solid;

	for (size_t i = 0; i < m_caloStruct.size(); i++) {
		const unsigned nEle = m_caloStruct[i].n_elements;
		for (unsigned ie(0); ie < nEle; ++ie) {

			std::string eleName = m_caloStruct[i].ele_name[ie];
			G4double thick = m_caloStruct[i].sublayer_thick[ie];

			eleName = "CFMix";
			sprintf(nameBuf, "%s%d_%d", eleName.c_str(), int(sectorNum),
					int(i + 1));
			std::string baseName(nameBuf);
			if (thick > 0) {
				solid = constructSolid(baseName, thick, zOffset + zOverburden,
						minL, width, i, nEle);
				std::cout << "Constructing layer " << baseName << "phys" << std::endl;
				G4LogicalVolume *logi = new G4LogicalVolume(solid,
						m_materials[eleName], baseName + "log");
				G4double xpvpos = -m_ECALSizeXY / 2. + minL + width / 2;


				G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(
						G4Colour::Magenta);
				simpleBoxVisAtt->SetVisibility(true);
				simpleBoxVisAtt->SetForceSolid(true);
				logi->SetVisAttributes(simpleBoxVisAtt);
				zOverburden = zOverburden + thick;
			}
		}	  //loop on elements
	}	  //loop on layers
}	  //fill intersector space

void DetectorConstruction::SetMagField(G4double fieldValue) {

	if (fieldValue <= 0)
		return;

	//apply a global uniform magnetic field along Z axis
	G4FieldManager* fieldMgr =
			G4TransportationManager::GetTransportationManager()->GetFieldManager();
	if (m_magField)
		delete m_magField;                //delete the existing magn field
	m_magField = new G4UniformMagField(G4ThreeVector(0., 0., fieldValue));
	fieldMgr->SetDetectorField(m_magField);
	fieldMgr->CreateChordFinder(m_magField);
	fieldMgr->SetDetectorField(m_magField);
}

void DetectorConstruction::SetMagField(char *fileName, G4double zOffset)
{
   // G4cout << "!!@! ENTERING SETMAGFIELD !!@!" << G4endl;
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    if(p_magField) delete p_magField; // Delete the existing magnetic field
    p_magField = new PurgMagTabulatedField3D(fileName, zOffset);
    for (int i = 0; i < 350;i++){
    	 G4double disp = i * 10 * mm;
    	 double magField[3] = {0,0,0};
    	 //G4double point[4] = {0,0,-(0.5*GetCalorSizeZ() - 607.5*mm - disp),0};
    	 const G4double point[4] = {0,0,-GetCalorSizeZ()/2.0 + disp,0};
    	 p_magField->GetFieldValue(point,magField);
    	 //G4cout << "The magnetic field(z = " << -(0.5*GetCalorSizeZ() - 607.5*mm - disp) << "*mm)" << "= "  << magField[2] << G4endl;
    	 G4cout << "The magnetic field(z = " << -GetCalorSizeZ()/2.0 + disp << "*mm)" << "= "  << magField[1]/tesla << G4endl;

    }
    fieldMgr->SetDetectorField(p_magField);
    fieldMgr->CreateChordFinder(p_magField);
    fieldMgr->SetDetectorField(p_magField);
    //G4cout << "!!@! EXITING SETMAGFIELD !!@!" << G4endl;

}


G4VSolid *DetectorConstruction::constructSolid(std::string baseName,
		G4double thick, G4double zpos, const G4double & minL,
		const G4double & width, size_t which_ele, const unsigned nEle) {
		G4VSolid *solid;
	if (which_ele == initLayer_ - 1) {
                solid = new G4Box(baseName + "box", m_HCALSizeXY / 2, m_HCALSizeXY / 2, thick / 2);
		return solid;
	}
	if (nEle >= 6){ // this means that it's ECAL!
		G4double hexaRad = 85.;
		G4double a[2] = {0,thick},b[2] ={0,0},c[2] = {hexaRad,hexaRad};


		G4VSolid* sHexa = new G4Polyhedra(baseName + "box",
						0, 2 * pi,
						6, 2,
						a,b,c);

		G4ThreeVector trans = G4ThreeVector(0.,2*hexaRad,0);
		G4ThreeVector zAxis(0,0,1);
		G4RotationMatrix* rot = new G4RotationMatrix(0,0,0);

		G4UnionSolid* sUnion = new G4UnionSolid(baseName + "box",
				sHexa,
				sHexa,
				rot,
				trans);
		for (int i = 1; i < 6; i ++){
			G4ThreeVector trans = G4ThreeVector(0.,2*hexaRad,0);

			sUnion = new G4UnionSolid(baseName + "box",
				sUnion,
				sHexa,
				rot,
				trans.rotate(i*3.14/3,zAxis ));
		}

		solid = sUnion;
	}
	else{
		solid = new G4Box(baseName + "box", m_HCALSizeXY / 2, m_HCALSizeXY / 2, thick / 2);
/*		if (hcalThickness/60 < 20){
			G4double hexaRad = 85.;
			G4double hcalXYSide = hexaRad*6,hcalLen=m_ECALSizeZ,constShift = -15;
			G4ThreeVector trans = G4ThreeVector(hcalXYSide/2 + hcalThickness + thick/2.0+10,0,-hcalLen/2.0 -hcalThickness - thick/2.0 + constShift);


			G4VSolid *solid_1 = new G4Box(baseName + "box",hcalLen/2, hcalXYSide/2,  thick / 2);

			G4RotationMatrix* rot = new G4RotationMatrix(0,0,0);
			rot->rotateY(3.14/2);
			G4UnionSolid* sUnion = new G4UnionSolid(baseName + "box",
					solid,
					solid_1,
					rot,
					trans);

			rot = new G4RotationMatrix(0,0,0);
			rot->rotateY(3.14+3.14/2);
			trans = G4ThreeVector(-hcalXYSide/2 - hcalThickness - thick/2.0 - 10,0,-hcalLen/2.0 -hcalThickness - thick/2.0+constShift);
			sUnion = new G4UnionSolid(baseName + "box",
					sUnion,
					solid_1,
					rot,
					trans);
			rot = new G4RotationMatrix(0,0,0);
			rot->rotateX(3.14/2);
			trans = G4ThreeVector(0,hcalXYSide/2 + hcalThickness + thick/2.0 + 10,-hcalLen/2.0 - hcalThickness - thick/2.0+constShift);
			solid_1 = new G4Box(baseName + "box",m_HCALSizeXY/2, hcalLen/2,  thick / 2);
			sUnion = new G4UnionSolid(baseName + "box",
					sUnion,
					solid_1,
					rot,
					trans);
			solid = sUnion;
			rot = new G4RotationMatrix(0,0,0);
			rot->rotateX(3.14/2+3.14);
			trans = G4ThreeVector(0,-hcalXYSide/2 - hcalThickness - thick/2.0 - 10 ,-hcalLen/2.0 - hcalThickness - thick/2.0+constShift);
			solid_1 = new G4Box(baseName + "box",m_HCALSizeXY/2, hcalLen/2,  thick / 2);
			sUnion = new G4UnionSolid(baseName + "box",
					sUnion,
					solid_1,
					rot,
					trans);
			   hcalThickness += thick/mm;
			solid = sUnion;

		}*/
	}
	return solid;
}
