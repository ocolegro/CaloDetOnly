#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "SamplingSection.hh"

#include "PurgMagTabulatedField3D.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>
#include <string>

class G4VSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class G4Colour;

/**
 @class DetectorConstruction
 @short builds a simple detector
 */
class DetectorConstruction: public G4VUserDetectorConstruction {
public:

	enum DetectorVersion {
		HE			 = 1,
		TEH			 = 2,

	};




	/**
	 @short CTOR
	 */
	DetectorConstruction(G4int ver = DetectorConstruction::TEH);

	void buildHGCALFHE(const unsigned aVersion);
	void buildHGCALBHE(const unsigned aVersion);
	/**
	 @short calorimeter structure (sampling sections)
	 */
	std::vector<SamplingSection> m_caloStruct;
	std::vector<SamplingSection> *getStructure() {
		return &m_caloStruct;
	}


	int getVersion() const {
		return version_;
	}

	const std::vector<G4LogicalVolume*> & getSiLogVol() {
		return m_logicSi;
	}
	const std::vector<G4LogicalVolume*> & getAbsLogVol() {
		return m_logicAbs;
	}

	/**
	 @short define the calorimeter materials
	 */
	void DefineMaterials();
	std::map<std::string, G4Material *> m_materials;
	std::map<std::string, G4double> m_dEdx;
	std::map<std::string, G4Colour> m_colours;

	/**
	 @short set magnetic field
	 */

	void SetMagField(G4double fieldValue);
	G4UniformMagField* m_magField;      //pointer to the magnetic field
    void SetMagField(char *fileName, G4double zOffset);
	PurgMagTabulatedField3D *p_magField;

	/**
	 @short DTOR
	 */
	~DetectorConstruction();

	/**
	 @short getters
	 */
	G4double GetECALSizeXY() {
		return m_ECALSizeXY;
	}
	G4double GetHCALSizeXY() {
		return m_HCALSizeXY;
	}
	G4double GetCalorSizeZ() {
		return m_CalorSizeZ;
	}
	G4double GetWorldSizeXY() {
		return m_WorldSizeXY;
	}
	G4double GetWorldSizeZ() {
		return m_WorldSizeZ;
	}

	unsigned initLayer() {
		return initLayer_;
	}
	inline void initLayer(int aVal) {
		 initLayer_ = aVal;
	}
	/**
	 @short build the detector
	 */

	G4VPhysicalVolume* Construct();

private:

	//detector version
	int version_;
	//integer to define detector model
	int	hcalThickness;
	/**
	 @short compute the calor dimensions
	 */
	void UpdateCalorSize();
	void buildTracker();
	void buildECal();
	void buildHCal();

	/**
	 @short build the calorimeter
	 */
	G4VPhysicalVolume* ConstructCalorimeter();

	void buildSectorStack(const unsigned sectorNum, const G4double & minL,
			const G4double & width);

	void fillInterSectorSpace(const unsigned sectorNum, const G4double & minL,
			const G4double & width);

	G4VSolid *constructSolid(std::string baseName, G4double thick,
			G4double zpos, const G4double & minL, const G4double & width,
			const size_t which_ele, const unsigned nEle);

	std::vector<G4Material*> m_SensitiveMaterial;
	G4int initLayer_;
	G4double m_ECALSizeXY, m_HCALSizeXY, m_ECALSizeZ,m_CalorSizeZ;
	G4double m_z0pos;
	G4double m_WorldSizeXY, m_WorldSizeZ;
	G4double m_nSectors, m_sectorWidth, m_interSectorWidth;
	G4double steelThick_;
	G4VSolid* m_solidWorld;    //pointer to the solid World 
	G4LogicalVolume* m_logicWorld;    //pointer to the logical World
	G4VPhysicalVolume* m_physWorld;     //pointer to the physical World  

	std::vector<G4LogicalVolume*> m_logicSi; //pointer to the logical Si volumes
	std::vector<G4LogicalVolume*> m_logicAbs; //pointer to the logical absorber volumes situated just before the si

	DetectorMessenger* m_detectorMessenger;  //pointer to the Messenger
};

#endif

