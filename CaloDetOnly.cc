#include <iostream>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

// Geant4
#include "G4RunManager.hh"
#include "G4GDMLParser.hh"
#include "FTFP_BERT.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4GDMLParser.hh"
#include "G4UIExecutive.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"

int main(int argc, char** argv)  {

    std::cout << "Hello CaloDetOnly!" << std::endl;

    G4RunManager* runManager = new G4RunManager;

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    runManager->SetUserInitialization(new FTFP_BERT);
    runManager->SetUserAction(new PrimaryGeneratorAction());
    runManager->SetUserInitialization(new DetectorConstruction());
    UImanager->ApplyCommand("/N03/det/setModel 2");
    runManager->Initialize();

    G4GDMLParser* parser = new G4GDMLParser();
    G4VPhysicalVolume* world =
            G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume();
    parser->Write("calo_model2.gdml", world);

    ui->SessionStart();
    delete ui;

    std::cout << "Bye CaloDetOnly!" << std::endl;
}
