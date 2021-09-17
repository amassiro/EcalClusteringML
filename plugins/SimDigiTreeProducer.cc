// -*- C++ -*-
//
// Package:    EcalZee
// Class:      SimDigiTreeProducer
// 
/**\class SimDigiTreeProducer SimDigiTreeProducer.cc EcalZee/plugins/SimDigiTreeProducer.cc
 * 
 * Description: [one line class summary]
 * 
 * Implementation:
 * [Notes on implementation]
 */
//
// Original Author:  Andrea Massironi
//         Created:  Thu, 1 September 2021 10:09:05 GMT
//
//

//
// Sim energy part is based on Badder's dumper: 
//    https://github.com/bmarzocc/RecoSimStudies/
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/Framework/interface/EventSetup.h"



// Energy Sim dumper

#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"


// ECAL specific

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"




#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"


#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "TTree.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class SimDigiTreeProducer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit SimDigiTreeProducer(const edm::ParameterSet& ); // z, edm::ConsumesCollector& );
  ~SimDigiTreeProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
   
  edm::EDGetTokenT<EBDigiCollection> _token_ebdigi;
  edm::EDGetTokenT<EEDigiCollection> _token_eedigi;
  
  const EcalPedestals* _peds;
  
  edm::EDGetTokenT<std::vector<CaloParticle> > _caloPartToken;
  
  
  TTree *_outTree;
  
  UInt_t _run;
  UShort_t _lumi;
  UShort_t _bx;
  UShort_t _event;      
  
  
  float _simenergy_EB[61200*5]; // 5 because I have - 3 BX  [0], -2 BX  [1], - 1 BX   [2], nominal BX   [3], +1 BX    [4]
  float _digi_ped_subtracted_EB[61200*10];
  int   _ieta[61200];
  int   _iphi[61200];
  
  
  float _simenergy_EE[14648*5]; // 5 because I have - 3 BX  [0], -2 BX  [1], - 1 BX   [2], nominal BX   [3], +1 BX    [4]
  float _digi_ped_subtracted_EE[14648*10];
  int   _ix[14648];
  int   _iy[14648];
  int   _iz[14648];
    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SimDigiTreeProducer::SimDigiTreeProducer(const edm::ParameterSet& iConfig) //,  edm::ConsumesCollector& myConsumesCollector)

{
  //now do what ever initialization is needed
  usesResource("TFileService");
  
  
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  
  _token_ebdigi = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  _token_eedigi = consumes<EEDigiCollection>(iConfig.getParameter<edm::InputTag>("EEDigiCollection"));
  
  
  //   _pedsToken = myConsumesCollector.esConsumes<EcalPedestals, EcalPedestalsRcd>();
  
  _outTree = fs->make<TTree>("tree","tree");
  
  _outTree->Branch("run",               &_run,             "run/i");
  _outTree->Branch("lumi",              &_lumi,            "lumi/s");
  _outTree->Branch("bx",                &_bx,              "bx/s");
  _outTree->Branch("event",             &_event,           "event/i");
  
  _outTree->Branch("digi_ped_subtracted_EB",        _digi_ped_subtracted_EB,        "digi_ped_subtracted_EB[612000]/F"); // 61200*10
  _outTree->Branch("simenergy_EB",                  _simenergy_EB,        "simenergy_EB[306000]/F"); // 61200*5
  _outTree->Branch("ieta",                _ieta,                "ieta[61200]/I");
  _outTree->Branch("iphi",                _iphi,                "iphi[61200]/I");
  
  _outTree->Branch("digi_ped_subtracted_EE",        _digi_ped_subtracted_EE,        "digi_ped_subtracted_EE[146480]/F"); // 14648*10
  _outTree->Branch("simenergy_EE",                  _simenergy_EE,        "simenergy_EE[73240]/F"); // 14648*5
  _outTree->Branch("ix",                  _ix,                  "ix[14648]/I");
  _outTree->Branch("iy",                  _iy,                  "iy[14648]/I");
  _outTree->Branch("iz",                  _iz,                  "iz[14648]/I");
  
  
  
}


SimDigiTreeProducer::~SimDigiTreeProducer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
SimDigiTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  _run = iEvent.eventAuxiliary().run();
  _lumi = iEvent.eventAuxiliary().luminosityBlock();
  _bx = iEvent.eventAuxiliary().bunchCrossing();
  _event = iEvent.eventAuxiliary().event();
  
  
  
  //---- pedestals
  edm::ESHandle< EcalPedestals > ecalPedestals;
  iSetup.get<EcalPedestalsRcd>().get(ecalPedestals);
  _peds = ecalPedestals.product();
  
  
  //---- digis
  edm::Handle<EBDigiCollection> ebdigihandle;
  const EBDigiCollection *ebdigis = NULL;
  edm::Handle<EEDigiCollection> eedigihandle;
  const EEDigiCollection *eedigis = NULL;
  
  iEvent.getByToken(_token_ebdigi,ebdigihandle);
  ebdigis = ebdigihandle.product();
  iEvent.getByToken(_token_eedigi,eedigihandle);
  eedigis = eedigihandle.product();
  
  
  
  
  //---- setup default
  for (int ixtal=0; ixtal < 61200; ixtal++) {
    for (int i=0; i<10; i++) _digi_ped_subtracted_EB[ixtal*10+i] = -999;
    for (int i=0; i<5; i++) _simenergy_EB[ixtal*5+i] = 0.;
    _ieta[ixtal] = -999;
    _iphi[ixtal] = -999;
  }
  for (int ixtal=0; ixtal < 14648; ixtal++) {
    for (int i=0; i<10; i++) _digi_ped_subtracted_EE[ixtal*10+i] = -999;
    for (int i=0; i<5; i++) _simenergy_EB[ixtal*5+i] = 0.;
    _ix[ixtal] = -999;
    _iy[ixtal] = -999;
    _iz[ixtal] = -999;
  }
  
  
  
  //---- dump digis
  
  for (EBDigiCollection::const_iterator itdigi = ebdigis->begin(); itdigi != ebdigis->end(); itdigi++ ) {
    
    float pedestal = 0;    
    DetId id = (EBDetId&)((*itdigi));
    pedestal = float((_peds->find(id))->mean_x12);
    
    //                                                           0xFFF = 4095
    for (int iSample = 0; iSample < 10; iSample++) {
      float value = ( int( (*itdigi) [iSample] ) & 0xFFF );
      _digi_ped_subtracted_EB[((EBDetId&)((*itdigi))).hashedIndex() *10 + iSample] = value - pedestal;
    }
    _ieta[((EBDetId&)((*itdigi))).hashedIndex()] = ((EBDetId&)((*itdigi))).ieta();
    _iphi[((EBDetId&)((*itdigi))).hashedIndex()] = ((EBDetId&)((*itdigi))).iphi();      
  }
  
  
  for (EEDigiCollection::const_iterator itdigi = eedigis->begin(); itdigi != eedigis->end(); itdigi++ ) {
    
    float pedestal = 0;
    DetId id = (EEDetId&)((*itdigi));
    pedestal = float((_peds->find(id))->mean_x12);
    
    //                                                           0xFFF = 4095
    for (int iSample = 0; iSample < 10; iSample++) {
      float value = ( int( (*itdigi) [iSample] ) & 0xFFF );
      _digi_ped_subtracted_EE[((EEDetId&)((*itdigi))).hashedIndex() *10 + iSample] = value - pedestal;
    }
    _ix[((EEDetId&)((*itdigi))).hashedIndex()] = ((EEDetId&)((*itdigi))).ix();
    _iy[((EEDetId&)((*itdigi))).hashedIndex()] = ((EEDetId&)((*itdigi))).iy();      
    _iz[((EEDetId&)((*itdigi))).hashedIndex()] = ((EEDetId&)((*itdigi))).zside();      
  }
  
  
  
  
  //---- dump sim energy
  
  edm::Handle<std::vector<CaloParticle> > caloParticles;
  iEvent.getByToken(_caloPartToken,caloParticles);
  if (!caloParticles.isValid()) {
    std::cerr << "Analyze --> caloParticles not found" << std::endl; 
    return;
  }
  
  
  std::vector<CaloParticle> caloParts;
  std::vector<GlobalPoint> caloParts_position;
  int caloParticle_size = 0;
  for(const auto& iCalo : *(caloParticles.product())) {
    caloParticle_size++;
//     std::vector<std::pair<DetId, float> > caloParticle_hitsAndEnergies = *getHitsAndEnergiesCaloPart(&iCalo,-1.);
//     GlobalPoint caloParticle_position = calculateAndSetPositionActual(&caloParticle_hitsAndEnergies, 7.4, 3.1, 1.2, 4.2, 0.89, 0.,true);
//     if(caloParticle_position == GlobalPoint(-999999., -999999., -999999.)){
//       std::cout << "Invalid position for caloParticle, skipping caloParticle!" << std::endl;
//       continue;
//     }    
//     
//     hitsAndEnergies_CaloPart.push_back(caloParticle_hitsAndEnergies);
//     caloParts_position.push_back(caloParticle_position);
//     caloParts.push_back(iCalo); 
  }
  
  
  for(unsigned int iCalo=0; iCalo<caloParts.size(); iCalo++){
  
   
  }
  
  
  _outTree->Fill();  
  
}




// ------------ method called once each job just before starting event loop  ------------
void 
SimDigiTreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SimDigiTreeProducer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SimDigiTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimDigiTreeProducer);
