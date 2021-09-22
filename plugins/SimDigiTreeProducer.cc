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
  edm::EDGetTokenT<std::vector<CaloParticle> > _puCaloPartToken;
  edm::EDGetTokenT<std::vector<CaloParticle> > _ootpuCaloPartToken;
  
  TTree *_outTree;
  
  UInt_t _run;
  UShort_t _lumi;
  UShort_t _bx;
  UShort_t _event;      
  
  //
  // it's not possible to get information about the BX of the OOT pu unfortunately
  //      https://github.com/cms-sw/cmssw/blob/master/SimDataFormats/CaloAnalysis/interface/CaloParticle.h
  // unless I modify the code by Badder: ReducedCaloParticleProducer
  //
  
  float _simenergy_EB[61200*3]; // 3 because I have   signal [0],  PU [1], OOT PU [2]
  float _digi_ped_subtracted_EB[61200*10];
  int   _ieta[61200];
  int   _iphi[61200];
  
  
  float _simenergy_EE[14648*3]; // 3 because I have   signal [0],  PU [1], OOT PU [2]
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
  _caloPartToken = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticleCollection"));
  _puCaloPartToken     = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("puCaloParticleCollection"));
  _ootpuCaloPartToken  = consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("ootpuCaloParticleCollection"));
  
  
  //   _pedsToken = myConsumesCollector.esConsumes<EcalPedestals, EcalPedestalsRcd>();
  
  _outTree = fs->make<TTree>("tree","tree");
  
  _outTree->Branch("run",               &_run,             "run/i");
  _outTree->Branch("lumi",              &_lumi,            "lumi/s");
  _outTree->Branch("bx",                &_bx,              "bx/s");
  _outTree->Branch("event",             &_event,           "event/i");
  
  _outTree->Branch("digi_ped_subtracted_EB",        _digi_ped_subtracted_EB,        "digi_ped_subtracted_EB[612000]/F"); // 61200*10
  _outTree->Branch("simenergy_EB",                  _simenergy_EB,        "simenergy_EB[183600]/F"); // 61200*3
  _outTree->Branch("ieta",                _ieta,                "ieta[61200]/I");
  _outTree->Branch("iphi",                _iphi,                "iphi[61200]/I");
  
  _outTree->Branch("digi_ped_subtracted_EE",        _digi_ped_subtracted_EE,        "digi_ped_subtracted_EE[146480]/F"); // 14648*10
  _outTree->Branch("simenergy_EE",                  _simenergy_EE,        "simenergy_EE[43944]/F"); // 14648*3
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
    for (int i=0; i<3; i++) _simenergy_EB[ixtal*3+i] = 0.;
    _ieta[ixtal] = -999;
    _iphi[ixtal] = -999;
  }
  for (int ixtal=0; ixtal < 14648; ixtal++) {
    for (int i=0; i<10; i++) _digi_ped_subtracted_EE[ixtal*10+i] = -999;
    for (int i=0; i<3; i++) _simenergy_EE[ixtal*3+i] = 0.;
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
  
  
//   
//   for reference:
//   https://github.com/cms-sw/cmssw/blob/master/SimDataFormats/CaloAnalysis/interface/SimCluster.h
//   https://github.com/cms-sw/cmssw/blob/master/SimDataFormats/CaloAnalysis/interface/CaloParticle.h#L72
// 
  
  // signal calo particles
  
  for (const auto& iCalo : *(caloParticles.product())) {
    
    // for each calo particle, get all energy deposits, and save them ...
    
    //     _simenergy_EB[ixtal] += energy
    //     _simenergy_EE[ixtal] += energy

    const auto& simClusters = iCalo.simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
      auto simCluster = simClusters[iSC];  
      auto hits_and_energies = simCluster->hits_and_energies();
      for(unsigned int i = 0; i < hits_and_energies.size(); i++){   
//         det id = DetId(hits_and_energies[i].first)
//         energy = hits_and_energies[i].second
        
        DetId id(hits_and_energies[i].first);
                
        if(id.subdetId()==EcalBarrel) {
          EBDetId eb_id(id);
          _simenergy_EB[eb_id.hashedIndex()*3] += hits_and_energies[i].second;
        }
        else if(id.subdetId()==EcalEndcap) {
          EEDetId ee_id(id);
          _simenergy_EE[ee_id.hashedIndex()*3] += hits_and_energies[i].second;
        }
        
      }  
    }
  }
  
  
  
  edm::Handle<std::vector<CaloParticle> > puCaloParticle;
  iEvent.getByToken(_puCaloPartToken,puCaloParticle);
  if (!puCaloParticle.isValid()) {
    std::cerr << "Analyze --> puCaloParticle not found" << std::endl; 
    return;
  }
  
  
  
  // in time pu calo particles
  
  for (const auto& iCalo : *(puCaloParticle.product())) {
    
    // for each calo particle, get all energy deposits, and save them ...
    
    const auto& simClusters = iCalo.simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
      auto simCluster = simClusters[iSC];  
      auto hits_and_energies = simCluster->hits_and_fractions();
      float totalenergy = simCluster->simEnergy();
//       std::cout << " totalenergy PU iSC [ " << iSC << "] = " << totalenergy << std::endl;
      for(unsigned int i = 0; i < hits_and_energies.size(); i++){   
        //         det id = DetId(hits_and_energies[i].first)
        //         energy = hits_and_energies[i].second

//         std::cout << "    [" << i << "] = " << (hits_and_energies[i].second) << std::endl;
        
        DetId id(hits_and_energies[i].first);
                
//         fraction = hAndE.second / totalenergy;
        
        if(id.subdetId()==EcalBarrel) {
          EBDetId eb_id(id);
          _simenergy_EB[eb_id.hashedIndex()*3 + 1] += ( (hits_and_energies[i].second) );
          totalenergy += ( (hits_and_energies[i].second) );
        }
        else if(id.subdetId()==EcalEndcap) {
          EEDetId ee_id(id);
          _simenergy_EE[ee_id.hashedIndex()*3 + 1] += ( (hits_and_energies[i].second) );
          totalenergy += ( (hits_and_energies[i].second) );
        }
        
      } 
//       std::cout << "   sum fraction = " << totalenergy << std::endl;
    }
  }
  
  
  
  
  edm::Handle<std::vector<CaloParticle> > ootpuCaloParticle;
  iEvent.getByToken(_ootpuCaloPartToken,ootpuCaloParticle);
  if (!ootpuCaloParticle.isValid()) {
    std::cerr << "Analyze --> ootpuCaloParticle not found" << std::endl; 
    return;
  }
  
  
  
  
  // out of time pu calo particles
  
  for (const auto& iCalo : *(ootpuCaloParticle.product())) {
      
    // for each calo particle, get all energy deposits, and save them ...
    
    const auto& simClusters = iCalo.simClusters();
    for(unsigned int iSC = 0; iSC < simClusters.size() ; iSC++){
      auto simCluster = simClusters[iSC];  
      auto hits_and_energies = simCluster->hits_and_fractions();
      float totalenergy = simCluster->simEnergy();
//       std::cout << " totalenergy OOTPU iSC [ " << iSC << "] = " << totalenergy << std::endl;
      
      for(unsigned int i = 0; i < hits_and_energies.size(); i++){   
        //         det id = DetId(hits_and_energies[i].first)
        //         energy = hits_and_energies[i].second
        
        DetId id(hits_and_energies[i].first);
                
        //         fraction = hAndE.second / totalenergy;
        
        if(id.subdetId()==EcalBarrel) {
          EBDetId eb_id(id);
          _simenergy_EB[eb_id.hashedIndex()*3 + 2] += ( (hits_and_energies[i].second) );
          totalenergy += ( (hits_and_energies[i].second) );
        }
        else if(id.subdetId()==EcalEndcap) {
          EEDetId ee_id(id);
          _simenergy_EE[ee_id.hashedIndex()*3 + 2] += ( (hits_and_energies[i].second) );
          totalenergy += ( (hits_and_energies[i].second) );
        }
      }  
//       std::cout << "   sum fraction = " << totalenergy << std::endl;
    }
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
