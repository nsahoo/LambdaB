// -*- C++ -*-
//
// Package:    LambdaB
// Class:      LambdaB
// 
/**\class LambdaB LambdaB.cc BpHaNA/LambdaB/src/LambdaB.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//********************************************************************
//     Ntuplizer code for /\b -> /\(->p+pi) mu+ mu-                  *
//                                                                   *
// Original Author:  Niladribihari Sahoo,42 3-024,+41227662373,      *
//         Created:  Tue Dec  8 23:43:26 CET 2015                    *
//********************************************************************
// $Id$
//
//

//-----------------------
// system include files
//-----------------------
#include <memory>

//----------------------
// user include files
//----------------------
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TH1.h>


using namespace std;


const int MUONMINUS_PDG_ID = 13;
const int PIONPLUS_PDG_ID = 211;
const int PROTON_PDG_ID = 2212;
const int LAMBDA_PDG_ID = 3122;
const int LAMBDAB_PDG_ID = 5122;
const int JPSI_PDG_ID = 443;
const int PSI2S_PDG_ID = 100443;

const double PI = 3.141592653589793;

//-------------
// Structures  
//-------------                                                                                                                                              
struct HistArgs{
  char name[128];
  char title[128];
  int n_bins;
  double x_min;
  double x_max;
};

enum HistName{
  h_events,
  h_mupt,
  h_mueta,
  h_mumdcabs,
  h_mumutrkr,
  h_mumutrkz,

  h_mumudca,
  h_mumuvtxcl,
  h_mumupt,
  h_mumumass,
  h_mumulxybs,

  h_mumucosalphabs,
  h_trkpt,
  h_trkdcasigbs,
  h_lbvtxchisq,
  h_lbvtxcl,

  h_lzmass,   // lb = LambdaB, lz = Lambda0
  h_lbmass,

  kHistNameSize
};

//--------------------
// Global hist args                                                                                                                                    
//--------------------
HistArgs hist_args[kHistNameSize] = {
  // name, title, n_bins, x_min, x_max                                                                                                                       

  {"h_events", "Processed Events", 1, 0, 1},
  {"h_mupt", "Muon pT; [GeV]", 100, 0, 30},
  {"h_mueta", "Muon eta", 100, 0, 3},
  {"h_mumdcabs", "#mu^{-} DCA beam spot; DCA [cm]", 100, 0, 10},
  {"h_mumutrkr", "#mu^{+}#mu^{-} distance in phi-eta; [cm]", 100, 0, 50},
  {"h_mumutrkz", "#mu^{+}#mu^{-} distance in Z; [cm]", 100, 0, 100},

  {"h_mumudca",  "#mu^{+}#mu^{-} DCA; [cm]", 100, 0, 20},
  {"h_mumuvtxcl",  "#mu^{+}#mu^{-} vertex CL", 100, 0, 1},
  {"h_mumupt",    "#mu^{+}#mu^{-} pT ; pT [GeV]", 100, 0, 50},
  {"h_mumumass", "#mu^{+}#mu^{-} invariant mass; M(#mu^{+}#mu^{-}) [GeV/c^{2}]",
   100, 2, 20},
  {"h_mumulxybs", "#mu^{+}#mu^{-} Lxy #sigma beam spot", 100, 0, 100},

  {"h_mumucosalphabs", "#mu^{+}#mu^{-} cos #alpha beam spot", 100, 0, 1},
  {"h_trkpt", "Pion track pT; pT [GeV]", 100, 0, 20},
  {"h_trkdcasigbs", "Pion track DCA/#sigma beam spot; DCA/#sigma", 1000, 0, 100},
  {"h_lbvtxchisq", "#lambda_{b} decay vertex chisq", 100, 0, 1000},
  {"h_lbvtxcl", "#lambda_{b} decay vertex CL", 100, 0, 1},

  {"h_lzmass", "#lambda^{0} mass; M(#lambda^{0}) [GeV/^{2}]", 100, 0, 20},   // Lambda0 mass

  {"h_lbmass", "#lambda_{b} mass; M(#lambda_{b}) [GeV]", 100, 0, 20},     // LambdaB mass

};

//--------------------
// Define histograms  
//--------------------                                                                                                                                           
TH1F *histos[kHistNameSize];


//---------------------
// class declaration
//---------------------

class LambdaB : public edm::EDAnalyzer {
   public:
      explicit LambdaB(const edm::ParameterSet&);
      ~LambdaB();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);


  bool buildLbToLzMuMu(const edm::Event &);

  void calLS (double, double, double, double, double, double, double,
              double, double,  double, double, double, double, double,
              double, double, double, double, double*, double*);

  void calCosAlpha (double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double, double,
                    double, double, double, double,
                    double*, double*);

  void calCosAlpha2d (double, double, double, double, double,   
                      double, double, double, double, double,
                      double, double, double, double,
                      double, double, double, double,
                      double*, double*);

  void calCtau(RefCountedKinematicTree, double &, double &);
  double calEta(double, double, double);
  double calPhi(double, double, double);
  double calEtaPhiDistance (double, double, double, double, double, double);
  void clearVariables();

  bool hasBeamSpot(const edm::Event&);

  bool calClosestApproachTracks(const reco::TransientTrack,
                                const reco::TransientTrack,
                                double&, double &, double &);

  bool hasGoodLzVertex(const reco::TransientTrack, const reco::TransientTrack, double &);  // Lambda0 vertex

  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaPoint (const reco::TransientTrack, const GlobalPoint,
                             double, double &, double &);

  bool hasGoodLbMass(RefCountedKinematicTree, double &);  // LambdaB mass

  bool hasGoodLbVertex(const reco::TransientTrack, const reco::TransientTrack,   // LambdaB vertex
                       const reco::TransientTrack, const reco::TransientTrack,
                       double &, double &, double &, RefCountedKinematicTree &);

  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
                          reco::TransientTrack &, reco::TransientTrack &,
                          double &, double &, double &, double &, double &,
                          double &, double &, double &);

  bool hasGoodTrack(const edm::Event&, const pat::GenericParticle, double &);

  bool hasPrimaryVertex(const edm::Event &);

  void hltReport(const edm::Event&);

  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);
  bool matchMuonTracks (const edm::Event&, const vector<reco::TrackRef>);
  bool matchPrimaryVertexTracks ();

  void saveLbToLzMuMu(const RefCountedKinematicTree);
  void saveLbVertex(RefCountedKinematicTree);
  void saveLbCosAlpha(RefCountedKinematicTree);
  void saveLbCosAlpha2d(RefCountedKinematicTree);    
  void saveLbLsig(RefCountedKinematicTree);
  void saveLbCtau(RefCountedKinematicTree);

  void saveGenInfo(const edm::Event&);
  void saveSoftMuonVariables(pat::Muon, pat::Muon, reco::TrackRef, reco::TrackRef);
  void saveDimuVariables(double, double, double, double, double, double,
                         double, double, double, double, double, double,
                         double, double);
  void saveMuonTriggerMatches(const pat::Muon, const pat::Muon);
  void saveTruthMatch(const edm::Event& iEvent);


      // ----------member data ---------------------------


  // --- begin input from python file ---                                                                                                                     
  string OutputFileName_;
  bool BuildLbToLzMuMu_;

  //----------------------
  // particle properties  
  //----------------------                                                                                                                                
  ParticleMass MuonMass_;
  float MuonMassErr_;
  ParticleMass PionMass_;
  float PionMassErr_;
  ParticleMass ProtonMass_;
  float ProtonMassErr_;
  double LbMass_;

  //----------
  // labels   
  //----------                                                                                                                                            
  edm::InputTag GenParticlesLabel_;
  edm::InputTag TriggerResultsLabel_;
  edm::InputTag BeamSpotLabel_;
  edm::InputTag VertexLabel_;
  edm::InputTag MuonLabel_;
  //edm::InputTag KshortLabel_;
  edm::InputTag TrackLabel_;
  vector<string> TriggerNames_;
  vector<string> LastFilterNames_;

  //---------------
  // gen particle  
  //---------------                                                                                                                                       
  bool   IsMonteCarlo_;
  bool   KeepGENOnly_;
  double TruthMatchMuonMaxR_;
  double TruthMatchPionMaxR_;
  double TruthMatchProtonMaxR_;

  //---------------------
  // pre-selection cuts  
  //---------------------                                                                                                                                       
  double MuonMinPt_;
  double MuonMaxEta_;
  double MuonMaxDcaBs_;
  double TrkMinPt_;
  double TrkMinDcaSigBs_;
  double TrkMaxR_;
  double TrkMaxZ_;
  double MuMuMaxDca_;
  double MuMuMinVtxCl_;
  double MuMuMinPt_;
  double MuMuMinInvMass_;
  double MuMuMaxInvMass_;
  double MuMuMinLxySigmaBs_;
  double MuMuMinCosAlphaBs_;
  double LzMinMass_;
  double LzMaxMass_;
  double LbMinVtxCl_;
  double LbMinMass_;
  double LbMaxMass_;

  //--------------------
  // Across the event   
  //--------------------                                                                                                                                      
  map<string, string> mapTriggerToLastFilter_;
  reco::BeamSpot beamSpot_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  //-----------------
  // Root Variables                                                                                                                                        
  //-----------------
  TFile* fout_;
  TTree* tree_;

  unsigned int run, event, lumiblock, nprivtx;
  vector<string> *triggernames;
  vector<int> *triggerprescales;

  //----------
  // dimuon   
  //----------                                                                                                                                                     
  vector<double> *mumdcabs, *mumdcabserr, *mumpx, *mumpy, *mumpz;
  vector<double> *mupdcabs, *mupdcabserr, *muppx, *muppy, *muppz;
  vector<double> *mumutrkr, *mumutrkz , *mumudca;
  vector<double> *mumuvtxcl, *mumulsbs, *mumulsbserr;
  vector<double> *mumucosalphabs, *mumucosalphabserr;
  vector<double> *mumumass, *mumumasserr;

  //-----------------------
  // soft muon variables   
  //-----------------------                                                                                                                                    
  vector<bool>   *mumisgoodmuon, *mupisgoodmuon ;
  vector<int>    *mumnpixhits, *mupnpixhits, *mumnpixlayers, *mupnpixlayers;
  vector<int>    *mumntrkhits, *mupntrkhits, *mumntrklayers, *mupntrklayers;
  vector<double> *mumnormchi2, *mupnormchi2;
  vector<double> *mumdxyvtx, *mupdxyvtx, *mumdzvtx, *mupdzvtx;
  vector<string> *mumtriglastfilter, *muptriglastfilter;
  vector<double> *mumpt, *muppt, *mumeta, *mupeta;

  //---------------
  // pion track                                                                                                                                                       
  //---------------
  vector<int> *trkchg; // +1 for pi+, -1 for pi-                                                                                                                
  vector<double> *trkpx, *trkpy, *trkpz, *trkpt;
  vector<double> *trkdcabs, *trkdcabserr;

  vector<double> *lzpx, *lzpy, *lzpz;
  vector<double> *lzvtxx, *lzvtxy, *lzvtxz;

  //--------------------
  // proton, pion track 
  //--------------------                                                                                                                                       
  vector<double> *prpx, *prpy, *prpz;
  vector<double> *pipx, *pipy, *pipz;

  //-----------
  // Lambda0, anti-Lambda0 
  //-----------
  vector<double> *lzmass, *lzmasserr, *lzbarmass, *lzbarmasserr;

  //-----------
  // LambdaB
  //-----------
  int nb;
  vector<double> *lbpx, *lbpxerr, *lbpy, *lbpyerr, *lbpz, *lbpzerr, *lbmass, *lbmasserr;
  vector<double> *lbvtxcl, *lbvtxx, *lbvtxxerr, *lbvtxy, *lbvtxyerr, *lbvtxz, *lbvtxzerr;
  vector<double> *lbcosalphabs, *lbcosalphabserr, *lbcosalphabs2d, *lbcosalphabs2derr, *lblsbs, *lblsbserr, *lbctau, *lbctauerr; 

  //---------------
  // anti-LambdaB
  //---------------                                                                                                                                                
  vector<double> *lbbarmass, *lbbarmasserr;

  //---------
  // For MC  
  //---------                                                                                                                                               
  double genlbpx, genlbpy, genlbpz;
  double genlzpx, genlzpy, genlzpz;
  double genlzvtxx, genlzvtxy, genlzvtxz;

  int genprchg;
  double genprpx, genprpy, genprpz;
  int genpichg;
  double genpipx, genpipy, genpipz;

  double genmumpx, genmumpy, genmumpz;
  double genmuppx, genmuppy, genmuppz;
  double genpippx, genpippy, genpippz;
  double genpimpx, genpimpy, genpimpz;

  string decname;

  vector<bool> *istruemum, *istruemup, *istruepr, *istruepi, *istruelb;

  //-----------------------
  // variables to monitor  
  //-----------------------                                                                                                                                    
  TDatime t_begin_ , t_now_ ;
  int n_processed_, n_selected_;


};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//-------------------------------
// constructors and destructor
//-------------------------------
LambdaB::LambdaB(const edm::ParameterSet& iConfig):
  OutputFileName_(iConfig.getParameter<string>("OutputFileName")),
  BuildLbToLzMuMu_(iConfig.getUntrackedParameter<bool>("BuildLbToLzMuMu")),

  //----------------------
  // particle properties  
  //----------------------                                                                                                                                     
  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),
  PionMass_(iConfig.getUntrackedParameter<double>("PionMass")),
  PionMassErr_(iConfig.getUntrackedParameter<double>("PionMassErr")),
  ProtonMass_(iConfig.getUntrackedParameter<double>("ProtonMass")),
  ProtonMassErr_(iConfig.getUntrackedParameter<double>("ProtonMassErr")),
  LbMass_(iConfig.getUntrackedParameter<double>("LbMass")),

  //---------
  // labels  
  //---------                                                                                                                                                  
  GenParticlesLabel_(iConfig.getParameter<edm::InputTag>("GenParticlesLabel")),
  TriggerResultsLabel_(iConfig.getParameter<edm::InputTag>("TriggerResultsLabel")),
  BeamSpotLabel_(iConfig.getParameter<edm::InputTag>("BeamSpotLabel")),
  VertexLabel_(iConfig.getParameter<edm::InputTag>("VertexLabel")),
  MuonLabel_(iConfig.getParameter<edm::InputTag>("MuonLabel")),
  TrackLabel_(iConfig.getParameter<edm::InputTag>("TrackLabel")),
  TriggerNames_(iConfig.getParameter< vector<string> >("TriggerNames")),
  LastFilterNames_(iConfig.getParameter< vector<string> >("LastFilterNames")),

  //----------------
  // gen particle   
  //----------------                                                                                                                                          
  IsMonteCarlo_(iConfig.getUntrackedParameter<bool>("IsMonteCarlo")),
  KeepGENOnly_(iConfig.getUntrackedParameter<bool>("KeepGENOnly")),
  TruthMatchMuonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchMuonMaxR")),
  TruthMatchPionMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchPionMaxR")),
  TruthMatchProtonMaxR_(iConfig.getUntrackedParameter<double>("TruthMatchProtonMaxR")),
  
  //--------------------
  // pre-selection cuts 
  //--------------------                                                                                                                                      
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),

  TrkMinPt_(iConfig.getUntrackedParameter<double>("TrkMinPt")),
  TrkMinDcaSigBs_(iConfig.getUntrackedParameter<double>("TrkMinDcaSigBs")),
  TrkMaxR_(iConfig.getUntrackedParameter<double>("TrkMaxR")),
  TrkMaxZ_(iConfig.getUntrackedParameter<double>("TrkMaxZ")),

  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinPt_(iConfig.getUntrackedParameter<double>("MuMuMinPt")),
  MuMuMinInvMass_(iConfig.getUntrackedParameter<double>("MuMuMinInvMass")),
  MuMuMaxInvMass_(iConfig.getUntrackedParameter<double>("MuMuMaxInvMass")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")),
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")),

  LzMinMass_(iConfig.getUntrackedParameter<double>("LzMinMass")),
  LzMaxMass_(iConfig.getUntrackedParameter<double>("LzMaxMass")),
  LbMinVtxCl_(iConfig.getUntrackedParameter<double>("LbMinVtxCl")),
  LbMinMass_(iConfig.getUntrackedParameter<double>("LbMinMass")),
  LbMaxMass_(iConfig.getUntrackedParameter<double>("LbMaxMass")),

  tree_(0),
  triggernames(0), triggerprescales(0),
  mumdcabs(0), mumdcabserr(0), mumpx(0), mumpy(0), mumpz(0),
  mupdcabs(0),  mupdcabserr(0), muppx(0),  muppy(0), muppz(0),
  mumutrkr(0), mumutrkz(0), mumudca(0),  mumuvtxcl(0),  mumulsbs(0),
  mumulsbserr(0), mumucosalphabs(0),  mumucosalphabserr(0),
  mumumass(0), mumumasserr(0),
  mumisgoodmuon(0), mupisgoodmuon(0),
  mumnpixhits(0), mupnpixhits(0), mumnpixlayers(0), mupnpixlayers(0),
  mumntrkhits(0), mupntrkhits(0), mumntrklayers(0), mupntrklayers(0),
  mumnormchi2(0), mupnormchi2(0), mumdxyvtx(0), mupdxyvtx(0),
  mumdzvtx(0), mupdzvtx(0), mumtriglastfilter(0), muptriglastfilter(0),
  mumpt(0), muppt(0), mumeta(0), mupeta(0),

  trkchg(0), trkpx(0), trkpy(0), trkpz(0), trkpt(0),
  trkdcabs(0), trkdcabserr(0),

  lzpx(0), lzpy(0), lzpz(0),
  lzvtxx(0), lzvtxy(0), lzvtxz(0),

  prpx(0), prpy(0), prpz(0),
  pipx(0), pipy(0), pipz(0),

  lzmass(0), lzmasserr(0), lzbarmass(0), lzbarmasserr(0),

  nb(0), lbpx(0), lbpxerr(0), lbpy(0), lbpyerr(0), lbpz(0), lbpzerr(0), lbmass(0), lbmasserr(0),
  lbvtxcl(0), lbvtxx(0), lbvtxxerr(0), lbvtxy(0), lbvtxyerr(0), lbvtxz(0), lbvtxzerr(0),
  lbcosalphabs(0), lbcosalphabserr(0), lbcosalphabs2d(0), lbcosalphabs2derr(0), lblsbs(0), lblsbserr(0), lbctau(0), lbctauerr(0),

  lbbarmass(0), lbbarmasserr(0),

  genlbpx(0), genlbpy(0), genlbpz(0),
  genlzpx(0), genlzpy(0), genlzpz(0), genlzvtxx(0), genlzvtxy(0), genlzvtxz(0),

  genprchg(0),
  genprpx(0), genprpy(0), genprpz(0),

  genpichg(0),
  genpipx(0), genpipy(0), genpipz(0),

  genmumpx(0), genmumpy(0), genmumpz(0),
  genmuppx(0), genmuppy(0), genmuppz(0),

  genpippx(0), genpippy(0), genpippz(0),
  genpimpx(0), genpimpy(0), genpimpz(0),

  decname(""),
  istruemum(0), istruemup(0), istruepr(0), istruepi(0), istruelb(0)


{
   //now do what ever initialization is needed
  assert(TriggerNames_.size() == LastFilterNames_.size());
  for (size_t i = 0; i < TriggerNames_.size(); ++i)
    mapTriggerToLastFilter_[TriggerNames_[i]] = LastFilterNames_[i];

}


LambdaB::~LambdaB()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
LambdaB::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  clearVariables();
  
  run = iEvent.id().run() ;
  event = iEvent.id().event() ;
  lumiblock = iEvent.luminosityBlock();

  n_processed_ += 1;
  histos[h_events]->Fill(0);

  if (IsMonteCarlo_) saveGenInfo(iEvent);

  hltReport(iEvent);

  if ( KeepGENOnly_){
    tree_->Fill();
    n_selected_ += 1;
  }else{
    if ( hasBeamSpot(iEvent) ) {
      iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle_);
      if ( bFieldHandle_.isValid() && hasPrimaryVertex(iEvent) ) {
	buildLbToLzMuMu(iEvent) ;
	if (IsMonteCarlo_) saveTruthMatch(iEvent);
	n_selected_ += 1;               
      }
    }
   
    if (IsMonteCarlo_ || nb > 0){ // Keep failed events for MC to calculate reconstruction efficiency.
      tree_->Fill();
    }
  }


  clearVariables();

}


// ------------ method called once each job just before starting event loop  ------------
void 
LambdaB::beginJob()
{

  t_begin_.Set();
  printf("\n ---------- Begin Job ---------- \n");
  t_begin_.Print();

  n_processed_ = 0;
  n_selected_ = 0;


  fout_ = new TFile(OutputFileName_.c_str(), "RECREATE");
  fout_->cd();

  for(int i=0; i<kHistNameSize; i++) {
    histos[i] = new TH1F(hist_args[i].name, hist_args[i].title,
			 hist_args[i].n_bins,
			 hist_args[i].x_min, hist_args[i].x_max);

  }


  tree_ = new TTree ("tree", "LambdaB");

  tree_->Branch("run", &run, "run/i");
  tree_->Branch("event", &event, "event/i");
  tree_->Branch("lumiblock", &lumiblock, "lumiblock/i");
  tree_->Branch("nprivtx", &nprivtx, "nprivtx/i");
  tree_->Branch("triggernames", &triggernames);
  tree_->Branch("triggerprescales", &triggerprescales);
  tree_->Branch("mumdcabs", &mumdcabs);
  tree_->Branch("mumdcabserr", &mumdcabserr);
  tree_->Branch("mumpx", &mumpx);
  tree_->Branch("mumpy", &mumpy);
  tree_->Branch("mumpz", &mumpz);
  tree_->Branch("mupdcabs", &mupdcabs);
  tree_->Branch("mupdcabserr", &mupdcabserr);
  tree_->Branch("muppx", &muppx);
  tree_->Branch("muppy", &muppy);
  tree_->Branch("muppz", &muppz);
  tree_->Branch("mumutrkr", &mumutrkr);
  tree_->Branch("mumutrkz", &mumutrkz);
  tree_->Branch("mumudca", &mumudca);
  tree_->Branch("mumuvtxcl", &mumuvtxcl);
  tree_->Branch("mumulsbs", &mumulsbs);
  tree_->Branch("mumulsbserr", &mumulsbserr);
  tree_->Branch("mumucosalphabs", &mumucosalphabs);
  tree_->Branch("mumucosalphabserr", &mumucosalphabserr);
  tree_->Branch("mumumass", &mumumass);
  tree_->Branch("mumumasserr", &mumumasserr);
  tree_->Branch("mumisgoodmuon", &mumisgoodmuon);
  tree_->Branch("mupisgoodmuon", &mupisgoodmuon);
  tree_->Branch("mumnpixhits", &mumnpixhits);
  tree_->Branch("mupnpixhits", &mupnpixhits);
  tree_->Branch("mumnpixlayers", &mumnpixlayers);
  tree_->Branch("mupnpixlayers", &mupnpixlayers);
  tree_->Branch("mumntrkhits", &mumntrkhits);
  tree_->Branch("mupntrkhits", &mupntrkhits);
  tree_->Branch("mumntrklayers", &mumntrklayers);
  tree_->Branch("mupntrklayers", &mupntrklayers);
  tree_->Branch("mumnormchi2", &mumnormchi2);
  tree_->Branch("mupnormchi2", &mupnormchi2);
  tree_->Branch("mumdxyvtx", &mumdxyvtx);
  tree_->Branch("mupdxyvtx", &mupdxyvtx);
  tree_->Branch("mumdzvtx", &mumdzvtx);
  tree_->Branch("mupdzvtx", &mupdzvtx);
  tree_->Branch("mumtriglastfilter", &mumtriglastfilter);
  tree_->Branch("muptriglastfilter", &muptriglastfilter);
  tree_->Branch("mumpt", &mumpt);
  tree_->Branch("muppt", &muppt);
  tree_->Branch("mumeta", &mumeta);
  tree_->Branch("mupeta", &mupeta);
  tree_->Branch("trkchg", &trkchg);
  tree_->Branch("trkpx", &trkpx);
  tree_->Branch("trkpy", &trkpy);
  tree_->Branch("trkpz", &trkpz);
  tree_->Branch("trkpt", &trkpt);
  tree_->Branch("trkdcabs", &trkdcabs);
  tree_->Branch("trkdcabserr", &trkdcabserr);
  tree_->Branch("lzpx", &lzpx);
  tree_->Branch("lzpy", &lzpy);
  tree_->Branch("lzpz", &lzpz);
  tree_->Branch("lzvtxx", &lzvtxx);
  tree_->Branch("lzvtxy", &lzvtxy);
  tree_->Branch("lzvtxz", &lzvtxz);
  tree_->Branch("prpx", &prpx);
  tree_->Branch("prpy", &prpy);
  tree_->Branch("prpz", &prpz);
  tree_->Branch("pipx", &pipx);
  tree_->Branch("pipy", &pipy);
  tree_->Branch("pipz", &pipz);
  tree_->Branch("lzmass", &lzmass);
  tree_->Branch("lzmasserr", &lzmasserr);
  tree_->Branch("lzbarmass", &lzbarmass);
  tree_->Branch("lzbarmasserr", &lzbarmasserr);
  tree_->Branch("nb", &nb, "nb/I");
  tree_->Branch("lbpx", &lbpx);
  tree_->Branch("lbpxerr", &lbpxerr);
  tree_->Branch("lbpy", &lbpy);
  tree_->Branch("lbpyerr", &lbpyerr);
  tree_->Branch("lbpz", &lbpz);
  tree_->Branch("lbpzerr", &lbpzerr);
  tree_->Branch("lbmass", &lbmass);
  tree_->Branch("lbmasserr", &lbmasserr);
  tree_->Branch("lbvtxcl", &lbvtxcl);
  tree_->Branch("lbvtxx", &lbvtxx);
  tree_->Branch("lbvtxxerr", &lbvtxxerr);
  tree_->Branch("lbvtxy", &lbvtxy);
  tree_->Branch("lbvtxyerr", &lbvtxyerr);
  tree_->Branch("lbvtxz", &lbvtxz);
  tree_->Branch("lbvtxzerr", &lbvtxzerr);
  tree_->Branch("lbcosalphabs", &lbcosalphabs);
  tree_->Branch("lbcosalphabserr", &lbcosalphabserr);
  tree_->Branch("lbcosalphabs2d", &lbcosalphabs2d);           
  tree_->Branch("lbcosalphabs2derr", &lbcosalphabs2derr);    
  tree_->Branch("lblsbs", &lblsbs);
  tree_->Branch("lblsbserr", &lblsbserr);
  tree_->Branch("lbctau", &lbctau);
  tree_->Branch("lbctauerr", &lbctauerr);
  tree_->Branch("lbbarmass", &lbbarmass);
  tree_->Branch("lbbarmasserr", &lbbarmasserr);

  if (IsMonteCarlo_) {

    tree_->Branch("genlbpx",      &genlbpx     , "genlbpx/D"    );
    tree_->Branch("genlbpy",      &genlbpy     , "genlbpy/D"    );
    tree_->Branch("genlbpz",      &genlbpz     , "genlbpz/D"    );
    tree_->Branch("genlzpx",      &genlzpx   , "genlzpx/D"  );
    tree_->Branch("genlzpy",      &genlzpy   , "genlzpy/D"  );
    tree_->Branch("genlzpz",      &genlzpz   , "genlzpz/D"  );
    tree_->Branch("genlzvtxx",     &genlzvtxx    , "genlzvtxx/D"   );
    tree_->Branch("genlzvtxy",     &genlzvtxy    , "genlzvtxy/D"   );
    tree_->Branch("genlzvtxz",     &genlzvtxz    , "genlzvtxz/D"   );
    tree_->Branch("genprchg",    &genprchg   , "genprchg/I"   );
    tree_->Branch("genprpx",     &genprpx    , "genprpx/D"   );
    tree_->Branch("genprpy",     &genprpy    , "genprpy/D"   );
    tree_->Branch("genprpz",     &genprpz    , "genprpz/D"   );
    tree_->Branch("genpichg",    &genpichg   , "genpichg/I"   );
    tree_->Branch("genpipx",     &genpipx    , "genpipx/D"   );
    tree_->Branch("genpipy",     &genpipy    , "genpipy/D"   );
    tree_->Branch("genpipz",     &genpipz    , "genpipz/D"   );
    tree_->Branch("genmumpx",    &genmumpx   , "genmumpx/D"  );
    tree_->Branch("genmumpy",    &genmumpy   , "genmumpy/D"  );
    tree_->Branch("genmumpz",    &genmumpz   , "genmumpz/D"  );
    tree_->Branch("genmuppx",    &genmuppx   , "genmuppx/D"  );
    tree_->Branch("genmuppy",    &genmuppy   , "genmuppy/D"  );
    tree_->Branch("genmuppz",    &genmuppz   , "genmuppz/D"  );
    tree_->Branch("genpimpx",    &genpimpx   , "genpimpx/D"  );
    tree_->Branch("genpimpy",    &genpimpy   , "genpimpy/D"  );
    tree_->Branch("genpimpz",    &genpimpz   , "genpimpz/D"  );
    tree_->Branch("genpippx",    &genpippx   , "genpippx/D"  );
    tree_->Branch("genpippy",    &genpippy   , "genpippy/D"  );
    tree_->Branch("genpippz",    &genpippz   , "genpippz/D"  );
    tree_->Branch("decname",  &decname);
    tree_->Branch("istruemum",  &istruemum );
    tree_->Branch("istruemup",  &istruemup );
    tree_->Branch("istruepr",   &istruepr  );
    tree_->Branch("istruepi",   &istruepi  );
    tree_->Branch("istruelb",   &istruelb  );


  }


}

// ------------ method called once each job just after ending the event loop  ------------
void 
LambdaB::endJob() 
{

  fout_->cd();
  tree_->Write();

  for(int i = 0; i < kHistNameSize; i++) {
    histos[i]->Write();
    histos[i]->Delete();
  }
  fout_->Close();

  t_now_.Set();
  printf(" \n ---------- End Job ---------- \n" ) ;
  t_now_.Print();
  printf(" processed: %i \n selected: %i \n \
 duration: %i sec \n rate: %g evts/sec\n",
         n_processed_, n_selected_,
         t_now_.Convert() - t_begin_.Convert(),
         float(n_processed_)/(t_now_.Convert()-t_begin_.Convert()) );

}

// ------------ method called when starting to processes a run  ------------
void 
LambdaB::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
LambdaB::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
LambdaB::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
LambdaB::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LambdaB::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
LambdaB::clearVariables(){
  run = 0;
  event = 0;
  lumiblock = 0;
  nprivtx = 0;
  triggernames->clear();
  triggerprescales->clear();
  mumdcabs->clear();  mumdcabserr->clear();  mumpx->clear();   mumpy->clear();  mumpz->clear();
  mupdcabs->clear();  mupdcabserr->clear();  muppx->clear();   muppy->clear();  muppz->clear();
  mumutrkr->clear(); mumutrkz->clear();
  mumudca->clear();  mumuvtxcl->clear();   mumulsbs->clear();  mumulsbserr->clear();
  mumucosalphabs->clear();  mumucosalphabserr->clear();
  mumumass->clear(); mumumasserr->clear();
  mumisgoodmuon->clear();  mupisgoodmuon->clear();
  mumnpixhits->clear();  mupnpixhits->clear();  mumnpixlayers->clear();  mupnpixlayers->clear();
  mumntrkhits->clear();  mupntrkhits->clear();  mumntrklayers->clear();  mupntrklayers->clear();
  mumnormchi2->clear(); mupnormchi2->clear(); 
  mumdxyvtx->clear(); mupdxyvtx->clear();
  mumdzvtx->clear(); mupdzvtx->clear();  
  mumtriglastfilter->clear(); muptriglastfilter->clear();
  mumpt->clear(); muppt->clear();
  mumeta->clear(); mupeta->clear();

  trkchg->clear(); trkpx->clear(); trkpy->clear(); trkpz->clear(); trkpt->clear();
  trkdcabs->clear(); trkdcabserr->clear();

  lzpx->clear(); lzpy->clear(); lzpz->clear();
  lzvtxx->clear(); lzvtxy->clear(); lzvtxz->clear();

  prpx->clear(); prpy->clear(); prpz->clear();
  pipx->clear(); pipy->clear(); pipz->clear();

  lzmass->clear(); lzmasserr->clear();
  lzbarmass->clear(); lzbarmasserr->clear();

  nb = 0;

  lbpx->clear(); lbpxerr->clear(); lbpy->clear();  lbpyerr->clear();
  lbpz->clear(); lbpzerr->clear();

  lbmass->clear(); lbmasserr->clear();
  lbvtxcl->clear(); lbvtxx->clear(); lbvtxxerr->clear(); lbvtxy->clear();
  lbvtxyerr->clear();
  lbvtxz->clear(); lbvtxzerr->clear(); lbcosalphabs->clear(); lbcosalphabserr->clear();
  lbcosalphabs2d->clear(); lbcosalphabs2derr->clear();
  lblsbs->clear(); lblsbserr->clear();  lbctau->clear(); lbctauerr->clear();

  lbbarmass->clear(); lbbarmasserr->clear();

  if (IsMonteCarlo_) {

    genlbpx = 0;    genlbpy = 0;    genlbpz = 0;
    genlzpx = 0;  genlzpy = 0;  genlzpz = 0;
    genlzvtxx = 0; genlzvtxy = 0; genlzvtxz = 0;

    genprchg = 0;
    genprpx = 0;  genprpy = 0;  genprpz = 0;
    genpichg = 0;
    genpipx = 0;  genpipy = 0;  genpipz = 0;

    genmumpx = 0;  genmumpy = 0;  genmumpz = 0;
    genmuppx = 0;  genmuppy = 0;  genmuppz = 0;

    genpippx = 0; genpippy = 0; genpippz = 0;
    genpimpx = 0; genpimpy = 0; genpimpz = 0;

    decname = "";
    istruemum->clear(); istruemup->clear(); istruepr->clear();
    istruepi->clear(); istruelb->clear();


  }

}


void 
LambdaB::hltReport(const edm::Event& iEvent)
{

  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try {iEvent.getByLabel( TriggerResultsLabel_, hltTriggerResults ); }
  catch ( ... ) { edm::LogInfo("myHLT")
      << __LINE__ << " : couldn't get handle on HLT Trigger" ; }

  HLTConfigProvider hltConfig_;
  if (hltTriggerResults.isValid()) {
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);

    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++){

      // Only consider the triggered case.                                                                                                                                   
      if ((*hltTriggerResults)[itrig].accept() == 1){

        string triggername = triggerNames_.triggerName(itrig);
        int triggerprescale = hltConfig_.prescaleValue(itrig, triggername);

        // Loop over our interested HLT trigger names to find if this event contains.                                                                                        
        for (unsigned int it=0; it<TriggerNames_.size(); it++){
          if (triggername.find(TriggerNames_[it]) != string::npos) {
            // save the no versioned case                                                                                                                                    
            triggernames->push_back(TriggerNames_[it]);
            triggerprescales->push_back(triggerprescale);

          }}}}}

}


bool
LambdaB::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByLabel(BeamSpotLabel_, beamSpotHandle);

  if ( ! beamSpotHandle.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ;
    return false;
  }

  beamSpot_ = *beamSpotHandle;
  return true;
}


bool
LambdaB::hasPrimaryVertex(const edm::Event& iEvent)
{
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel(VertexLabel_, recVtxs);
  nprivtx = recVtxs->size();

  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtxs->begin();
       iVertex != recVtxs->end(); iVertex++) {
    primaryVertex_ = *(iVertex);
    if (primaryVertex_.isValid()) break;
  }

  if (!primaryVertex_.isValid()) return false;

  return true;
}



//define this as a plug-in
DEFINE_FWK_MODULE(LambdaB);
