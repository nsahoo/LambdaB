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
//*         Ntuplizer code for /\b -> /\(->p+pi) mu+ mu-             *
//*----   ---------------   ----------------   --------------   ---- *
//*   original author:  Niladribihari Sahoo,42 3-024,+41227662373,   *
//*      date created:  Tue Dec  8 23:43:26 CET 2015                 *
//*          modified:  dec 21 (fixed vertexing error)               *
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
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
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

  bool hasGoodLzVertex(const vector<reco::TrackRef>,
			   RefCountedKinematicTree &);

  bool hasGoodLzVertexMKC(const vector<reco::TrackRef>,
			      RefCountedKinematicTree &);


  bool hasGoodMuonDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaBs (const reco::TransientTrack, double &, double &);
  bool hasGoodTrackDcaPoint (const reco::TransientTrack, const GlobalPoint,
                             double, double &, double &);

  bool hasGoodLbMass(RefCountedKinematicTree, double &);  // LambdaB mass


  bool hasGoodLbVertex(const reco::TrackRef, const reco::TrackRef,   // LambdaB vertex
                       const vector<reco::TrackRef>, double &, double &, double &,      
		       RefCountedKinematicTree & , RefCountedKinematicTree & );


  bool hasGoodMuMuVertex (const reco::TransientTrack, const reco::TransientTrack,
                          reco::TransientTrack &, reco::TransientTrack &,
                          double &, double &, double &, double &, double &,
                          double &, double &, double &);

  bool hasGoodTrack(const edm::Event&, const pat::GenericParticle, double &);

  bool hasPrimaryVertex(const edm::Event &);

  void hltReport(const edm::Event&);

  bool matchMuonTrack (const edm::Event&, const reco::TrackRef);
  bool matchMuonTracks (const edm::Event&, const vector<reco::TrackRef>);
  bool matchLzTrack (const edm::Event&, const reco::TrackRef, double &, double &);     // added on dec 22 , 2015
  bool matchPrimaryVertexTracks ();

  void saveLbToLzMuMu(const RefCountedKinematicTree);
  void saveLbVertex(RefCountedKinematicTree);
  void saveLbCosAlpha(RefCountedKinematicTree);
  void saveLbCosAlpha2d(RefCountedKinematicTree);    
  void saveLbLsig(RefCountedKinematicTree);
  void saveLbCtau(RefCountedKinematicTree);

  void saveGenInfo(const edm::Event&);
  void saveLzVariables(RefCountedKinematicTree, reco::VertexCompositeCandidate);
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
  ParticleMass LzMass_;
  float LzMassErr_;
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
  edm::InputTag LambdaLabel_;
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



  vector<int> *trkchg; // +1 for pi+, -1 for pi-                                                                                                                
  vector<double> *trkpx, *trkpy, *trkpz, *trkpt;
  vector<double> *trkdcabs, *trkdcabserr;


  //--------------------
  // Lambda0  
  //--------------------                                                                                                                                       
  vector<int> *prchg;
  vector<double> *prpx, *prpy, *prpz ;
  vector<double>  *prd0, *prd0err ;
  vector<int> *pichg;
  vector<double> *pipx, *pipy, *pipz ;
  vector<double> *pid0, *pid0err ;
  vector<double> *lzpx, *lzpy, *lzpz;
  vector<double> *lzvtxx, *lzvtxy, *lzvtxz, *lzvtxcl, *lzlsbs, *lzlsbserr;

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
  LambdaLabel_(iConfig.getParameter<edm::InputTag>("LambdaLabel")),
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

  prchg(0),
  prpx(0), prpy(0), prpz(0),
  prd0(0), prd0err(0),
  pichg(0),
  pipx(0), pipy(0), pipz(0),
  pid0(0), pid0err(0),

  lzpx(0), lzpy(0), lzpz(0),
  lzvtxx(0), lzvtxy(0), lzvtxz(0),

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
  tree_->Branch("prchg", &prchg);
  tree_->Branch("prpx", &prpx);
  tree_->Branch("prpy", &prpy);
  tree_->Branch("prpz", &prpz);
  tree_->Branch("prd0", &prd0);
  tree_->Branch("prd0err", &prd0err);
  tree_->Branch("pichg", &pichg);
  tree_->Branch("pipx", &pipx);
  tree_->Branch("pipy", &pipy);
  tree_->Branch("pipz", &pipz);
  tree_->Branch("pid0", &pid0);
  tree_->Branch("pid0err", &pid0err);

  tree_->Branch("lzpx", &lzpx);
  tree_->Branch("lzpy", &lzpy);
  tree_->Branch("lzpz", &lzpz);
  tree_->Branch("lzvtxx", &lzvtxx);
  tree_->Branch("lzvtxy", &lzvtxy);
  tree_->Branch("lzvtxz", &lzvtxz);
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

  prchg->clear();
  prpx->clear(); prpy->clear(); prpz->clear();
  prd0->clear(); prd0err->clear();
  pichg->clear();
  pipx->clear(); pipy->clear(); pipz->clear();
  pid0->clear(); pid0err->clear();

  lzpx->clear(); lzpy->clear(); lzpz->clear();
  lzvtxx->clear(); lzvtxy->clear(); lzvtxz->clear();

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


bool
LambdaB::hasGoodTrack(const edm::Event& iEvent,
			   const pat::GenericParticle iTrack,
			   double & trk_pt)
{
  reco::TrackRef theTrackRef = iTrack.track();
  if ( theTrackRef.isNull() ) return false;

  // veto muon tracks                                                                                                                       
  if ( matchMuonTrack(iEvent, theTrackRef) ) return false;

  // check the track kinematics                                                                                                 
  trk_pt = theTrackRef->pt();

  if ( theTrackRef->pt() < TrkMinPt_ ) return false;

  return true;
}



//-----------------------------------------------------------------------------------
// main function to build/reconstruct LambdaB --> Lambda0 mu+ mu- candidate
//-----------------------------------------------------------------------------------
bool
LambdaB::buildLbToLzMuMu(const edm::Event& iEvent)
{

  // init variables                                                                                                                               
  edm::Handle< vector<pat::Muon> > patMuonHandle;
  iEvent.getByLabel(MuonLabel_, patMuonHandle);
  if( patMuonHandle->size() < 2 ) return false;

  edm::Handle<reco::VertexCompositeCandidateCollection> theLambdas;
  iEvent.getByLabel(LambdaLabel_, theLambdas);
  if ( theLambdas->size() <= 0) return false;

  //edm::Handle< vector<pat::GenericParticle> >thePATTrackHandle;
  //iEvent.getByLabel(TrackLabel_, thePATTrackHandle);

  bool passed;
  double DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr;
  double mumutrk_R, mumutrk_Z, DCAmumu;
  //double trk_R, trk_Z, trk_DCA;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  double MuMuLSBS, MuMuLSBSErr;
  double MuMuCosAlphaBS, MuMuCosAlphaBSErr;
  double lz_mass, lzbar_mass; 
  double lb_vtx_chisq, lb_vtx_cl, lb_mass, lbbar_mass;
  double DCALzTrkBS, DCALzTrkBSErr;
  vector<reco::TrackRef> LambdaDaughterTracks;
  RefCountedKinematicTree vertexFitTree, barVertexFitTree, LzvertexFitTree ;

  // --------------------                               
  // loop 1: mu-                                                                                                                          
  // --------------------                                                                                                                      
  for (vector<pat::Muon>::const_iterator iMuonM = patMuonHandle->begin();
       iMuonM != patMuonHandle->end(); iMuonM++){

    reco::TrackRef muTrackm = iMuonM->innerTrack();
    if ( muTrackm.isNull() ) continue;

    histos[h_mupt]->Fill(muTrackm->pt());
    histos[h_mueta]->Fill(muTrackm->eta());

    if ( (muTrackm->charge() != -1) ||
         (muTrackm->pt() < MuonMinPt_) ||
         (fabs(muTrackm->eta()) > MuonMaxEta_)) continue;

    // check mu- DCA to beam spot                                                                                                                                 
    const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle_));
    passed = hasGoodMuonDcaBs(muTrackmTT, DCAmumBS, DCAmumBSErr) ;
    histos[h_mumdcabs]->Fill(DCAmumBS);
    if ( ! passed ) continue;

    // --------------------- 
    // loop 2: mu+                                                                                                                      
    // ---------------------                                                                                                                                  
    for (vector<pat::Muon>::const_iterator iMuonP = patMuonHandle->begin();
         iMuonP != patMuonHandle->end(); iMuonP++){

      reco::TrackRef muTrackp = iMuonP->innerTrack();
      if ( muTrackp.isNull() ||
           (muTrackp->charge() != 1) ||
           (muTrackp->pt() < MuonMinPt_) ||
           (fabs(muTrackp->eta()) > MuonMaxEta_)) continue;

      // check mu+ DCA to beam spot                                                                                                             
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle_));
      passed = hasGoodMuonDcaBs(muTrackpTT, DCAmupBS, DCAmupBSErr);
      if ( ! passed ) continue;


      // check goodness of muons closest approach and the 3D-DCA                                                                                             
      // passed = hasGoodClosestApproachTracks(muTrackpTT, muTrackmTT,                                                                                              
      //                                            mumutrk_R, mumutrk_Z, DCAmumu);                                                                      
      if ( !calClosestApproachTracks(muTrackpTT, muTrackmTT,
                                     mumutrk_R, mumutrk_Z, DCAmumu)) continue;
      histos[h_mumutrkr]->Fill(mumutrk_R);
      histos[h_mumutrkz]->Fill(mumutrk_Z);
      histos[h_mumudca]->Fill(DCAmumu);
      // if ( !passed ) continue;                                                                                                                              
      if ( mumutrk_R > TrkMaxR_ ||
           mumutrk_Z > TrkMaxZ_ ||
           DCAmumu > MuMuMaxDca_ ) continue;

      // check dimuon vertex                                                                                                                            
      passed = hasGoodMuMuVertex(muTrackpTT, muTrackmTT, refitMupTT, refitMumTT,
                                 mu_mu_vtx_cl, mu_mu_pt,
                                 mu_mu_mass, mu_mu_mass_err,
                                 MuMuLSBS, MuMuLSBSErr,
                                 MuMuCosAlphaBS, MuMuCosAlphaBSErr);

      histos[h_mumuvtxcl]->Fill(mu_mu_vtx_cl);
      histos[h_mumupt]->Fill(mu_mu_pt);
      histos[h_mumumass]->Fill(mu_mu_mass);
      histos[h_mumulxybs]->Fill(MuMuLSBS/MuMuLSBSErr);
      histos[h_mumucosalphabs]->Fill(MuMuCosAlphaBS);
      if ( !passed) continue;

      //----------------
      // loop 3: /\0
      //----------------
      for ( reco::VertexCompositeCandidateCollection::const_iterator iLambda
	      = theLambdas->begin(); iLambda != theLambdas->end(); ++iLambda) {
	
	LambdaDaughterTracks.push_back((dynamic_cast<const
					reco::RecoChargedCandidate *>
					(iLambda->daughter(0)))->track());
	LambdaDaughterTracks.push_back((dynamic_cast<const
					reco::RecoChargedCandidate *>
					(iLambda->daughter(1)))->track());

	
	// check that the tracks of Lambda0(bar) are not *muons*
	if ( matchMuonTracks(iEvent, LambdaDaughterTracks) ) continue;





	  // fit /\b vertex  mu- mu+ /\(pi- p+)
          if ( ! hasGoodLbVertex(muTrackm, muTrackp, LambdaDaughterTracks,
                                 lb_vtx_chisq, lb_vtx_cl, lb_mass,
                                 vertexFitTree, LzvertexFitTree ) ) continue;

          histos[h_lbvtxchisq]->Fill(lb_vtx_chisq);
          histos[h_lbvtxcl]->Fill(lb_vtx_cl);

          if ( lb_vtx_cl < LbMinVtxCl_ || lb_mass < LbMinMass_ || lb_mass > LbMaxMass_ ) continue;


	  // fit /\bbar vertex mu- mu+ /\bar(pi+ p-)                                                                                                                
          if ( ! hasGoodLbVertex(muTrackm, muTrackp, LambdaDaughterTracks,
                                 lb_vtx_chisq, lb_vtx_cl, lbbar_mass,
                                 barVertexFitTree, LzvertexFitTree ) ) continue;

          if ( lb_vtx_cl < LbMinVtxCl_ || lbbar_mass < LbMinMass_ || lbbar_mass > LbMaxMass_ ) continue;


	  // need to check with primaryVertex tracks ???                                                                                                           

          nb++;

          // save the tree variables                                                                                            
          saveDimuVariables(DCAmumBS, DCAmumBSErr, DCAmupBS, DCAmupBSErr,
                            mumutrk_R, mumutrk_Z, DCAmumu, mu_mu_vtx_cl,
                            MuMuLSBS, MuMuLSBSErr,
                            MuMuCosAlphaBS, MuMuCosAlphaBSErr,
                            mu_mu_mass, mu_mu_mass_err);

	  ///saveLzVariables(LzvertexFitTree, *iLambda);
          saveSoftMuonVariables(*iMuonM, *iMuonP, muTrackm, muTrackp);
     
          trkdcabs->push_back(DCALzTrkBS);
          trkdcabserr->push_back(DCALzTrkBSErr);
          lzmass->push_back(lz_mass);
	  lzbarmass->push_back(lzbar_mass);
	  lbvtxcl->push_back(lb_vtx_cl);
	  

          saveLbToLzMuMu(vertexFitTree);
          saveLbVertex(vertexFitTree);
          saveLbCosAlpha(vertexFitTree);
          saveLbCosAlpha2d(vertexFitTree);
          saveLbLsig(vertexFitTree);
          saveLbCtau(vertexFitTree);



      }// close Lambda0 loop 
    } // close mu+ loop                                                                                                                                        
  } // close mu- loop                                                                                                                                         

  if ( nb > 0) {
    edm::LogInfo("myLb") << "Found " << nb << " /\b -> /\0 mu+ mu- ";
    return true;
  }
  return false;

}


void
LambdaB::calLS (double Vx, double Vy, double Vz,
                double Wx, double Wy, double Wz,
                double VxErr2, double VyErr2, double VzErr2,
                double VxyCov, double VxzCov, double VyzCov,
                double WxErr2, double WyErr2, double WzErr2,
                double WxyCov, double WxzCov, double WyzCov,
                double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt(  (Vx-Wx) * (Vx-Wx) * VxErr2 +
                        (Vy-Wy) * (Vy-Wy) * VyErr2 +
                        (Vz-Wz) * (Vz-Wz) * VzErr2 +

                        (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
                        (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
                        (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +

                        (Vx-Wx) * (Vx-Wx) * WxErr2 +
                        (Vy-Wy) * (Vy-Wy) * WyErr2 +
                        (Vz-Wz) * (Vz-Wz) * WzErr2 +

                        (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
                        (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
                        (Vy-Wy) * (Vz-Wz) * 2.*WyzCov ) / *deltaD;

  else *deltaDErr = 0.;

}


void
LambdaB::calCosAlpha (double Vx, double Vy, double Vz,
                      double Wx, double Wy, double Wz,
                      double VxErr2, double VyErr2, double VzErr2,
                      double VxyCov, double VxzCov, double VyzCov,
                      double WxErr2, double WyErr2, double WzErr2,
                      double WxyCov, double WxzCov, double WyzCov,
                      double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( (
			  (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			  (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			  (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			  (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			  (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			  (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /(Wnorm*Wnorm*Wnorm*Wnorm) +

              		  ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			  (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			  (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +

			  (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			  (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			  (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /(Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);

  }  else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }

}


void
LambdaB::calCosAlpha2d (double Vx, double Vy, double Vz,
                        double Wx, double Wy, double Wz,
                        double VxErr2, double VyErr2, double VzErr2,
                        double VxyCov, double VxzCov, double VyzCov,
                        double WxErr2, double WyErr2, double WzErr2,
                        double WxyCov, double WxzCov, double WyzCov,
                        double* cosAlpha2d, double* cosAlpha2dErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha2d = VdotW / (Vnorm * Wnorm);
    *cosAlpha2dErr = sqrt( (
			    (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			    (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			    (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			    (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			    (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			    (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			   (Wnorm*Wnorm*Wnorm*Wnorm) +

			   ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			    (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			    (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +

			    (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			    (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			    (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			   (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha2d = 0.;
    *cosAlpha2dErr = 0.;
  }
}



bool
LambdaB::hasGoodMuonDcaBs (const reco::TransientTrack muTrackTT,
                                double &muDcaBs, double &muDcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    muTrackTT.trajectoryStateClosestToPoint(
					    GlobalPoint(beamSpot_.position().x(),
							beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() )  return false;

  muDcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  muDcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(muDcaBs) > MuonMaxDcaBs_ )   return false;
  return true;
}


bool
LambdaB::hasGoodTrackDcaBs (const reco::TransientTrack TrackTT,
                                 double &DcaBs, double &DcaBsErr)
{
  TrajectoryStateClosestToPoint theDCAXBS =
    TrackTT.trajectoryStateClosestToPoint(
					  GlobalPoint(beamSpot_.position().x(),
						      beamSpot_.position().y(),beamSpot_.position().z()));

  if ( !theDCAXBS.isValid() )  return false;

  DcaBs = theDCAXBS.perigeeParameters().transverseImpactParameter();
  DcaBsErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  if ( fabs(DcaBs/DcaBsErr) < TrkMinDcaSigBs_ )   return false;
  return true;
}


bool
LambdaB::hasGoodTrackDcaPoint (const reco::TransientTrack track,
                               const GlobalPoint p,
                               double maxdca, double &dca, double &dcaerr)
{
  TrajectoryStateClosestToPoint theDCAX = track.trajectoryStateClosestToPoint(p);
  if ( !theDCAX.isValid() ) return false;

  dca = theDCAX.perigeeParameters().transverseImpactParameter();
  dcaerr = theDCAX.perigeeError().transverseImpactParameterError();
  if ( dca > maxdca ) return false;

  return true;
}

bool
LambdaB::calClosestApproachTracks (const reco::TransientTrack trackpTT,
                                        const reco::TransientTrack trackmTT,
                                        double & trk_R,
                                        double & trk_Z,
                                        double & trk_DCA)
{
  ClosestApproachInRPhi ClosestApp;
  ClosestApp.calculate(trackpTT.initialFreeState(),
                       trackmTT.initialFreeState());
  if (! ClosestApp.status() )  return false ;

  GlobalPoint XingPoint = ClosestApp.crossingPoint();

  trk_R = sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y());
  trk_Z = fabs(XingPoint.z());

  // if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) >                                                                            
  //      TrkMaxR_) || (fabs(XingPoint.z()) > TrkMaxZ_))  return false;                                                                              

  trk_DCA = ClosestApp.distance();
  // if (DCAmumu > MuMuMaxDca_) return false;                                                                                                                 

  return true;
}


bool
LambdaB::hasGoodMuMuVertex (const reco::TransientTrack muTrackpTT,
                            const reco::TransientTrack muTrackmTT,
                            reco::TransientTrack &refitMupTT,
                            reco::TransientTrack &refitMumTT,
                            double & mu_mu_vtx_cl, double & mu_mu_pt,
                            double & mu_mu_mass, double & mu_mu_mass_err,
                            double & MuMuLSBS, double & MuMuLSBSErr,
                            double & MuMuCosAlphaBS,
                            double & MuMuCosAlphaBSErr)
{
  KinematicParticleFactoryFromTransientTrack partFactory;
  KinematicParticleVertexFitter PartVtxFitter;

  vector<RefCountedKinematicParticle> muonParticles;
  double chi = 0.;
  double ndf = 0.;
  muonParticles.push_back(partFactory.particle(muTrackmTT, MuonMass_,chi,ndf,MuonMassErr_));
  muonParticles.push_back(partFactory.particle(muTrackpTT, MuonMass_,chi,ndf,MuonMassErr_));

  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);

  if ( !mumuVertexFitTree->isValid())  return false;

  mumuVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
  RefCountedKinematicVertex mumu_KV = mumuVertexFitTree->currentDecayVertex();

  if ( !mumu_KV->vertexIsValid()) return false;

  mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
                             int(rint(mumu_KV->degreesOfFreedom())));

  if (mu_mu_vtx_cl < MuMuMinVtxCl_)  return false;

  // extract the re-fitted tracks                                                                                                               
  mumuVertexFitTree->movePointerToTheTop();

  mumuVertexFitTree->movePointerToTheFirstChild();
  RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
  refitMumTT = refitMum->refittedTransientTrack();

  mumuVertexFitTree->movePointerToTheNextChild();
  RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
  refitMupTT = refitMup->refittedTransientTrack();

  TLorentzVector mymum, mymup, mydimu;

  mymum.SetXYZM(refitMumTT.track().momentum().x(),
                refitMumTT.track().momentum().y(),
                refitMumTT.track().momentum().z(), MuonMass_);

  mymup.SetXYZM(refitMupTT.track().momentum().x(),
                refitMupTT.track().momentum().y(),
                refitMupTT.track().momentum().z(), MuonMass_);

  mydimu = mymum + mymup;
  mu_mu_pt = mydimu.Perp();

  mu_mu_mass = mumu_KP->currentState().mass();
  mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
                        matrix()(6,6));

  if ((mu_mu_pt < MuMuMinPt_) || (mu_mu_mass < MuMuMinInvMass_) ||
      (mu_mu_mass > MuMuMaxInvMass_))  return false;

  // compute the distance between mumu vtx and beam spot                                                                                                         
  calLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
         beamSpot_.position().x(),beamSpot_.position().y(),0.0,
         mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
         mumu_KV->error().matrix()(0,1),0.0,0.0,
         beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
         beamSpot_.covariance()(0,1),0.0,0.0,
         &MuMuLSBS,&MuMuLSBSErr);

  if (MuMuLSBS/MuMuLSBSErr < MuMuMinLxySigmaBs_)  return false;

  calCosAlpha(mumu_KP->currentState().globalMomentum().x(),
              mumu_KP->currentState().globalMomentum().y(),
              0.0,
              mumu_KV->position().x() - beamSpot_.position().x(),
              mumu_KV->position().y() - beamSpot_.position().y(),
              0.0,
              mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
              mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
              0.0,
              mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
              0.0,
              0.0,
              mumu_KV->error().cxx() + beamSpot_.covariance()(0,0),
              mumu_KV->error().cyy() + beamSpot_.covariance()(1,1),
              0.0,
              mumu_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
              0.0,
              0.0,
              &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);

  if (MuMuCosAlphaBS < MuMuMinCosAlphaBs_)  return false;

  return true;

}


bool
LambdaB::matchMuonTrack (const edm::Event& iEvent,
                         const reco::TrackRef theTrackRef)
{
  if ( theTrackRef.isNull() ) return false;

  edm::Handle< vector<pat::Muon> > thePATMuonHandle;
  iEvent.getByLabel(MuonLabel_, thePATMuonHandle);

  reco::TrackRef muTrackRef;
  for (vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin();
       iMuon != thePATMuonHandle->end(); iMuon++){

    muTrackRef = iMuon->innerTrack();
    if ( muTrackRef.isNull() ) continue;

    if (muTrackRef == theTrackRef) return true;
  }

  return false;
}


bool
LambdaB::matchMuonTracks (const edm::Event& iEvent,
                          const vector<reco::TrackRef> theTracks)
{
  reco::TrackRef theTrackRef;
  for(unsigned int j = 0; j < theTracks.size(); ++j) {
    theTrackRef = theTracks[j];
    if ( matchMuonTrack(iEvent, theTrackRef) ) return true;
  }
  return false;
}


// added on dec 22, 2015
bool
LambdaB::matchLzTrack (const edm::Event& iEvent,
		       const reco::TrackRef theTrackRef, double & dau1_pt, double & dau2_pt)
{
  if ( theTrackRef.isNull() ) return false;

  edm::Handle<reco::VertexCompositeCandidateCollection> theLambdas;
  iEvent.getByLabel(LambdaLabel_, theLambdas);

  reco::TrackRef theDau1TrackRef;
  reco::TrackRef theDau2TrackRef;

  for ( reco::VertexCompositeCandidateCollection::const_iterator iLambda
	  = theLambdas->begin(); iLambda != theLambdas->end(); ++iLambda) {

    theDau1TrackRef = (*(dynamic_cast<const reco::RecoChargedCandidate *>
			 (iLambda->daughter(0)))).track();

    dau1_pt = theDau1TrackRef->pt();

    if ( ! theDau1TrackRef.isNull() && theDau1TrackRef == theTrackRef) return true;


    theDau2TrackRef = (*(dynamic_cast<const reco::RecoChargedCandidate *>
			 (iLambda->daughter(1)))).track();

    dau2_pt = theDau2TrackRef->pt();

    if ( ! theDau2TrackRef.isNull() && theDau2TrackRef == theTrackRef) return true;
  }

  return false;

}


//------------------------
// added on 13 dec 2015 //
//------------------------
bool
LambdaB::hasGoodLzVertex(const vector<reco::TrackRef> theDaughterTracks,
				  RefCountedKinematicTree &lzVertexFitTree)
{
  reco::TransientTrack dau1TT(theDaughterTracks[0], &(*bFieldHandle_) );
  reco::TransientTrack dau2TT(theDaughterTracks[1], &(*bFieldHandle_) );

  KinematicParticleFactoryFromTransientTrack pFactory;

  float chi = 0.;
  float ndf = 0.;
  vector<RefCountedKinematicParticle> lzdauParticles;

  // assign proton mass to the track w/ higher momentum
  if ( dau1TT.track.momentum() > dau2TT.track.momentum() ) {
    lzdauParticles.push_back(pFactory.particle(dau1TT,ProtonMass_,chi,ndf,ProtonMassErr_));
    lzdauParticles.push_back(pFactory.particle(dau2TT,  PionMass_,chi,ndf,  PionMassErr_));
  } else{
    lzdauParticles.push_back(pFactory.particle(dau2TT,ProtonMass_,chi,ndf,ProtonMassErr_));
    lzdauParticles.push_back(pFactory.particle(dau1TT,  PionMass_,chi,ndf,  PionMassErr_));
  }

  KinematicParticleVertexFitter fitter;
  lzVertexFitTree = fitter.fit(lzdauParticles);
  if (!lzVertexFitTree->isValid()) return false;

  return true;
}


bool
LambdaB::hasGoodLzVertexMKC(const vector<reco::TrackRef> theDaughterTracks,
				     RefCountedKinematicTree &lzVertexFitTree)
{
  if ( !hasGoodLzVertex(theDaughterTracks, lzVertexFitTree) ) return false;
  KinematicParticleFitter csFitterLz;
  KinematicConstraint * lz_c = new MassKinematicConstraint(LzMass_,
							   LzMassErr_);
  lzVertexFitTree = csFitterLz.fit(lz_c,lzVertexFitTree);

  delete lz_c;
  if (!lzVertexFitTree->isValid()) return false;
  return true;
}



bool
LambdaB::hasGoodLbVertex(const reco::TrackRef mu1Track,
			 const reco::TrackRef mu2Track,
			 const vector<reco::TrackRef> LambdaDaughterTracks,
			 double & lb_vtx_chisq, double & lb_vtx_cl,
			 double & lb_mass,
			 RefCountedKinematicTree & vertexFitTree, RefCountedKinematicTree & LzVertexFitTree )

{

  if ( ! hasGoodLzVertexMKC(LambdaDaughterTracks, LzVertexFitTree) )
    return false;

  LzVertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle lz_KP = LzVertexFitTree->currentParticle();

  KinematicParticleFactoryFromTransientTrack pFactory;

  reco::TransientTrack mu1TT(mu1Track, &(*bFieldHandle_) );
  reco::TransientTrack mu2TT(mu2Track, &(*bFieldHandle_) );
  //reco::TransientTrack pionTT(pionTrack, &(*bFieldHandle_) );
  //reco::TransientTrack protonTT(protonTrack, &(*bFieldHandle_) );

  float chi = 0.;
  float ndf = 0.;

  // /\b -> mu+ mu- /\0 (p+ pi-)                                                                                                                        
  vector<RefCountedKinematicParticle> vFitMCParticles;
  vFitMCParticles.push_back(pFactory.particle(   mu1TT,   MuonMass_, chi, ndf,   MuonMassErr_));
  vFitMCParticles.push_back(pFactory.particle(   mu2TT,   MuonMass_, chi, ndf,   MuonMassErr_));
  //vFitMCParticles.push_back(pFactory.particle(  pionTT,   PionMass_, chi, ndf,   PionMassErr_));
  //vFitMCParticles.push_back(pFactory.particle(protonTT, ProtonMass_, chi, ndf, ProtonMassErr_));
  vFitMCParticles.push_back(lz_KP);


  KinematicParticleVertexFitter fitter;
  vertexFitTree = fitter.fit(vFitMCParticles);
  if (!vertexFitTree->isValid()) return false;

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex LbDecayVertexMC = vertexFitTree->currentDecayVertex();
  if ( !LbDecayVertexMC->vertexIsValid()) return false;

  lb_vtx_chisq = LbDecayVertexMC->chiSquared();

  if ( LbDecayVertexMC->chiSquared()<0
       || LbDecayVertexMC->chiSquared()>1000 ) return false;

  RefCountedKinematicVertex lb_KV = vertexFitTree->currentDecayVertex();
  lb_vtx_cl = ChiSquaredProbability((double)(lb_KV->chiSquared()),
                                    (double)(lb_KV->degreesOfFreedom()));

  if ( lb_vtx_cl < LbMinVtxCl_ ) return false;

  RefCountedKinematicParticle lb_KP = vertexFitTree->currentParticle();
  lb_mass = lb_KP->currentState().mass();

  return true;
}


//---------------------------------
//   need further fixing ??
//---------------------------------

void
LambdaB::saveLbToLzMuMu(const RefCountedKinematicTree vertexFitTree){

  vertexFitTree->movePointerToTheTop(); // /\b --> /\(p+ pi-) mu+ mu-                                                                                         
  RefCountedKinematicParticle lb_KP = vertexFitTree->currentParticle();

  lbpx->push_back(lb_KP->currentState().globalMomentum().x());
  lbpxerr->push_back( sqrt( lb_KP->currentState().kinematicParametersError().matrix()(3,3) ) );
  lbpy->push_back(lb_KP->currentState().globalMomentum().y());
  lbpyerr->push_back( sqrt( lb_KP->currentState().kinematicParametersError().matrix()(4,4) ) );
  lbpz->push_back(lb_KP->currentState().globalMomentum().z());
  lbpzerr->push_back( sqrt( lb_KP->currentState().kinematicParametersError().matrix()(5,5) ) );
  lbmass->push_back(lb_KP->currentState().mass());
  lbmasserr->push_back( sqrt( lb_KP->currentState().kinematicParametersError().matrix()(6,6) ) );

  vertexFitTree->movePointerToTheFirstChild(); // mu1                                                                                                      
  RefCountedKinematicParticle mu1_KP = vertexFitTree->currentParticle();
  vertexFitTree->movePointerToTheNextChild();  // mu2                                                                                                         
  RefCountedKinematicParticle mu2_KP = vertexFitTree->currentParticle();

  RefCountedKinematicParticle mup_KP, mum_KP ;

  if ( mu1_KP->currentState().particleCharge() > 0 ) mup_KP = mu1_KP;
  if ( mu1_KP->currentState().particleCharge() < 0 ) mum_KP = mu1_KP;
  if ( mu2_KP->currentState().particleCharge() > 0 ) mup_KP = mu2_KP;
  if ( mu2_KP->currentState().particleCharge() < 0 ) mum_KP = mu2_KP;

  muppx->push_back(mup_KP->currentState().globalMomentum().x());
  muppy->push_back(mup_KP->currentState().globalMomentum().y());
  muppz->push_back(mup_KP->currentState().globalMomentum().z());

  mumpx->push_back(mum_KP->currentState().globalMomentum().x());
  mumpy->push_back(mum_KP->currentState().globalMomentum().y());
  mumpz->push_back(mum_KP->currentState().globalMomentum().z());


  vertexFitTree->movePointerToTheNextChild();  // pion track                                                                   
  RefCountedKinematicParticle pi_KP = vertexFitTree->currentParticle();
  pichg->push_back(pi_KP->currentState().particleCharge());
  pipx->push_back(pi_KP->currentState().globalMomentum().x());
  pipy->push_back(pi_KP->currentState().globalMomentum().y());
  pipz->push_back(pi_KP->currentState().globalMomentum().z());

  vertexFitTree->movePointerToTheNextChild();  // proton track
  RefCountedKinematicParticle pr_KP = vertexFitTree->currentParticle();
  prchg->push_back(pr_KP->currentState().particleCharge());
  prpx->push_back(pr_KP->currentState().globalMomentum().x());
  prpy->push_back(pr_KP->currentState().globalMomentum().y());
  prpz->push_back(pr_KP->currentState().globalMomentum().z());



}

void
LambdaB::saveLbVertex(RefCountedKinematicTree vertexFitTree){
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex lb_KV = vertexFitTree->currentDecayVertex();
  lbvtxx->push_back((*lb_KV).position().x());
  lbvtxxerr->push_back(sqrt( abs(lb_KV->error().cxx()) ));
  lbvtxy->push_back((*lb_KV).position().y());
  lbvtxyerr->push_back(sqrt( abs(lb_KV->error().cyy()) ));
  lbvtxz->push_back((*lb_KV).position().z());
  lbvtxzerr->push_back(sqrt( abs(lb_KV->error().czz()) ));

}

void 
LambdaB::saveLbCosAlpha(RefCountedKinematicTree vertexFitTree)
{
  // alpha is the angle in the transverse plane between the B0 momentum                                                                             
  // and the seperation between the B0 vertex and the beamspot                                                                                           

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle lb_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   lb_KV = vertexFitTree->currentDecayVertex();

  double cosAlphaBS, cosAlphaBSErr;
  calCosAlpha(lb_KP->currentState().globalMomentum().x(),
              lb_KP->currentState().globalMomentum().y(),
              lb_KP->currentState().globalMomentum().z(),
              lb_KV->position().x() - beamSpot_.position().x(),
              lb_KV->position().y() - beamSpot_.position().y(),
              lb_KV->position().z() - beamSpot_.position().z(),
              lb_KP->currentState().kinematicParametersError().matrix()(3,3),
              lb_KP->currentState().kinematicParametersError().matrix()(4,4),
              lb_KP->currentState().kinematicParametersError().matrix()(5,5),
              lb_KP->currentState().kinematicParametersError().matrix()(3,4),
              lb_KP->currentState().kinematicParametersError().matrix()(3,5),
              lb_KP->currentState().kinematicParametersError().matrix()(4,5),
              lb_KV->error().cxx() + beamSpot_.covariance()(0,0),
              lb_KV->error().cyy() + beamSpot_.covariance()(1,1),
              lb_KV->error().czz() + beamSpot_.covariance()(2,2),
              lb_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),
              lb_KV->error().matrix()(0,2) + beamSpot_.covariance()(0,2),
              lb_KV->error().matrix()(1,2) + beamSpot_.covariance()(1,2),
              &cosAlphaBS,&cosAlphaBSErr);


  lbcosalphabs->push_back(cosAlphaBS);
  lbcosalphabserr->push_back(cosAlphaBSErr);


}

void
LambdaB::saveLbCosAlpha2d(RefCountedKinematicTree vertexFitTree)
{

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle lb_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   lb_KV = vertexFitTree->currentDecayVertex();

  double cosAlphaBS2d, cosAlphaBS2dErr;
  calCosAlpha2d(lb_KP->currentState().globalMomentum().x(),
                lb_KP->currentState().globalMomentum().y(),0.0,                   
                lb_KV->position().x() - beamSpot_.position().x(),
                lb_KV->position().y() - beamSpot_.position().y(),0.0,
                lb_KP->currentState().kinematicParametersError().matrix()(3,3),
                lb_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                lb_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                lb_KV->error().cxx() + beamSpot_.covariance()(0,0),
                lb_KV->error().cyy() + beamSpot_.covariance()(1,1),0.0,
                lb_KV->error().matrix()(0,1) + beamSpot_.covariance()(0,1),0.0,0.0,
                &cosAlphaBS2d,&cosAlphaBS2dErr);


  lbcosalphabs2d->push_back(cosAlphaBS2d);
  lbcosalphabs2derr->push_back(cosAlphaBS2dErr);

}

bool
LambdaB::matchPrimaryVertexTracks ()
{
  vector<reco::TransientTrack> vertexTracks;
  for (vector<reco::TrackBaseRef>::const_iterator iTrack =
         primaryVertex_.tracks_begin();
       iTrack != primaryVertex_.tracks_end(); iTrack++){
    reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
  }

  return false;

}

void 
LambdaB::saveLbLsig(RefCountedKinematicTree vertexFitTree)
{
  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex lb_KV = vertexFitTree->currentDecayVertex();
  double LSBS, LSBSErr;

  calLS (lb_KV->position().x(), lb_KV->position().y(), 0.0,
         beamSpot_.position().x(), beamSpot_.position().y(), 0.0,
         lb_KV->error().cxx(), lb_KV->error().cyy(), 0.0,
         lb_KV->error().matrix()(0,1), 0.0, 0.0,
         beamSpot_.covariance()(0,0), beamSpot_.covariance()(1,1), 0.0,
         beamSpot_.covariance()(0,1), 0.0, 0.0,
         &LSBS,&LSBSErr);

  lblsbs->push_back(LSBS);
  lblsbserr->push_back(LSBSErr);

}


void
LambdaB::calCtau(RefCountedKinematicTree vertexFitTree,
                      double &lbctau, double &lbctauerr)
{
  //calculate ctau = (mB*(Bvtx-Pvtx)*pB)/(|pB|**2)                                                                                                     

  vertexFitTree->movePointerToTheTop();
  RefCountedKinematicParticle lb_KP = vertexFitTree->currentParticle();
  RefCountedKinematicVertex   lb_KV = vertexFitTree->currentDecayVertex();

  double betagamma = (lb_KP->currentState().globalMomentum().mag()/LbMass_);

  // calculate ctau error. Momentum error is negligible compared to                                                                                       
  // the vertex errors, so don't worry about it                                                                                                         

  GlobalPoint BVP = GlobalPoint( lb_KV->position() );
  GlobalPoint PVP = GlobalPoint( primaryVertex_.position().x(),
                                 primaryVertex_.position().y(),
                                 primaryVertex_.position().z() );
  GlobalVector sep3D = BVP-PVP;
  GlobalVector pBV = lb_KP->currentState().globalMomentum();
  lbctau = (LbMass_* (sep3D.dot(pBV)))/(pBV.dot(pBV));

  GlobalError BVE = lb_KV->error();
  GlobalError PVE = GlobalError( primaryVertex_.error() );
  VertexDistance3D theVertexDistance3D;
  Measurement1D TheMeasurement = theVertexDistance3D.distance( VertexState(BVP, BVE), VertexState(PVP, PVE) );
  double myError = TheMeasurement.error();

  //  ctau is defined by the portion of the flight distance along                                                                  
  //  the compoenent of the B momementum, so only consider the error                                                                                      
  //  of that component, too, which is accomplished by scaling by                                                                                        
  //  ((VB-VP)(dot)PB)/|VB-VP|*|PB|                                                                                                                       

  double scale = abs( (sep3D.dot(pBV))/(sep3D.mag()*pBV.mag()) );
  lbctauerr =  (myError*scale)/betagamma;

}

double
LambdaB::calEta (double Px, double Py, double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double
LambdaB::calPhi (double Px, double Py, double Pz)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double
LambdaB::calEtaPhiDistance (double Px1, double Py1, double Pz1,
				 double Px2, double Py2, double Pz2)
{
  double phi1 = calPhi (Px1,Py1,Pz1);
  double eta1 = calEta (Px1,Py1,Pz1);
  double phi2 = calPhi (Px2,Py2,Pz2);
  double eta2 = calEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}


void 
LambdaB::saveLbCtau(RefCountedKinematicTree vertexFitTree)
{
  double lbctau_temp, lbctauerr_temp;
  calCtau(vertexFitTree, lbctau_temp, lbctauerr_temp);
  lbctau->push_back(lbctau_temp);
  lbctauerr->push_back(lbctauerr_temp);
}

void
LambdaB::saveGenInfo(const edm::Event& iEvent){
  edm::Handle<reco::GenParticleCollection> genparticles;
  iEvent.getByLabel(GenParticlesLabel_, genparticles );

  // loop over all gen particles                                                                                                                      
  for( size_t i = 0; i < genparticles->size(); ++i ) {
    const reco::GenParticle & b = (*genparticles)[i];

    // only select /\b candidate                                                                                                                 
    if ( abs(b.pdgId()) != LAMBDAB_PDG_ID ) continue;

    int imum(-1), imup(-1), ilambz(-1), ipr(-1), ipi(-1);
    int ijpsi(-1), ipsi2s(-1);

    // loop over all /\b daughters                                                                                                       
    for ( size_t j = 0; j < b.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(b.daughter(j));

      if (dau.pdgId() == MUONMINUS_PDG_ID) imum = j;
      if (dau.pdgId() == -MUONMINUS_PDG_ID) imup = j;
      if (abs(dau.pdgId()) == LAMBDA_PDG_ID) ilambz = j;
      if (dau.pdgId() == JPSI_PDG_ID ) ijpsi = j;
      if (dau.pdgId() == PSI2S_PDG_ID ) ipsi2s = j;

    }

    if ( ilambz == -1 ) continue;

    const reco::Candidate & lambz = *(b.daughter(ilambz));

    for ( size_t j = 0; j < lambz.numberOfDaughters(); ++j){
      const reco::Candidate  &dau = *(lambz.daughter(j));
      if (abs(dau.pdgId()) == PROTON_PDG_ID) ipr = j;
      if (abs(dau.pdgId()) == PIONPLUS_PDG_ID) ipi = j;

    }

    if (ipr == -1 || ipi == -1) continue;


    // store the /\b and /\0 vars                                                                                                                        
    const reco::Candidate & pr = *(lambz.daughter(ipr));
    const reco::Candidate & pi = *(lambz.daughter(ipi));


    const reco::Candidate *mum = NULL;
    const reco::Candidate *mup = NULL;


    //--------------------
    // /\b -> /\ mu mu    
    //--------------------                                                                                                                       
    if (imum != -1 && imup != -1) {
      // cout << "Found GEN /\b-> /\ mu mu " << endl;                                                                                                    
      mum = b.daughter(imum);
      mup = b.daughter(imup);
      decname = "LbToLzMuMu";
    }

    //----------------------------
    // /\b -> /\ J/psi(->mu mu)   
    //----------------------------                                                                                                            
    else if ( ijpsi != -1 ) {
      // cout << "Found GEN /\b --> /\ J/psi " << endl;                                                                                                  
      const reco::Candidate & jpsi = *(b.daughter(ijpsi));
      for ( size_t j = 0; j < jpsi.numberOfDaughters(); ++j){
        const reco::Candidate  &dau = *(jpsi.daughter(j));
        if ( dau.pdgId() == MUONMINUS_PDG_ID) imum = j;
        if ( dau.pdgId() == -MUONMINUS_PDG_ID) imup = j;
      }
      if (imum != -1 && imup != -1) {
        mum = jpsi.daughter(imum);
        mup = jpsi.daughter(imup);
        decname = "LbToLzJPsi";
      }
    }

    //------------------------------
    // /\b -> /\ psi(2S)(->mu mu)   
    //------------------------------                                                                                                                    
    else if ( ipsi2s != -1) {
      // cout << "Found GEN /\b --> /\ psi(2S) " << endl;                                                                                                 
      const reco::Candidate & psi2s = *(b.daughter(ipsi2s));
      for ( size_t j = 0; j < psi2s.numberOfDaughters(); ++j){
        const reco::Candidate  &dau = *(psi2s.daughter(j));
        if ( dau.pdgId() == MUONMINUS_PDG_ID) imum = j;
        if ( dau.pdgId() == -MUONMINUS_PDG_ID) imup = j;
      }
      if (imum != -1 && imup != -1) {
        mum = psi2s.daughter(imum);
        mup = psi2s.daughter(imup);
        decname = "LbToLzPsi2S";
      }
    }

    if ( mum == NULL || mup == NULL) continue;

    // save gen info                                                                                                                                     
    genlbpx = b.px();
    genlbpy = b.py();
    genlbpz = b.pz();

    genlzpx = lambz.px();
    genlzpy = lambz.py();
    genlzpz = lambz.pz();

    genlzvtxx = lambz.vx();
    genlzvtxy = lambz.vy();
    genlzvtxz = lambz.vz();

    genprchg = pr.charge();
    genprpx = pr.px();
    genprpy = pr.py();
    genprpz = pr.pz();

    genpichg = pi.charge();
    genpipx = pi.px();
    genpipy = pi.py();
    genpipz = pi.pz();

    genmumpx = mum->px();
    genmumpy = mum->py();
    genmumpz = mum->pz();

    genmuppx = mup->px();
    genmuppy = mup->py();
    genmuppz = mup->pz();
  }
}


// added on dec 21, 2015

void 
LambdaB::saveLzVariables(RefCountedKinematicTree LzvertexFitTree,
			     reco::VertexCompositeCandidate iLambda)
{

  LzvertexFitTree->movePointerToTheTop();
  RefCountedKinematicVertex lz_vertex = LzvertexFitTree->currentDecayVertex();

  lzvtxcl->push_back( ChiSquaredProbability(
    (double)(lz_vertex->chiSquared()),
    (double)(lz_vertex->degreesOfFreedom())) );

  lzvtxx->push_back(lz_vertex->position().x());
  lzvtxy->push_back(lz_vertex->position().y());
  lzvtxz->push_back(lz_vertex->position().z());

  // Lambda0 vertex distance to the beam spot
  double LzLSBS, LzLSBSErr;

  calLS (lz_vertex->position().x(),lz_vertex->position().y(),0.0,
	 beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	 lz_vertex->error().cxx(),lz_vertex->error().cyy(),0.0,
	 lz_vertex->error().matrix()(0,1),0.0,0.0,
	 beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	 beamSpot_.covariance()(0,1),0.0,0.0,
	 &LzLSBS,&LzLSBSErr);

  lzlsbs->push_back(LzLSBS);
  lzlsbserr->push_back(LzLSBSErr);

  LzvertexFitTree->movePointerToTheFirstChild(); // Lambda0 proton
  RefCountedKinematicParticle LzPr1 = LzvertexFitTree->currentParticle();

  LzvertexFitTree->movePointerToTheNextChild(); // Lambda0 pion
  RefCountedKinematicParticle LzPi2 = LzvertexFitTree->currentParticle();


  KinematicParameters LzPr1KP = LzPr1->currentState().kinematicParameters();
  KinematicParameters LzPi2KP = LzPi2->currentState().kinematicParameters();
  KinematicParameters LzPrpKP;
  KinematicParameters LzPrmKP;
  KinematicParameters LzPipKP;
  KinematicParameters LzPimKP;


  // xcheck again !!
  if ( LzPr1->currentState().particleCharge() > 0 ) LzPrpKP = LzPr1KP;
  if ( LzPr1->currentState().particleCharge() < 0 ) LzPrmKP = LzPr1KP;
  if ( LzPi2->currentState().particleCharge() > 0 ) LzPipKP = LzPi2KP;
  if ( LzPi2->currentState().particleCharge() < 0 ) LzPimKP = LzPi2KP;


  pippx->push_back(ksPipKP.momentum().x());
  pippy->push_back(ksPipKP.momentum().y());
  pippz->push_back(ksPipKP.momentum().z());

  pimpx->push_back(ksPimKP.momentum().x());
  pimpy->push_back(ksPimKP.momentum().y());
  pimpz->push_back(ksPimKP.momentum().z());



  if ( iKshort.daughter(0)->charge() < 0) {
    pimd0->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
		      (iKshort.daughter(0)))->track()->d0());
    pimd0err->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
			 (iKshort.daughter(0)))->track()->d0Error());
    pipd0->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
		      (iKshort.daughter(1)))->track()->d0());
    pipd0err->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
			 (iKshort.daughter(1)))->track()->d0Error());

  } else {
    pimd0->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
		      (iKshort.daughter(1)))->track()->d0());
    pimd0err->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
			 (iKshort.daughter(1)))->track()->d0Error());

    pipd0->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
		      (iKshort.daughter(0)))->track()->d0());
    pipd0err->push_back((dynamic_cast<const reco::RecoChargedCandidate *>
			 (iKshort.daughter(0)))->track()->d0Error());
  }

}


void
LambdaB::saveSoftMuonVariables(pat::Muon iMuonM, pat::Muon iMuonP,
                                    reco::TrackRef muTrackm, reco::TrackRef muTrackp)
{

  mumisgoodmuon->push_back(muon::isGoodMuon(iMuonM, muon::TMOneStationTight));
  mupisgoodmuon->push_back(muon::isGoodMuon(iMuonP, muon::TMOneStationTight));
  mumnpixhits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
  mupnpixhits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
  mumnpixlayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());
  mupnpixlayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());

  mumntrkhits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
  mupntrkhits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
  mumntrklayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
  mupntrklayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());

  mumnormchi2->push_back(muTrackm->normalizedChi2());
  mupnormchi2->push_back(muTrackp->normalizedChi2());

  mumdxyvtx->push_back(muTrackm->dxy(primaryVertex_.position()));
  mupdxyvtx->push_back(muTrackp->dxy(primaryVertex_.position()));

  mumdzvtx->push_back(muTrackm->dz(primaryVertex_.position()));
  mupdzvtx->push_back(muTrackp->dz(primaryVertex_.position()));

  mumpt->push_back(muTrackm->pt());
  muppt->push_back(muTrackp->pt());

  mumeta->push_back(muTrackm->eta());
  mupeta->push_back(muTrackp->eta());

}


void
LambdaB::saveDimuVariables(double DCAmumBS, double DCAmumBSErr,
                                double DCAmupBS, double DCAmupBSErr,
                                double mumutrk_R, double mumutrk_Z,
                                double DCAmumu,  double mu_mu_vtx_cl,
                                double MuMuLSBS, double MuMuLSBSErr,
                                double MuMuCosAlphaBS, double MuMuCosAlphaBSErr,
                                double mu_mu_mass, double mu_mu_mass_err)

{
  mumdcabs->push_back(DCAmumBS);
  mumdcabserr->push_back(DCAmumBSErr);

  mupdcabs->push_back(DCAmupBS);
  mupdcabserr->push_back(DCAmupBSErr);

  mumutrkr->push_back(mumutrk_R);
  mumutrkz->push_back(mumutrk_Z);
  mumudca->push_back(DCAmumu);
  mumuvtxcl->push_back(mu_mu_vtx_cl);
  mumulsbs->push_back(MuMuLSBS);
  mumulsbserr->push_back(MuMuLSBSErr);
  mumucosalphabs->push_back(MuMuCosAlphaBS);
  mumucosalphabserr->push_back(MuMuCosAlphaBSErr);

  mumumass->push_back(mu_mu_mass);
  mumumasserr->push_back(mu_mu_mass_err);
}

void
LambdaB::saveMuonTriggerMatches(const pat::Muon iMuonM, const pat::Muon iMuonP)
{
  string mum_matched_lastfilter_name = "";
  string mup_matched_lastfilter_name = "";

  for(vector<string>::iterator it = triggernames->begin();
      it != triggernames->end(); ++it) {

    string hltLastFilterName = mapTriggerToLastFilter_[*it] ;

    const pat::TriggerObjectStandAloneCollection mumHLTMatches
      = iMuonM.triggerObjectMatchesByFilter( hltLastFilterName );
    const pat::TriggerObjectStandAloneCollection mupHLTMatches
      = iMuonP.triggerObjectMatchesByFilter( hltLastFilterName );

    if ( mumHLTMatches.size() > 0 )
      mum_matched_lastfilter_name.append(hltLastFilterName+" ") ;

    if ( mupHLTMatches.size() > 0 )
      mup_matched_lastfilter_name.append(hltLastFilterName+" ") ;
  }

  mumtriglastfilter->push_back(mum_matched_lastfilter_name);
  muptriglastfilter->push_back(mup_matched_lastfilter_name);
}

void
LambdaB::saveTruthMatch(const edm::Event& iEvent){
  double deltaEtaPhi;

  for (vector<int>::size_type i = 0; i < lbmass->size(); i++) {

    //-----------------------
    // truth match with mu-  
    //-----------------------                                                                                                                        
    deltaEtaPhi = calEtaPhiDistance(genmumpx, genmumpy, genmumpz,
                                    mumpx->at(i), mumpy->at(i), mumpz->at(i));
    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemum->push_back(true);
    } else {
      istruemum->push_back(false);
    }

    //-----------------------
    // truth match with mu+  
    //-----------------------                                                                                                                       
    deltaEtaPhi = calEtaPhiDistance(genmuppx, genmuppy, genmuppz,
                                    muppx->at(i), muppy->at(i), muppz->at(i));

    if (deltaEtaPhi < TruthMatchMuonMaxR_) {
      istruemup->push_back(true);
    }
    else {
      istruemup->push_back(false);
    }

    //---------------------------------
    // truth match with proton track   
    //---------------------------------                                                                                                                       
    deltaEtaPhi = calEtaPhiDistance(genprpx, genprpy, genprpz,
                                    prpx->at(i), prpy->at(i), prpz->at(i));
    if (deltaEtaPhi < TruthMatchProtonMaxR_){
      istruepr->push_back(true);
    } else {
      istruepr->push_back(false);
    }

    //--------------------------------
    // truth match with pion track    
    //--------------------------------                                                                                                                      
    deltaEtaPhi = calEtaPhiDistance(genpipx, genpipy, genpipz,
                                    pipx->at(i), pipy->at(i), pipz->at(i));
    if (deltaEtaPhi < TruthMatchPionMaxR_){
      istruepi->push_back(true);
    } else {
      istruepi->push_back(false);
    }

    //---------------------------------------
    // truth match with /\b or /\b bar 
    //---------------------------------------                                                                                                
    if ( istruemum->back() && istruemup->back()
         && istruepr->back() && istruepi->back()) {
      istruelb->push_back(true);
    } else {
      istruelb->push_back(false);
    }
  }
}



//define this as a plug-in
DEFINE_FWK_MODULE(LambdaB);
