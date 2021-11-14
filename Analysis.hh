#ifndef Analysis_h
#define Analysis_h 1

#include <vector>
#include "math.h"
#include <cmath>

#include "globals.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4TrajectoryContainer.hh"
#include "G4SteppingManager.hh"
#include "G4SteppingVerbose.hh"
#include "G4SDManager.hh"
#include "G4THitsCollection.hh"
#include "G4HCofThisEvent.hh"
#include "G4Timer.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"

#include "Analysis_Messenger.hh"
#include "ROOTRecorder.hh"
#include "TrackerIonHit.hh"
#include "TrackerGammaHit.hh"
#include "G4ios.hh"

#include "DetectorConstruction.hh"
#include "Incoming_Beam.hh"
#include "Outgoing_Beam.hh"

class ROOTRecorder;

class Analysis
{


public:

  static Analysis* getInstance();

  //Constructor (takes pointers to the following objects)
  Analysis(DetectorConstruction*,Incoming_Beam*,Outgoing_Beam*);
  //destructor
  ~Analysis();

  // ! Recorder
  ROOTRecorder* Recorder;
  // ! Detector Construction object
  DetectorConstruction* detector;
  // ! Incomming beam
  Incoming_Beam* BeamIn;
  // ! Outgoing beam
  Outgoing_Beam* BeamOut;

  // ! Event number
  G4int			fEvtNum;
  //!  Number of Events processed in the run
  G4int			Total_Evts;
  // ! Run number
  G4int			fRunNum;	  
  
  /**
	* @brief Function called when G4Run starts
	* @param aRun : G4Run
	*/
  void BeginOfRun(const G4Run *aRun);
  
  /**
	* @brief Function called when G4Run stops
	* @param aRun : G4Run
	*/
  void EndOfRun(const G4Run *aRun);
  
  /**
	* @brief Function called when G4Event starts
	* @param aRun : G4Event
	*/
  void BeginOfEvent(const G4Event *anEvent);
  
  /**
	* @brief Function called when G4Event stops
	* @param anEvent : G4Event
	*/
  void EndOfEvent(const G4Event *anEvent);
     
  //! Initialization of the Analysis
  void Init();
  
  //! Initialization of variables at each event
  void InitEventData();

  //! Treatement of the hits
  void Treat(const G4Event *anEvent);
  
  //! Treatement of Ions 
  void TreatIons(TrackerIonHitsCollection* ionCollection); 

  //! Treatement of SeGa Gamma rays 
  void TreatSeGA(TrackerGammaHitsCollection* gammaCollection); 

  //! Treatement of Gretina Gamma rays 
  void TreatGretina(GretinaHitDetectorCollection* GretinaCollection); 
    
  //! Sorr Gamma-rays Energy to SeGA Segement
  void AddESegment(G4int r, G4int d, G4int s, G4int q, G4double E);
  
  //! Sort Gamma-rays Energy to Gretina Cristal
  void AddGretinaGeEnergy( G4double energy, G4int detnum, G4int clunum );

  //! Report the analysis parameters
  void Report(void);
  
  //! Set Ions axis method for Doppler Correction
  inline void Set_DC_IonAxis(G4int i){DC_IonAxis= i;};
  //! Set Ion Position method for Doppler Correction
  inline void Set_DC_IonPos(G4int i){DC_IonPos = i;};
  //! Set Gamma Position method for Doppler Correction
  inline void Set_DC_GammaPos(G4int i){DC_Gamma_Pos = i;};
  //! Set beta method for Doppler Correction
  inline void Set_DC_Beta(G4int i){DC_Beta = i;};
  //! Set FWHM simulation flag
  inline void Set_FWHM(G4bool f){FWHM_FLAG = f;};
  //! Set Sega Tracking flag
  inline void Set_SeGATracking(G4bool f){SeGA_TRACK_FLAG = f;};
  //! Set Gretina Tracking flag
  inline void Set_GretinaTracking(G4bool f){GRETINA_TRACK_FLAG = f;};
  //! Set Gretina Gamma Position method for Doppler Correction
  inline void Set_GretinaPosSigma(G4double s){DC_Gretina_Pos_Sigma = s;};
  //new in v4.2
  //! Set Max Distance for AddbackSphere
  inline void Set_GretinaAddBackSphereDist(G4double s){ MaxDist = s;};  
  //! Set S800 ACCEPTANCE flag
  void Set_S800_Accept(G4bool f){S800_ACCEPTANCE_FLAG = f;};  
  //! Set S800 Angular Acceptance lower limit
  void SetAcceptanceThetaLow(G4double x){AccThL=x;};
  //! Set S800 Angular Acceptance higher limit
  void SetAcceptanceThetaHigh(G4double x){AccThH=x;};
  //! Set S800 momnetum Acceptance lower limit
  void SetAcceptanceDLow(G4double x){AccDL=x;};
  //! Set S800 momentum Acceptance higher limit
  void SetAcceptanceDHigh(G4double x){AccDH=x;};

  //AR New in v4.3
  //! Set Disabled Cryst in Gretina
  void Set_GretDisCrys(G4int i){GretDisCrys = i;};
  //! Set Rob's Cascade Analysis
  void Set_CascadeAnalysisOn(G4bool f){CASCADE_ANALYSIS_FLAG = f;};
  //! Rob's Cascade Analysis ELow
  void Set_CascadeAnalysisElow(G4double x){Casc_Ana_ELow = x;};
  //! Rob's Cascade Analysis EHigh
  void Set_CascadeAnalysisEHigh(G4double x){Casc_Ana_EHigh = x;};
  //! Rob's Cascade Analysis ETrue
  void Set_CascadeAnalysisETrue(G4double x){Casc_Ana_ETrue = x;};
  
  //TM inlcuded in v4.2 for sphere addback. All these are just standard headers.
  void AddbackMode1(void); //CRL333 
  void AddbackMode3Simple(void); //CRL333
  void AddbackMode3Simple(GretinaHitDetectorCollection* ghdc); //CRL333
  float GetDistance(int i, int j); //CRL333 
  
  G4double GetMaxDist(){return MaxDist;};

  G4double GetDepth(G4ThreeVector pos, G4double tarz); //AR New in v4.3 -> Function to derive Depth of the interaction in GRETINA
  
  inline int GetMaxDecay(){return MaxDecay;};

  int GetLargestRemaining(vector< vector <int> > tested);//CRL333 This is to find largest hit not yet tested

  //AR New in v4.3 -> To include Rob's Cascade Analysis
  float GetDistance(G4ThreeVector p1, G4ThreeVector p2);
  int GetLargestRemaining(vector< vector <int> > tested, GretinaHitDetectorCollection* ghdc);
  double GetCosTheta(G4ThreeVector gpos);
  G4double CalcSpan(vector<int> abPoint, GretinaHitDetectorCollection* ghdc);
  double CosComptAngle(double enDep, double enRem);
  double CasZangle(G4ThreeVector posCompt, G4ThreeVector posNext, const G4ThreeVector& decay);
  
private:
  static Analysis* fManager;
  
  G4bool   SeGA_TRACK_FLAG;  


  static const G4int MaxDecay = 10;

  static const G4int MaxDetG = 4;
  static const G4int MaxClusG = 32;
  static const G4int MaxDet = 8;
  static const G4int MaxRing = 2;
  static const G4int MaxTrack = 200;
  G4double ESeg[MaxRing][MaxDet][8][4];

  // SeGA energy resoltution flag
  G4bool FWHM_FLAG;
  //! S800 Acceptance flag
  G4bool S800_ACCEPTANCE_FLAG;
  //! S800 Ions Acceptance flag depending on plunger configuration
  G4int RECO_FLAG;
  
  //! Doppler Correction Flags for ions axis
  G4int DC_IonAxis;
  //! Doppler Correction Flags for ions position
  G4int DC_IonPos;
  //! Doppler Correction Flags for gamma position
  G4int DC_Gamma_Pos;
  //! Doppler Correction Flags for beta of the ion
  G4int DC_Beta;

  //AR New in v4.3
  // Gretina Disable Crystal Number
  G4int GretDisCrys;
  // Set Rob's Cascade Analysis
  G4bool CASCADE_ANALYSIS_FLAG;
  // Rob's Cascade Analysis ELow
  G4double Casc_Ana_ELow;
  // Rob's Cascade Analysis EHigh
  G4double Casc_Ana_EHigh;
  // Rob's Cascade Analysis ETrue
  G4double Casc_Ana_ETrue;
  
public:  
  //! Missing Tracking Events 
  G4int MissCtr;

  G4int DetM;	 
  // Sega Variables
  G4int RingNr[MaxDet];
  G4int DetNr[MaxDet];
  G4int ExtDetNr[MaxDet];
  G4double DetE[MaxDet];
  G4double DetEDC[MaxDet];
  G4double DetEDC511[MaxDet];
 
  G4int MaxSegNr[MaxDet];
  G4int MaxS[MaxDet];
  G4int MaxQ[MaxDet];
  G4int MaxSegE;
  G4double SegTheta[MaxDet];
  G4double CosThetaCM[MaxDet];
  G4double CosThetaCMTrue[MaxDet];
  G4double SegPhi[MaxDet];

  G4int SegM;
  G4int SegMul[MaxDet*32];
  G4int SegNr[MaxDet*32];
  G4int SegQ[MaxDet*32];
  G4int SegS[MaxDet*32];
  G4int SegDetPtr[MaxDet*32];
  G4double SegE[MaxDet*32];
  
  //TM added in v4.2 for sphere addback
  G4double addbackbeta;
  G4double addbackcos;//CRL333
  /////
  
  //! Segment Position defined as middle of the volume
  G4ThreeVector Seg_Pos;
  //! Gamma Positon 
  G4ThreeVector Gamma_Pos;
  //! Gamma True Positon 
  G4ThreeVector Gamma_TruePos;
  //! Ion Positon 
  G4ThreeVector Ion_Pos;
  //! Ion Vector
  G4ThreeVector Ion_V;
  //! Gamma vector
  G4ThreeVector Gamma_V;
  //! Gamma true vector direction
  G4ThreeVector Gamma_VTrue;

  //! Target center  position
  G4ThreeVector Target_Pos;
  //! Degrader center  position
  G4ThreeVector Degrader_Pos;
  //! Stopper center  position
  G4ThreeVector Stopper_Pos;
  //! Plunger Type
  G4int PlungerType;



  G4double dTheta[MaxDet];
  G4double dPhi[MaxDet];
  G4double dR[MaxDet];
  G4double Theta[MaxDet];
  G4double Phi[MaxDet];
  G4double R[MaxDet];
  G4double dX[MaxDet];
  G4double dY[MaxDet];
  G4double dZ[MaxDet];


  G4double beta_r;
  G4double gamma_r;
	
  G4int    track_N;
  G4int    track_D[MaxTrack];
  G4int    track_R[MaxTrack];
  G4int    track_S[MaxTrack];
  G4int    track_Q[MaxTrack];
  G4double track_x[MaxTrack];
  G4double track_y[MaxTrack];
  G4double track_z[MaxTrack];
  G4double track_E[MaxTrack];

  // ROOT Tree stuff
  struct ion_event_t{
      double      beta;
      double      theta;
      double      phi;
      double      time;
      double      labtime;
      double      globaltime;
      double      KE;
      double      x;
      double      y;
      double      z;
  };

  struct ion_decay_t{
	 int         LevelID[MaxDecay];
	 double      LevelE[MaxDecay];
	 double      beta[MaxDecay];
	 double      theta[MaxDecay];
	 double      phi[MaxDecay];
	 double      time[MaxDecay];
	 double      labtime[MaxDecay];
	 double      globaltime[MaxDecay];
	 double      KE[MaxDecay];
	 double      x[MaxDecay];
	 double      y[MaxDecay];
	 double      z[MaxDecay];
  };

  ion_event_t ion_ev[MAX_FLAGS];

  //! Decay multiplicity
  G4int DecayM;	 

  //! Decay storage
  ion_decay_t ion_dec;

  struct deposits_t{
    double TarDep;
    double TarLen;
    double DegDep;
    double DegLen;
    double StoDep;
    double StoLen;
    int    ReactionFlag;
    double ReactionDTheta;
    double ReactionDPhi;
  };
  deposits_t deposits;
  

  G4double     AccThL,AccDL;
  G4double     AccThH,AccDH;
  G4double     sig_a,sig_b,sig_y,sig_d;
  
  G4double d_ta;
  G4double y_ta;
  G4double theta_ta;
  G4double phi_ta;
  G4double a_ta;
  G4double b_ta;
  G4double beta_ta;
  

  // Gretina Variable

  //! Total Cristal multiplicity
  G4int DetM_G;	 
  //! Cristal Energies
  G4double ECr_G[MaxClusG][MaxDetG];
  //! Cristal Energies
  G4double EBackCr_G[MaxClusG][MaxDetG];  
  //! Array of internal cristal numbers sorted by multiplicity
  G4int DetNr_G[MaxDetG*MaxClusG];
  //! Array of cluster numbers sorted by multiplicity
  G4int DetClNr_G[MaxDetG*MaxClusG];
  //! Array of cristal energies in laboratory frame sorted by multiplicity
  G4double DetE_G[MaxDetG*MaxClusG];
  G4double DetE_Background_G[MaxDetG*MaxClusG];
  //! Array of cristal of doppler Corrected energies sorted by multiplicity
  G4double DetEDC_G[MaxDetG*MaxClusG];
  G4double DetEDC_Background_G[MaxDetG*MaxClusG];
  //! Array Cos Angle between gamma interaction point and 
  G4double CosSeg_G[MaxDetG*MaxClusG];
  //! Array radial distance of gamma interaction point
  G4double R_G[MaxDetG*MaxClusG];
  //! Array theta angle of gamma interaction point
  G4double Theta_G[MaxDetG*MaxClusG];
  //! Array phi angle of gamma interaction point
  G4double Phi_G[MaxDetG*MaxClusG];
  //! Array radial distance deviation from true gamma interaction point
  G4double dR_G[MaxDetG*MaxClusG];
  //! Array theta angle deviation from true gamma interaction point
  G4double dTheta_G[MaxDetG*MaxClusG];
  //! Array phi angle deviation from true gamma interaction point
  G4double dPhi_G[MaxDetG*MaxClusG];
  //! Array X coordinate of  gamma interaction point
  G4double X_G[MaxDetG*MaxClusG];
  //! Array Y coordinate of gamma interaction point
  G4double Y_G[MaxDetG*MaxClusG];
  //! Array Z coordinate of gamma interaction point
  G4double Z_G[MaxDetG*MaxClusG];
  //! Total energy in the Gretina Array
  G4double TotalE_G;
  //! Array depth of gamma interaction point
  G4double Depth_G[MaxDetG*MaxClusG];
  
  //! AddbackEnergy added in v4.2 for sphere addback
  G4double AddbackEdc;//Energy for the addback Mode, to make more simply add in AddbackEdcModex and then put that into the form.
  G4double AddbackEdcMode1;
  vector<G4double> AddbackEdcModeSphere;
  vector<G4double> AddbackMult;
  vector<G4double> AddbackAngleModeSphere;
  //! AR New for v4.3
  vector<G4double> AddbackEModeSphere; 
  vector<G4ThreeVector> AddbackPos;
  
  G4bool   GRETINA_TRACK_FLAG;  
  G4int    trackG_N;
  G4int    trackG_C[MaxTrack];
  G4int    trackG_D[MaxTrack];
  G4int    trackG_S[MaxTrack];
  G4double trackG_x[MaxTrack];
  G4double trackG_y[MaxTrack];
  G4double trackG_z[MaxTrack];
  G4double trackG_E[MaxTrack];
  
  G4int highenmult;
  ////
  
  //! Doppler Correction Flags for beta of the ion
  G4double DC_Gretina_Pos_Sigma;
  G4double MaxDist; 

  //--- All the static run parameters to put in the inputs tree. TM added for v4.2.
  G4double     RatioTarget;
  G4double     RatioDegrader;
  G4double     TargetThickness;
  G4double     DegraderThickness;
  G4double     StopperThickness;
  G4double     TarScaleDensity;
  G4double     DegScaleDensity;
  G4double     StopScaleDensity;   
  G4double     Distance1;
  G4double     Distance2;
  
  G4int        NumOfStates;
  G4double     StateEnergies[20];
  G4double     GammaEnergies[20];
  G4double     GammaTau[20];
  G4double     GammaFrac[20];  
  
  
  G4int        NumOfBackgroundStates;
  G4double     Energy[20];
  G4double     BackgroundEnergies[20];
  G4double     BackgroundFrac[20];    
  G4double     BackgroundFracTotal;    
  
  // AR Added for v4.3 - Variables used to calculate the values that are filled in CascadeZ and Corrected0pE histograms.  
  vector<G4double> CascadeZ;
  vector<vector<G4double> > CorrectedE;
  vector<G4double> CascadeZ_AB;
  vector<vector<G4double> > CorrectedE_AB;
  vector<G4double> ABsphereRad;
  vector<G4double> ABcomptAngDiff;
  G4double Theta_GammaKnown;
  G4double Theta_GammaOfInterest;
  G4double TrueE_GammaKnown; 
  G4double RelGamma;
  G4double LabE_GammaKnown;
  G4double GammaKnown_Low;
  G4double GammaKnown_High;
  G4double LabE_GammaOfInterest;
  G4double Z_Tar_GammaKnown;
  G4double R_Tar_GammaKnown;
  G4double Z_Tar_GammaOfInterest;
  G4double R_Tar_GammaOfInterest;
  G4double Z_Origin_GammaKnown;
  G4double R_Origin_GammaKnown;
  G4double Z_Origin_GammaOfInterest;
  G4double R_Origin_GammaOfInterest;
  G4double Z_Origin_Tar;
  G4double Rho_Tar_GammaKnown;
  G4double Rho_Tar_GammaOfInterest;
  G4double Z_Decay_GammaKnown;
  G4double Z_Tar_Decay;
  G4double Z_Decay_GammaOfInterest;
  
public :

  //! Simulate S800 acceptance
  G4bool S800Accept(ion_event_t evt);


};



#endif
