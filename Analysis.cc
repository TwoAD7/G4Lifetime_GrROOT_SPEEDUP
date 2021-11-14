#include "Analysis.hh"

Analysis::Analysis(DetectorConstruction* DE, Incoming_Beam* BI, Outgoing_Beam* BO):detector(DE),BeamIn(BI),BeamOut(BO)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__ << G4endl;
#endif
  Recorder = new ROOTRecorder(this);

  // Set Default flags
  SeGA_TRACK_FLAG = FALSE; 
  // Gretina Tracking ... Has to be set on for now as no infos on segments is available
  GRETINA_TRACK_FLAG = TRUE; 
  FWHM_FLAG = FALSE;
  S800_ACCEPTANCE_FLAG = FALSE; 
  RECO_FLAG = TARGET_BACK_FLAG;
  // Set Default Doppler Correction Flags
  DC_IonAxis =  0; // z Axis
  DC_IonPos =  0; // Middle of target
  DC_Gamma_Pos = 0; // Largest Deposit segment
  DC_Beta = 0; // Tr

  //AR New in v4.3
  //By default, every detector works
  GretDisCrys = -1;
  //By default, no Cascade Analysis
  CASCADE_ANALYSIS_FLAG = FALSE;
  //Default values for Rob's Cascde Analsysis Parameters
  Casc_Ana_ELow = -1.;
  Casc_Ana_EHigh = -1.;
  Casc_Ana_ETrue = -1.;
    
  // Default Position resolution for Gretina 
  DC_Gretina_Pos_Sigma = 0*mm; // true position
  
  //Default distance for Addback Sphere
  MaxDist = 8.0*cm;  
  
  // Default Values for the S800 accpetance  
  AccThL = 0*rad ;
  AccThH = 1*rad;
  AccDL = -100;
  AccDH = 100; 

  sig_a=0.;
  sig_b=0.;
  sig_y=0.;
  sig_d=0.;

  Target_Pos.set(0,0,0);
  Degrader_Pos.set(0,0,0);
  Stopper_Pos.set(0,0,0);
 
  //
  MissCtr= 0;

}


Analysis::~Analysis()
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__ << G4endl;
#endif
  delete Recorder;
}

void Analysis::BeginOfRun(const G4Run *aRun)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__ << G4endl;
#endif

  Recorder->BeginOfRun();

  Init();

  fRunNum = aRun->GetRunID();
  Total_Evts = aRun->GetNumberOfEventToBeProcessed();
  
  
  
#ifdef DEBUG
  G4cout << "Run num. " << fRunNum << G4endl;
#endif
}


void Analysis::EndOfRun(const G4Run *aRun)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  <<  aRun->GetRunID() << G4endl;
#endif

  Recorder->EndOfRun();
  G4cout << "\n\t\033[1m\033[31m------------------------------------------------------------  \033[0m" << G4endl;
  G4cout << "\t\033[1m\033[31m End of Simulation Notes  : \033[0m" << G4endl;
  G4cout << "\t\033[1m\033[31m------------------------------------------------------------ \033[0m" << G4endl;
  G4cout << "\t \033[1m\033[31m" << aRun->GetNumberOfEventToBeProcessed() << "\033[0m events simulated " <<  G4endl;
  G4cout << "\t Simulation saved to the file : \033[1m\033[31m" << Recorder->outTreeName  << "\033[0m" <<  G4endl;
  if(MissCtr >0)  G4cout << "\t Strange (too long) / Missing  \033[1m\033[31m" << MissCtr << "\033[0m tracks " <<  G4endl;
  G4cout << "\t\033[1m\033[31m------------------------------------------------------------ \033[0m" << G4endl;
}  

void Analysis::BeginOfEvent(const G4Event *anEvent)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif

  InitEventData();
  fEvtNum = anEvent->GetEventID();

#ifdef DEBUG
  G4cout << "Event num. " << fEvtNum << G4endl;
#endif
}


void Analysis::EndOfEvent(const G4Event *anEvent)
{

#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif

  Treat(anEvent);
   
  Recorder->FillTree();
  Recorder->FillHistograms();
  
  //TM v4.2 added for Charlie's addback
  AddbackEdcModeSphere.clear();
  AddbackAngleModeSphere.clear();
  AddbackMult.clear();

  //AR New for v4.3 -> Rob's Cascade Analysis
  AddbackEModeSphere.clear();
  AddbackPos.clear();
  ABcomptAngDiff.clear();
  ABsphereRad.clear();
  AddbackEdcModeSphere.clear();
  AddbackAngleModeSphere.clear();
  AddbackMult.clear();
  CascadeZ.clear();
  CascadeZ_AB.clear();
  CorrectedE.clear();
  CorrectedE_AB.clear();
  ////////
  
  if((fEvtNum+1)%100 == 0||fEvtNum%(Total_Evts-1) == 0)
	 {
		G4cout <<"\rProcessing: Events: " <<
		  Total_Evts << ", Done: " << (fEvtNum+1) <<
		  ", (%): " << ((float) (fEvtNum+1))/((float) Total_Evts)*100.;
		G4cout.flush();
	 }
	
}

void Analysis::InitEventData()
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif

  // Initialize Ions variables
  DecayM=0;

  for(G4int f=0;f<MAX_FLAGS;f++) {
	 ion_ev[f].beta=-1;
	 ion_ev[f].theta=-1;
	 ion_ev[f].phi=-1;
	 ion_ev[f].time=-1;
	 ion_ev[f].labtime=-1;
	 ion_ev[f].globaltime=-1;
	 ion_ev[f].KE=-1;
	 ion_ev[f].x=-9999.*mm;
	 ion_ev[f].y=-9999.*mm;
	 ion_ev[f].z=-9999.*mm;
  }
  for (G4int d=0;d<MaxDecay;d++){
	 ion_dec.LevelID[d]= -1;
	 ion_dec.LevelE[d]= -1;
	 ion_dec.beta[d]= -1;
	 ion_dec.theta[d]=-1;
	 ion_dec.phi[d]=-1;
	 ion_dec.time[d]=-1;
	 ion_dec.labtime[d]=-1;
	 ion_dec.globaltime[d]=-1;
	 ion_dec.KE[d]=-1;
	 ion_dec.x[d]=-9999.*mm;
	 ion_dec.y[d]=-9999.*mm;
	 ion_dec.z[d]=-9999.*mm;
  }

  
  y_ta=-9999.*mm;
  d_ta=-1000;
  theta_ta = -1;
  phi_ta = -1;
  a_ta = -100;
  b_ta = -100;

  // Initialize Gamma variables
  DetM = 0;
  SegM = 0;
  for (G4int d=0;d<MaxDet;d++){
	 RingNr[d] = -1;
	 DetNr[d] = -1;
	 ExtDetNr[d] = -1;
	 DetE[d] = 0.;
	 DetEDC[d] = 0.;
	 SegTheta[d] = 0.;
	 SegPhi[d] = 0.;
	 SegMul[d] = 0;
	 CosThetaCM[d] = -10.;
	 CosThetaCMTrue[d] = -10;
	 dTheta[d] = 0.;
	 dPhi[d] = 0.;
	 R[d] = 0.;
	 Theta[d] = 0.;
	 Phi[d] = 0.;
	 dR[d] = 0.;
	 dX[d] = 0.;
	 dY[d] = 0.;
	 dZ[d] = 0.;
	 MaxSegNr[d] = -1;
	 MaxS[d] = -1;
	 MaxQ[d] = -1;
	 for(G4int s=0;s<8;s++){
		for(G4int q=0;q<4;q++){
		  ESeg[0][d][s][q] = 0.;
		  ESeg[1][d][s][q] = 0.;
		  SegE[d*32+(s+q*8)] = 0.;
		  SegNr[d*32+(s+q*8)] = -1;
		  SegS[d*32+(s+q*8)] = -1;
		  SegQ[d*32+(s+q*8)] = -1;
		  SegDetPtr[d*32+(s+q*8)]=-1;	
		}
	 }
  }

  if(SeGA_TRACK_FLAG == TRUE) {
	 track_N = 0;
	 for(int n=0;n<MaxTrack;n++){
		track_D[n]=0;
		track_R[n]=0;
		track_S[n]=0;
		track_Q[n]=0;
		track_x[n]=0;
		track_y[n]=0;
		track_z[n]=0;
		track_E[n]=0;
	 }
  }	


  // For Gretina
  DetM_G = 0;
  for (G4int d=0;d<MaxDetG;d++){
	 for(G4int c=0;c<MaxClusG;c++){
		ECr_G[c][d] = 0.;
		DetClNr_G[4*c+d] = -1;
		DetNr_G[4*c+d] = -1;
		DetE_G[4*c+d] = 0.;	
		CosSeg_G[4*c+d] = -10;
		DetEDC_G[4*c+d] = 0.;
		dTheta_G[4*c+d] = 0.;
		dPhi_G[4*c+d] = 0.;
		R_G[4*c+d] = 0.;
		Theta_G[4*c+d] = 0.;
		Phi_G[4*c+d] = 0.;
		Depth_G[4*c+d] = 0.;
		X_G[4*c+d] = 0.;
		Y_G[4*c+d] = 0.;
		Z_G[4*c+d] = 0.;

	 }
	 
  }
  TotalE_G=0.;

  if(GRETINA_TRACK_FLAG == TRUE) {
	 trackG_N = 0;
	 for(int n=0;n<MaxTrack;n++){
		trackG_C[n]=0;
		trackG_D[n]=0;
		trackG_S[n]=0;
		trackG_x[n]=0;
		trackG_y[n]=0;
		trackG_z[n]=0;
		trackG_E[n]=0;
	 }
  }	


  Gamma_Pos.set(0,0,0);
  Gamma_TruePos.set(0,0,0);
  Ion_Pos.set(0,0,0);
  Seg_Pos.set(0,0,0);
  Ion_V.set(0,0,0);
  Gamma_V.set(0,0,0);
  Gamma_VTrue.set(0,0,0);
 

}

void Analysis::Init()
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif
  PlungerType = detector->GetPlungerType();

  // Select the Reconstruction for acceptance cut in S800 depending of
  // the Plunger Configuration
  switch(PlungerType){
  case 0 : // Target Only
	 RECO_FLAG = TARGET_BACK_FLAG;
	 break;
  case 1 : // Standard Plunger
	 RECO_FLAG = DEGRADER_BACK_FLAG;
	 break;
  case 2 : // Diff Plunger
	 RECO_FLAG = STOPPER_BACK_FLAG;
	 break;
  }


  // Get the position of Target, Degrader and Stopper
  
  Target_Pos = detector->GetTargetPlacement()->GetTranslation();
  if(detector->GetPlungerType() > 0){
	 Degrader_Pos = detector->GetDegraderPlacement()->GetTranslation();
	 if(detector->GetPlungerType() == 2){
		Stopper_Pos = detector->GetStopperPlacement()->GetTranslation();
	 }
  }
  
  //TM added in v4.2 for static parameters, input histograms in root file
  RatioTarget = detector->GetRatioT();
  RatioDegrader = detector->GetRatioD();  
  TargetThickness = detector->GetTargetThick();
  DegraderThickness = detector->GetDegraderThick();
  StopperThickness = detector->GetStopperThick();
  TarScaleDensity = detector->GetTargetScaleDensity();
  DegScaleDensity = detector->GetDegraderScaleDensity();
  StopScaleDensity = detector->GetStopperScaleDensity();
  Distance1 = detector->GetDis1();
  Distance2 = detector->GetDis2();   		
  
  vector<Energy_Level*> EnergyInfo = BeamOut->GetStateInfo();
  NumOfStates = EnergyInfo.size();
  for( int i = 0; i < NumOfStates; i++ )
  {
      StateEnergies[i] = EnergyInfo[i]->GetEx()*1.0e3;
      GammaEnergies[i] = EnergyInfo[i]->GetEg()*1.0e3;
      GammaTau[i] = EnergyInfo[i]->GetTau()*1.0e3;
      GammaFrac[i] = EnergyInfo[i]->GetFrac();
  }  
  
  //
  MissCtr = 0;

}

void Analysis::Treat(const G4Event *Event)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif

  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4int ionCollectionID=SDman->GetCollectionID("ionCollection");
  
  G4HCofThisEvent * HCE = Event->GetHCofThisEvent();

  // Treat Ions First
  TreatIons((TrackerIonHitsCollection*)(HCE->GetHC(ionCollectionID)));

  // Then Gamma 
  if(detector->IsDefinedSeGA()){
	 G4int gammaCollectionID=SDman->GetCollectionID("gammaCollection");
	 TreatSeGA((TrackerGammaHitsCollection*)(HCE->GetHC(gammaCollectionID)));
  }
  if(detector->IsDefinedGretina()) {
	 G4int GretinaCollectionID=SDman->GetCollectionID("GretinaTracker");
	 TreatGretina((GretinaHitDetectorCollection*)(HCE->GetHC(GretinaCollectionID)));
	 if(!CASCADE_ANALYSIS_FLAG){
	   AddbackMode3Simple();
	 }else{
	   AddbackMode3Simple((GretinaHitDetectorCollection*)HCE->GetHC(GretinaCollectionID));// This is a similar approximation as the mode 1 described above. The primary difference is that this does not create a connected region but just adds up all points within some distance
	 }
  }

}

void Analysis::TreatIons(TrackerIonHitsCollection* ionCollection)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif

  G4int  NIons=ionCollection->entries();
  G4int f;
  
  for(G4int i=0;i<NIons;i++){
	 f=(*ionCollection)[i]->GetFlag();
	 ion_ev[f].beta=(*ionCollection)[i]->GetBeta();
	 ion_ev[f].theta=(*ionCollection)[i]->GetTheta()/rad;
	 ion_ev[f].phi=(*ionCollection)[i]->GetPhi()/rad;
	 ion_ev[f].time=(*ionCollection)[i]->GetTime();
	 ion_ev[f].labtime=(*ionCollection)[i]->GetLabTime();
	 ion_ev[f].globaltime=(*ionCollection)[i]->GetGlobTime();
	 ion_ev[f].KE=(*ionCollection)[i]->GetKE()/MeV;
	 ion_ev[f].x=(*ionCollection)[i]->GetPos().getX()/mm;
	 ion_ev[f].y=(*ionCollection)[i]->GetPos().getY()/mm;
	 ion_ev[f].z=(*ionCollection)[i]->GetPos().getZ()/mm;


	 if(f==DECAY_FLAG){
		int b = (*ionCollection)[i]->GetParticleID().find_first_of('[',  1);
		int e = (*ionCollection)[i]->GetParticleID().find_first_of(']', b);
		ion_dec.LevelID[DecayM]= (*ionCollection)[i]->GetParentID();
		ion_dec.LevelE[DecayM]=  atof(((*ionCollection)[i]->GetParticleID().substr(b+1,e-b-1)).c_str());
		ion_dec.beta[DecayM]= (*ionCollection)[i]->GetBeta();
		ion_dec.theta[DecayM]=(*ionCollection)[i]->GetTheta()/rad;
		ion_dec.phi[DecayM]=(*ionCollection)[i]->GetPhi()/rad;
		ion_dec.time[DecayM]=(*ionCollection)[i]->GetTime();
		ion_dec.labtime[DecayM]=(*ionCollection)[i]->GetLabTime();
		ion_dec.globaltime[DecayM]=(*ionCollection)[i]->GetGlobTime();
		ion_dec.KE[DecayM]=(*ionCollection)[i]->GetKE()/MeV;
		ion_dec.x[DecayM]=(*ionCollection)[i]->GetPos().getX()/mm;
		ion_dec.y[DecayM]=(*ionCollection)[i]->GetPos().getY()/mm;
		ion_dec.z[DecayM]=(*ionCollection)[i]->GetPos().getZ()/mm;
		DecayM++;
	 }
  }
  
  G4bool tlog,dlog,slog;
  G4ThreeVector u,v,ez,au,av;
  u.setX(0.);u.setY(0.);u.setZ(1.);
  v.setX(0.);v.setY(0.);v.setZ(1.);
  ez.setX(0.);ez.setY(0.);ez.setZ(1.);
  deposits.TarLen=0.;
  deposits.TarDep=0.;
  deposits.DegLen=0;
  deposits.DegDep=0.;
  deposits.StoLen=0;
  deposits.StoDep=0.;
  tlog=false;dlog=false;slog=false;
		
  for(G4int i=1;i<NIons;i++)
	 {
		if((*ionCollection)[i]->GetFlag()==TARGET_FACE_FLAG)
		  {tlog=true;
			 deposits.TarLen=-(*ionCollection)[i-1]->GetLength()/mm;}
		if((*ionCollection)[i-1]->GetFlag()==TARGET_BACK_FLAG)
		  tlog=false;
		if(tlog)
		  {
			 deposits.TarDep+=(*ionCollection)[i-1]->GetEdep()/MeV;
			 deposits.TarLen+=(*ionCollection)[i-1]->GetLength()/mm;
		  }
			 
		if((*ionCollection)[i]->GetFlag()==DEGRADER_FACE_FLAG)
		  {dlog=true;
			 deposits.DegLen=-(*ionCollection)[i-1]->GetLength()/mm;}
		if((*ionCollection)[i-1]->GetFlag()==DEGRADER_BACK_FLAG)
		  dlog=false;
		if(dlog)
		  {
			 deposits.DegDep+=(*ionCollection)[i-1]->GetEdep()/MeV;
			 deposits.DegLen+=(*ionCollection)[i-1]->GetLength()/mm;
		  }
			 
		if((*ionCollection)[i]->GetFlag()==STOPPER_FACE_FLAG)
		  {slog=true;
			 deposits.StoLen=-(*ionCollection)[i-1]->GetLength()/mm;}
		if((*ionCollection)[i-1]->GetFlag()==STOPPER_BACK_FLAG)
		  slog=false;
		if(slog)
		  {
			 deposits.StoDep+=(*ionCollection)[i-1]->GetEdep()/MeV;
			 deposits.StoLen+=(*ionCollection)[i-1]->GetLength()/mm;
		  }
		if((*ionCollection)[i]->GetFlag()==REACTION_FLAG)
		  {
			 u.setTheta((*ionCollection)[i]->GetTheta());
			 u.setPhi((*ionCollection)[i]->GetPhi());
			 v.setTheta((*ionCollection)[i+1]->GetTheta());
			 v.setPhi((*ionCollection)[i+1]->GetPhi());
		  }
	 }
		
  //	setTargetEdepAverage(deposits.TarDep);
		
  deposits.ReactionFlag=BeamOut->GetReactionFlag();
  deposits.ReactionDTheta=acos(u.dot(v));
  if(ion_ev[REACTION_FLAG].beta>0)
	 {
		au=u.cross(ez);au.setMag(1.);
		av=u.cross(v);av.setMag(1.);
		if(au.dot(av)>=0) // to be ckeced
		  deposits.ReactionDPhi=acos(au.dot(av));
		else
		  deposits.ReactionDPhi=twopi-acos(au.dot(av));
	 }
  else deposits.ReactionDPhi=0.;
  

  

  if(S800Accept(ion_ev[RECO_FLAG])){
    d_ta = (ion_ev[RECO_FLAG].KE/BeamOut->GetURsetKE()-1.)*100.;
    y_ta = ion_ev[RECO_FLAG].y;
    theta_ta = ion_ev[RECO_FLAG].theta;
    phi_ta = ion_ev[RECO_FLAG].phi;
    a_ta=atan(tan(theta_ta)*cos(phi_ta))/mrad;
    b_ta=atan(tan(theta_ta)*sin(phi_ta))/mrad;
    theta_ta = theta_ta/mrad;
    phi_ta  = phi_ta/mrad;
    beta_ta = ion_ev[RECO_FLAG].beta;
	 
  }

}

void Analysis::TreatSeGA(TrackerGammaHitsCollection* gammaCollection)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif 
  G4int NbHits=gammaCollection->entries();
	 
  if(S800Accept(ion_ev[RECO_FLAG])){
	 if(NbHits>0) {
		if(NbHits > MaxTrack && SeGA_TRACK_FLAG == TRUE){
		  G4cerr << "\n Event # " << fEvtNum << "  with too many tracks (" << NbHits << ") -- Skipping Event" <<G4endl;
		  MissCtr++;
		  return;
		}
		
		G4int r,d,s,q;
		G4int HitPattern[MaxRing][MaxDet];

		memset(&HitPattern,0,sizeof(HitPattern));

		for (G4int i=0;i<NbHits;i++){
		  r=(*gammaCollection)[i]->GetRingID();
		  d=(*gammaCollection)[i]->GetDetNumb();
		  s=(*gammaCollection)[i]->GetSliceNbf();
		  q=(*gammaCollection)[i]->GetQuartNbf();
		  
		  if(strcmp((*gammaCollection)[i]->GetParticleID(),"gamma") == 0  || strcmp((*gammaCollection)[i]->GetParticleID(),"e-") == 0 )
				HitPattern[r-1][d] = 1;
			 
		
		  // Keep track of the individual interractions of gamma
		  if(SeGA_TRACK_FLAG == TRUE){
			 if(strcmp((*gammaCollection)[i]->GetParticleID(),"gamma") == 0 
				 || strcmp((*gammaCollection)[i]->GetParticleID(),"e-") == 0 ){
				track_D[track_N]=d+1;
				track_R[track_N]=r;
				track_S[track_N]=s;
				track_Q[track_N]=q;
				track_x[track_N]=(*gammaCollection)[i]->GetPos().getX()/mm;
				track_y[track_N]=(*gammaCollection)[i]->GetPos().getY()/mm;
				track_z[track_N]=(*gammaCollection)[i]->GetPos().getZ()/mm;
				track_E[track_N]=(*gammaCollection)[i]->GetEdep()/keV;
				track_N++;
			 }
			 else {
#ifdef DEBUG
				G4cout <<  " \n ========> " <<  G4endl; 
				G4cout <<  " \t Particle :" << (*gammaCollection)[i]->GetParticleID() <<  G4endl; 
				G4cout <<  " \t Energy   : " <<(*gammaCollection)[i]->GetEdep()/keV << "keV" <<  G4endl; 
				G4cout <<  " \t Pos X: " << (*gammaCollection)[i]->GetPos().getX()/mm  << " Y: " << (*gammaCollection)[i]->GetPos().getY()/mm  << " Z: " << (*gammaCollection)[i]->GetPos().getZ()/mm   <<  G4endl; 
				G4cout <<  " \t Ring : " << r  << " Detector: " <<  d <<  G4endl;
#endif		 	
			 }
		  }
		  // /!\ Counting Ring [1:2] Detectors [0,7]
		  if(strcmp((*gammaCollection)[i]->GetParticleID(),"gamma") == 0  || 
			  strcmp((*gammaCollection)[i]->GetParticleID(),"e-") == 0 )	
			 AddESegment(r-1,d,s,q,(*gammaCollection)[i]->GetEdep()/keV);
		}
		// Sort by Mutiplicty 
		for (G4int r=0;r<MaxRing;r++){
		  for (G4int d=0;d<MaxDet;d++){  
			 if(HitPattern[r][d] == 1){
				RingNr[DetM] = r+1;
				DetNr[DetM] = d+1;
				ExtDetNr[DetM] = detector->GetSeGA()->GetRing(RingNr[DetM])->GetDetector(DetNr[DetM]-1)->getID(); 
				MaxSegE = 0.;
				for(s=0;s<8;s++){
				  for(q=0;q<4;q++){
					 if(ESeg[r][d][s][q]>0){
						DetE[DetM]+=ESeg[r][d][s][q];
						SegNr[SegM] = detector->GetSeGA()->GetRing(RingNr[DetM])->GetDetector(DetNr[DetM]-1)->getSegID(s,q); 
						SegS[SegM] = s;
						SegQ[SegM] = q;
						SegE[SegM] =  ESeg[r][d][s][q];
						SegDetPtr[SegM] = d;
						// Max energy deposit 
						if(ESeg[r][d][s][q] > MaxSegE){
						  MaxSegNr[DetM] = SegNr[SegM];
						  MaxS[DetM] = s;
						  MaxQ[DetM] = q;
						  MaxSegE = ESeg[r][d][s][q];
						}
						SegM++;
						// Count Multiplicity per det
						SegMul[DetM]++;
					 }
				  }
				}
				
				// Apply resolution
				if(FWHM_FLAG){
				  DetE[DetM] = detector->GetSeGA()->GetRing(RingNr[DetM])->GetDetector(DetNr[DetM]-1)->FWHM_response(DetE[DetM]);
				}
 
				// Position of the Heavy Ion
				switch(DC_IonPos){
				case 0 : // Mid-target
				  Ion_Pos = Target_Pos;
				  break;
				case 1 : // True decay Position
				  Ion_Pos.setX(ion_dec.x[DecayM-1]*mm);
				  Ion_Pos.setY(ion_dec.y[DecayM-1]*mm);
				  Ion_Pos.setZ(ion_dec.z[DecayM-1]*mm);
				  break;
				case 2 : // User Defined Position
				  Ion_Pos=BeamOut->GetPosDopp();
				  break;
				case 3 : // Mid-Degrader
				  Ion_Pos = Degrader_Pos;
				  break;
				case 4 : // Mid-Stopper
				  Ion_Pos = Stopper_Pos;
				  break;						
				}
			 
				// Direction of the Heavy Ion
				switch(DC_IonAxis){
				case 0 : // z axis
				  Ion_V.setZ(1.);
				  Ion_V.setY(0.);
				  Ion_V.setX(0.);
				  break;
				case 1 : // Incomming Beam
				  Ion_V.setZ(1./sqrt(1.+tan(BeamIn->getAta0())*tan(BeamIn->getAta0())+tan(BeamIn->getBta0())*tan(BeamIn->getBta0())));
				  Ion_V.setY(tan(BeamIn->getBta0())*1./sqrt(1.+tan(BeamIn->getAta0())*tan(BeamIn->getAta0())+tan(BeamIn->getBta0())*tan(BeamIn->getBta0())));
				  Ion_V.setX(tan(BeamIn->getAta0())*1./sqrt(1.+tan(BeamIn->getAta0())*tan(BeamIn->getAta0())+tan(BeamIn->getBta0())*tan(BeamIn->getBta0())));
				  break;
				case 2 : // Decay
				  Ion_V.setZ(1.);
				  Ion_V.setTheta(ion_dec.theta[DecayM-1]);
				  Ion_V.setPhi(ion_dec.phi[DecayM-1]);
				  break;				
				case 3 : // Target back
				  Ion_V.setZ(1.);
				  Ion_V.setTheta(ion_ev[TARGET_BACK_FLAG].theta);
				  Ion_V.setPhi(ion_ev[TARGET_BACK_FLAG].phi);
				  break;	
				case 4 : // Degrader Back
				  Ion_V.setZ(1.);
				  Ion_V.setTheta(ion_ev[DEGRADER_BACK_FLAG].theta);
				  Ion_V.setPhi(ion_ev[DEGRADER_BACK_FLAG].phi);
				  break;				
				case 5 : // Stopper Back
				  Ion_V.setZ(1.);
				  Ion_V.setTheta(ion_ev[STOPPER_BACK_FLAG].theta);
				  Ion_V.setPhi(ion_ev[STOPPER_BACK_FLAG].phi);
				  break;	
				}
				Seg_Pos = detector->GetSeGA()->GetRing(RingNr[DetM])->GetDetector(DetNr[DetM]-1)->GetSegPos(MaxS[DetM],MaxQ[DetM]);		 
				
				// Gamma First Hit Position
				// True First Hit
				bool Trackfound = false;
				if(SeGA_TRACK_FLAG == TRUE){
				  for(G4int n=track_N-1;n>=0;n--){
#ifdef DEBUG
						G4cout << n << " --  track_D[n] " << track_D[n] << " DetNr " << DetNr[DetM] << "-  Track R " << track_R[n] << " - RingNr " <<  RingNr[DetM]<<  " X: " <<track_x[n]*mm << " Y : " << track_y[n]*mm << " Z : " << track_z[n]*mm <<  G4endl;
#endif			
					 if(track_D[n]==DetNr[DetM]&&track_R[n]==RingNr[DetM]){
						Gamma_TruePos.setX(track_x[n]*mm);
						Gamma_TruePos.setY(track_y[n]*mm);
						Gamma_TruePos.setZ(track_z[n]*mm);
						Trackfound = true;
					 }
				  }
				  
				  
				  if(!Trackfound){
					 MissCtr++;
					 G4cerr << "\n (1) Missing Track at event # " << fEvtNum<< " -- Skipping Event" <<endl;
					 continue;
				  }
				}
				switch(DC_Gamma_Pos){
				case 0 : // Largest Energy deposit segment Hit
				  Gamma_Pos= Seg_Pos;
				  break;
				case 1 : // Tracking : First Hit 
				  Gamma_Pos = Gamma_TruePos;
				  break;		
				}		
#ifdef DEBUG
 				G4cout << " Segment Position Segn " << Seg_Pos.getX()/cm  << "cm " << Seg_Pos.getY()/cm  << "cm " << Seg_Pos.getZ()/cm  << "cm - " << Seg_Pos.getR()/cm  << "cm " << Seg_Pos.getTheta()/deg  << "deg " << Seg_Pos.getPhi()/deg  << "deg "  << endl; 
 				G4cout << " Used Position          " << Gamma_Pos.getX()/cm << "cm " << Gamma_Pos.getY()/cm  << "cm " << Gamma_Pos.getZ()/cm << "cm - " << Gamma_Pos.getR()/cm  << "cm " << Gamma_Pos.getTheta()/deg  << "deg " << Gamma_Pos.getPhi()/deg  << "deg "  << endl; 
				if(SeGA_TRACK_FLAG == TRUE){
				  G4cout << " Real Position  (if tracking)        " << Gamma_TruePos.getX()/cm << "cm " << Gamma_TruePos.getY()/cm  << "cm " << Gamma_TruePos.getZ()/cm << "cm - " << Gamma_TruePos.getR()/cm  << "cm " << Gamma_TruePos.getTheta()/deg  << "deg " << Gamma_TruePos.getPhi()/deg  << "deg "  << endl; 
				}
#endif
		
				//direction of the gamma-ray
				Gamma_V=Gamma_Pos-Ion_Pos;

				//cos of the angle for doppler corrections
				Gamma_V.setMag(1.);
				Ion_V.setMag(1.);
				G4double cosSeg=Ion_V.dot(Gamma_V);

				//true direction of the gamma-ray
				if(SeGA_TRACK_FLAG == TRUE){
				  Gamma_VTrue = Gamma_TruePos-Ion_Pos;
				  Gamma_VTrue.setMag(1.);
				}
				
				// cos angle of the real position 			
				G4double	cosSegTrue = Ion_V.dot(Gamma_VTrue);
				
				// Velocity of the source for doppler corrections
				switch(DC_Beta){
				case 0 : // Use user defined Beta
				  beta_r=BeamOut->GetBetaDopp();
				  gamma_r=BeamOut->GetGammaDopp();
				  break;
				case 1 : // Use beta after target 
					beta_r=ion_ev[TARGET_BACK_FLAG].beta;
					gamma_r=1./sqrt(1-beta_r*beta_r);
					break;
				case 2 : // Use True beta at decay
				  beta_r=ion_dec.beta[DecayM-1];
				  gamma_r=1./sqrt(1-beta_r*beta_r);
				  break;
				case 3 : // Use beta after degrader
				  beta_r=ion_ev[DEGRADER_BACK_FLAG].beta;
				  gamma_r=1./sqrt(1-beta_r*beta_r);
				  break;
				case 4 : // Use beta after stopper
				  beta_r=ion_ev[STOPPER_BACK_FLAG].beta;
				  gamma_r=1./sqrt(1-beta_r*beta_r);
				  break;
				}

				
				 R[DetM] = (Gamma_Pos - Ion_Pos).getR()/cm;
				 Theta[DetM] = (Gamma_Pos - Ion_Pos).getTheta()/deg;
				 Phi[DetM] = (Gamma_Pos - Ion_Pos).getPhi()/deg;

				 if(SeGA_TRACK_FLAG == TRUE){
					dR[DetM] = (Gamma_Pos.getR()-Gamma_TruePos.getR())/mm;
					dTheta[DetM] = (Gamma_Pos.getTheta()-Gamma_TruePos.getTheta())/deg;
					dPhi[DetM] = (Gamma_Pos.getPhi()-Gamma_TruePos.getPhi())/deg;
					dX[DetM] = (Gamma_Pos.getX()-Gamma_TruePos.getX())/mm ;
					dY[DetM] = (Gamma_Pos.getY()-Gamma_TruePos.getY())/mm ;
 				  dZ[DetM] = (Gamma_Pos.getZ()-Gamma_TruePos.getZ())/mm ;
				 }
				 // Get angles information 
				 SegTheta[DetM] =  Seg_Pos.getTheta();
				 SegPhi[DetM] =  Seg_Pos.getPhi();

				 CosThetaCM[DetM] = (-beta_r*gamma_r + gamma_r * cosSeg)/(gamma_r - beta_r*gamma_r*cosSeg);
				 // True Angle in CM of the gamma
				 CosThetaCMTrue[DetM] = (-beta_r*gamma_r + gamma_r * cosSegTrue)/(gamma_r - beta_r*gamma_r*cosSegTrue);
				 
				 
				 DetEDC[DetM] = DetE[DetM] * gamma_r * (1-beta_r*cosSeg);
				 
				 // prompt 
				 //				DetEDC511[DetM] = detector->GetSeGA()->GetRing(RingNr[DetM])->GetDetector(DetNr[DetM]-1)->FWHM_response(511.) * gamma_r * (1-beta_r*cosSeg);
				DetM++;
			 }
		  }
		}
	 }
  }
}


void Analysis::TreatGretina( GretinaHitDetectorCollection* GretinaCollection)
{
  //TM added in v4.2 for sphere addback 

  G4float Highesten=0;//CRL333 This serves to keep track of the highest energy present in the collection to determine what point corresponds to the highest energy hit
  highenmult=0;
  //////////

#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif 

  G4int NbHits=GretinaCollection->entries();
	 
  if(S800Accept(ion_ev[RECO_FLAG]))
  {
    if(NbHits>0) {
		if(NbHits > MaxTrack && GRETINA_TRACK_FLAG == TRUE){
		  G4cerr << "\n Event # " << fEvtNum << "  with too many tracks (" << NbHits << ") -- Skipping Event" <<G4endl;
		  MissCtr++;
		  return;
		}
		
      G4int d,c,s;
      G4int HitPatternG[MaxClusG][MaxDetG];
		
      memset(&HitPatternG,0,sizeof(HitPatternG));
		
      for (G4int i=0;i<NbHits;i++){
	//AR New in v4.3 -> Condition on disabled crystal
	if((*GretinaCollection)[i]->GetDetNb()!=GretDisCrys){
	  c=(*GretinaCollection)[i]->GetCluNb();
	  d=(*GretinaCollection)[i]->GetDetNb()-c*4; 
	  s=(*GretinaCollection)[i]->GetSegNb();
	  //G4cout << "NCluster = " << c << " NDet = " << d+4*c << " NSeg = " << s << G4endl;
	  if(strcmp((*GretinaCollection)[i]->GetParticleID(),"gamma") == 0  || 
	     strcmp((*GretinaCollection)[i]->GetParticleID(),"e-") == 0 )
	    HitPatternG[c][d] = 1;	  
			 
		
#ifdef DEBUG
	  printf(" CluNb: %d, DetNb : %d (%d), SegNb %d, Edep : %f  %s %d\n",c,d,(*GretinaCollection)[i]->GetDetNb(),s,(*GretinaCollection)[i]->GetEdep()/keV,(*GretinaCollection)[i]->GetParticleID().c_str(),(*GretinaCollection)[i]->GetInterNb());			 
#endif
				
	  // Keep track of the individual interractions of gamma
	  if(GRETINA_TRACK_FLAG == TRUE){
	    if(strcmp((*GretinaCollection)[i]->GetParticleID(),"gamma") == 0  || 
	       strcmp((*GretinaCollection)[i]->GetParticleID(),"e-") == 0) {
	      trackG_C[trackG_N]=c;
	      trackG_D[trackG_N]=d;
	      trackG_S[trackG_N]=s;
	      trackG_x[trackG_N]=(*GretinaCollection)[i]->GetPos().getX()/mm;
	      trackG_y[trackG_N]=(*GretinaCollection)[i]->GetPos().getY()/mm;
	      trackG_z[trackG_N]=(*GretinaCollection)[i]->GetPos().getZ()/mm;
	      trackG_E[trackG_N]=(*GretinaCollection)[i]->GetEdep()/keV;
	      trackG_N++;
	    }else {
#ifdef DEBUG
	      G4cout <<  " \n ========> " <<  G4endl; 
	      G4cout <<  " \t Particle :" << (*GretinaCollection)[i]->GetParticleID() <<  G4endl; 
	      G4cout <<  " \t Energy   : " <<(*GretinaCollection)[i]->GetEdep()/keV << "keV" <<  G4endl; 
	      G4cout <<  " \t Pos X: " << (*GretinaCollection)[i]->GetPos().getX()/mm  << " Y: " << (*GretinaCollection)[i]->GetPos().getY()/mm  << " Z: " << (*GretinaCollection)[i]->GetPos().getZ()/mm   <<  G4endl; 
	      G4cout <<  " \t Cluster : " << c  << " Detector: " <<  d <<  G4endl;
#endif		 	
	    }
	  }
	  if(strcmp((*GretinaCollection)[i]->GetParticleID(),"gamma") == 0  || 
	     strcmp((*GretinaCollection)[i]->GetParticleID(),"e-") == 0 )	
	    AddGretinaGeEnergy((*GretinaCollection)[i]->GetEdep()/keV,c,d);

	  //AR New in v4.3 -> Finding interaction point with the highest energy for Rob's Cascade Analysis
	  if(CASCADE_ANALYSIS_FLAG){
	    if((*GretinaCollection)[i]->GetEdep()>Highesten){
	      Highesten=(*GretinaCollection)[i]->GetEdep();
	      highenmult=i;
	    }
	  }
	  
	} //AR end of condition on disabled Crystal  (New in v4.3)
      }
		
      // // Sort by Mutiplicty
      for (G4int cl=0;cl<MaxClusG;cl++){
		  for (G4int det=0;det<MaxDetG;det++){
			 if( HitPatternG[cl][det] != 0 ){
				DetClNr_G[DetM_G] = cl;
				DetNr_G[DetM_G] = det;
				DetE_G[DetM_G] = ECr_G[cl][det];
				
				if(FWHM_FLAG){
				
				  DetE_G[DetM_G]= detector->GetGretina()->FWHM_Response(DetE_G[DetM_G]);
				  
								
				}

				
				
				
				
				// Calorimeter
				TotalE_G += DetE_G[DetM_G];
				
				// Check if Doppler Reconstruction is required
				if(!detector->GetGretina()->DopplerOff()){
				  
				  // Position of the Heavy Ion
				  switch(DC_IonPos){
				  case 0 : // Mid-target
					 Ion_Pos = Target_Pos;
					 break;
				  case 1 : // True decay Position - In Case of Cascade decay, last decay position
					 Ion_Pos.setX(ion_dec.x[DecayM-1]*mm);
					 Ion_Pos.setY(ion_dec.y[DecayM-1]*mm);
					 Ion_Pos.setZ(ion_dec.z[DecayM-1]*mm);
					 break;
				  case 2 : // User Defined Position
					 Ion_Pos=BeamOut->GetPosDopp();
					 break;
				  case 3 : // Mid-Degrader
					 Ion_Pos = Degrader_Pos;
					 break;
				  case 4 : // Mid-Stopper
					 Ion_Pos = Stopper_Pos;
					 break;						
				  }
			 
				  // Direction of the Heavy Ion
				  switch(DC_IonAxis){
				  case 0 : // z axis
					 Ion_V.setZ(1.);
					 Ion_V.setY(0.);
					 Ion_V.setX(0.);
					 break;
				  case 1 : // Incomming Beam
					 Ion_V.setZ(1./sqrt(1.+tan(BeamIn->getAta0())*tan(BeamIn->getAta0())+tan(BeamIn->getBta0())*tan(BeamIn->getBta0())));
					 Ion_V.setY(tan(BeamIn->getBta0())*1./sqrt(1.+tan(BeamIn->getAta0())*tan(BeamIn->getAta0())+tan(BeamIn->getBta0())*tan(BeamIn->getBta0())));
					 Ion_V.setX(tan(BeamIn->getAta0())*1./sqrt(1.+tan(BeamIn->getAta0())*tan(BeamIn->getAta0())+tan(BeamIn->getBta0())*tan(BeamIn->getBta0())));
					 break;
				  case 2 : // Decay - In Case of Cascade decay, last decay position
					 Ion_V.setZ(1.);
					 Ion_V.setTheta(ion_dec.theta[DecayM-1]);
					 Ion_V.setPhi(ion_dec.phi[DecayM-1]);
					 break;				
				  case 3 : // Target back
					 Ion_V.setZ(1.);
					 Ion_V.setTheta(ion_ev[TARGET_BACK_FLAG].theta);
					 Ion_V.setPhi(ion_ev[TARGET_BACK_FLAG].phi);
					 break;	
				  case 4 : // Degrader Back
					 Ion_V.setZ(1.);
					 Ion_V.setTheta(ion_ev[DEGRADER_BACK_FLAG].theta);
					 Ion_V.setPhi(ion_ev[DEGRADER_BACK_FLAG].phi);
					 break;				
				  case 5 : // Stopper Back
					 Ion_V.setZ(1.);
					 Ion_V.setTheta(ion_ev[STOPPER_BACK_FLAG].theta);
					 Ion_V.setPhi(ion_ev[STOPPER_BACK_FLAG].phi);
					 break;	
				  }
				}
				
				bool Trackfound = false;
				if(GRETINA_TRACK_FLAG == TRUE){			  
				  for(G4int n=trackG_N-1;n>=0;n--){
#ifdef DEBUG
				    G4cout << n << " -- DetNr (track) " << trackG_D[n] << " DetNr (Curr) " << DetNr_G[DetM_G] << "-  Cluster (Track) " << trackG_C[n] << " - ClusterNr (Curr) " <<  DetClNr_G[DetM_G]<<  " X: " <<trackG_x[n]*mm << " Y : " << trackG_y[n]*mm << " Z : " << trackG_z[n]*mm <<  G4endl;
#endif		
				    if(trackG_D[n]==DetNr_G[DetM_G]&&trackG_C[n]==DetClNr_G[DetM_G]){
				      Gamma_TruePos.setX(trackG_x[n]*mm);
				      Gamma_TruePos.setY(trackG_y[n]*mm);
				      Gamma_TruePos.setZ(trackG_z[n]*mm);
				      Trackfound = true;
				    }
				  }
				  
				  if(!Trackfound){
					 MissCtr++;
					 G4cerr << "\n (2) Missing Track at event # " << fEvtNum<< " -- Skipping Event" <<endl;
					 continue;
				  }
				}
				
				// Position resolution ... 
				if(DC_Gretina_Pos_Sigma >0){
				  // Random X,Y,Z
				  Gamma_Pos.setX(G4RandGauss::shoot(Gamma_TruePos.getX()/mm, DC_Gretina_Pos_Sigma )*mm);
				  Gamma_Pos.setY(G4RandGauss::shoot(Gamma_TruePos.getY()/mm, DC_Gretina_Pos_Sigma )*mm);
				  Gamma_Pos.setZ(G4RandGauss::shoot(Gamma_TruePos.getZ()/mm, DC_Gretina_Pos_Sigma )*mm);
				}else{
				  Gamma_Pos = Gamma_TruePos;					 
				}
				
				
#ifdef DEBUG
				G4cout << " Simulated Position          " << Gamma_Pos.getX()/cm << "cm " << Gamma_Pos.getY()/cm  << "cm " << Gamma_Pos.getZ()/cm << "cm - " << Gamma_Pos.getR()/cm  << "cm " << Gamma_Pos.getTheta()/deg  << "deg " << Gamma_Pos.getPhi()/deg  << "deg "  << endl; 
				G4cout << " Real Position               " << Gamma_TruePos.getX()/cm << "cm " << Gamma_TruePos.getY()/cm  << "cm " << Gamma_TruePos.getZ()/cm << "cm - " << Gamma_TruePos.getR()/cm  << "cm " << Gamma_TruePos.getTheta()/deg  << "deg " << Gamma_TruePos.getPhi()/deg  << "deg "  << endl; 
#endif	
				if(!detector->GetGretina()->DopplerOff()){
				  //direction of the gamma-ray
				  Gamma_V=Gamma_Pos-Ion_Pos;
				 
				  //cos of the angle for doppler corrections
				  Gamma_V.setMag(1.);
				  Ion_V.setMag(1.);
				  G4double cosSeg=Ion_V.dot(Gamma_V);
					 
				  //direction of the gamma-ray
				  if(GRETINA_TRACK_FLAG == TRUE){
					 Gamma_VTrue = Gamma_TruePos-Ion_Pos;
					 Gamma_VTrue.setMag(1.);
				  }
					 
				  // cos angle of the real position				
				  G4double	cosSegTrue = Ion_V.dot(Gamma_VTrue);
					 
					 
				  // Velocity of the source for doppler corrections
				  switch(DC_Beta){
				  case 0 : // Use user defined Beta
				    beta_r=BeamOut->GetBetaDopp();
				    gamma_r=BeamOut->GetGammaDopp();
				    break;
				  case 1 : // Use beta after target 
				    beta_r=ion_ev[TARGET_BACK_FLAG].beta;
				    gamma_r=1./sqrt(1-beta_r*beta_r);
				    break;
				  case 2 : // Use True beta at decay - In Case of Cascade decay, last decay position
				    beta_r=ion_dec.beta[DecayM-1];
				    gamma_r=1./sqrt(1-beta_r*beta_r);
				    break;
				  case 3 : // Use beta after degrader
				    beta_r=ion_ev[DEGRADER_BACK_FLAG].beta;
				    gamma_r=1./sqrt(1-beta_r*beta_r);
				    break;
				  case 4 : // Use beta after stopper
				    beta_r=ion_ev[STOPPER_BACK_FLAG].beta;
				    gamma_r=1./sqrt(1-beta_r*beta_r);
				    break;
				  }

				  CosSeg_G[DetM_G] = cosSeg;
				  
				  // 				CosThetaCM[DetM] = (-beta_r*gamma_r + gamma_r * cosSeg)/(gamma_r - beta_r*gamma_r*cosSeg);
				  // 				// True Angle in CM of the gamma
				  // 				CosThetaCMTrue[DetM] = (-beta_r*gamma_r + gamma_r * cosSegTrue)/(gamma_r - beta_r*gamma_r*cosSegTrue);
				  DetEDC_G[DetM_G] = DetE_G[DetM_G] * gamma_r * (1-beta_r*cosSeg);
				  
				  //TM added in version 4.2 for background simulation			
				  G4double BackgroundGammaProb=G4UniformRand();				
				  vector<EnergyBackground_Level*> EnergyBackgroundInfo = BeamOut->GetBackgroundInfo();
				  NumOfBackgroundStates = EnergyBackgroundInfo.size();

				  for( int i = 0; i < NumOfBackgroundStates; i++ )
				    {
				      Energy[i] = EnergyBackgroundInfo[i]->GetEb()*1.0e3;
				      BackgroundFrac[i] = EnergyBackgroundInfo[i]->GetFracb()*0.01;	
				      DetE_Background_G[i] = Energy[i];
				      if(BackgroundGammaProb>0 and BackgroundGammaProb< BackgroundFrac[i]){	
					DetE_Background_G[i] = detector->GetGretina()->FWHM_Response( DetE_Background_G[i] ); 
					DetE_Background_G[DetM_G] =DetE_Background_G[i];

				      }	
				    }
				  DetEDC_Background_G[DetM_G] =   DetE_Background_G[DetM_G] * gamma_r * (1-beta_r*cosSeg);
				 
					
				  //TM added for version 4.2 Charlie's ad				
				  //				  if(DetM_G==0){Highesten=DetE_G[0];}
				  if(DetE_G[DetM_G]>Highesten){ //CRL333 This searches through all the events in a given multiplicity, and uses the global variables below to determine information for doppler shifting and angle cuts
				    if(!CASCADE_ANALYSIS_FLAG){
				      Highesten=DetE_G[DetM_G];     //As before the highest energy
				      addbackbeta=beta_r;           //This contains the beta needed for doppler shifting
				      addbackcos=CosSeg_G[DetM_G];  //This contains the cosine theta needed for doppler shifting
				      highenmult=DetM_G;            //This contains the array location of the event with highest energy
				    }else{
				      addbackbeta=beta_r;           //This contains the beta needed for doppler shifting
				    }
				  }
				  //////////////////////
				}
				
				
				
				if(GRETINA_TRACK_FLAG == TRUE){
				  // R_G[DetM_G] = (Gamma_Pos-Target_Pos).getR()/cm;
				  // Theta_G[DetM_G] = (Gamma_Pos-Target_Pos).getTheta()/deg;
				  // Phi_G[DetM_G] = (Gamma_Pos-Target_Pos).getPhi()/deg;
				  // dR_G[DetM] = (Gamma_Pos.getR()-Gamma_TruePos.getR())/mm;
				  // dTheta_G[DetM] = (Gamma_Pos.getTheta()-Gamma_TruePos.getTheta())/deg;
				  // dPhi_G[DetM] = (Gamma_Pos.getPhi()-Gamma_TruePos.getPhi())/deg;
				  // X_G[DetM_G] = (Gamma_Pos-Target_Pos).getX()/cm;
				  // Y_G[DetM_G] = (Gamma_Pos-Target_Pos).getY()/cm;
				  // Z_G[DetM_G] = (Gamma_Pos-Target_Pos).getZ()/cm;
				  
				  //AR New in v4.3
				  R_G[DetM_G] = (Gamma_Pos-Ion_Pos).getR()/cm;
				  Theta_G[DetM_G] = (Gamma_Pos-Ion_Pos).getTheta()/deg;
				  Phi_G[DetM_G] = (Gamma_Pos-Ion_Pos).getPhi()/deg;
				  Depth_G[DetM_G] = GetDepth(Gamma_Pos,Ion_Pos.getZ()); //AR New in v4.3 -> Computing Depth
				  dR_G[DetM] = (Gamma_Pos.getR()-Gamma_TruePos.getR())/mm;
				  dTheta_G[DetM] = (Gamma_Pos.getTheta()-Gamma_TruePos.getTheta())/deg;
				  dPhi_G[DetM] = (Gamma_Pos.getPhi()-Gamma_TruePos.getPhi())/deg;
				  X_G[DetM_G] = (Gamma_Pos-Ion_Pos).getX()/cm;
				  Y_G[DetM_G] = (Gamma_Pos-Ion_Pos).getY()/cm;
				  Z_G[DetM_G] = (Gamma_Pos-Ion_Pos).getZ()/cm;
				}
				

				DetM_G++;
			 }
		  }
      }
    }
  }

  //AR New in v4.3 -> To perform Rob's Cascade Analysis
  //Calculation of CascadeZ and Corrected0pE, the 
  // values for the new histograms for Rob's Analysis.
  if(CASCADE_ANALYSIS_FLAG){
    //CorrectedE will get a value for each gamma seen in coincidence with the known gamma.  Thus it is a vector in this system.
    vector<G4double> empty;

    GammaKnown_Low = Casc_Ana_ELow;
    GammaKnown_High = Casc_Ana_EHigh;
    TrueE_GammaKnown = Casc_Ana_ETrue;

    for(int i=0;i<DetM_G;i++){
    
      CorrectedE.push_back(empty);
      if(DetE_G[i]>GammaKnown_Low && DetE_G[i]<GammaKnown_High)
	{
	  LabE_GammaKnown = DetE_G[i];
	  Z_Tar_GammaKnown = Z_G[i]*10; //now in mm.
	  R_Tar_GammaKnown = R_G[i]*10; //now in mm.  
	  Rho_Tar_GammaKnown = sqrt(R_Tar_GammaKnown*R_Tar_GammaKnown-Z_Tar_GammaKnown*Z_Tar_GammaKnown);
	  Z_Origin_Tar = Target_Pos.getZ()/mm;  //in mm

	  RelGamma = 1/sqrt(1-beta_r*beta_r);
	  Theta_GammaKnown = acos((1-(TrueE_GammaKnown)/(RelGamma*LabE_GammaKnown))/beta_r);
	  Z_Decay_GammaKnown = Rho_Tar_GammaKnown/tan(Theta_GammaKnown);
	  Z_Tar_Decay = Z_Tar_GammaKnown-Z_Decay_GammaKnown;
	  CascadeZ.push_back(Z_Tar_Decay);

	  for(int j=0;j<DetM_G;j++)
	    {
	      if(j!=i)
		{
		  LabE_GammaOfInterest = DetE_G[j];
		  Z_Tar_GammaOfInterest = Z_G[j]*10;
		  R_Tar_GammaOfInterest = R_G[j]*10;
		  Rho_Tar_GammaOfInterest = sqrt(R_Tar_GammaOfInterest*R_Tar_GammaOfInterest-Z_Tar_GammaOfInterest*Z_Tar_GammaOfInterest);
		  Z_Decay_GammaOfInterest = -CascadeZ.back()+Z_Tar_GammaOfInterest;
		  Theta_GammaOfInterest = atan(Rho_Tar_GammaOfInterest/Z_Decay_GammaOfInterest);
		  // atan is arbitrary for its return value.  
		  // If it returns a value between -pi and 0, 
		  // we change it into the corresponding 
		  // value between 0 and +pi.  
		  if(Theta_GammaOfInterest<0)
		    {
		      Theta_GammaOfInterest = Theta_GammaOfInterest+pi;
		    }
		  CorrectedE[i].push_back(LabE_GammaOfInterest*RelGamma*(1-beta_r*cos(Theta_GammaOfInterest)));

		}
	      else
		{
		  CorrectedE[i].push_back(-1000.);
		}
	    }
	  
	}
      else
	{
	  CascadeZ.push_back(-1000);
	}
    }
  }
}


//AR New in v4.3 -> To determine the depth parameter
G4double Analysis::GetDepth(G4ThreeVector pos, G4double Tz){
  Tz=0.;
  G4double depth = 0;

  G4double R_G = 185.; //Gretina radius

  G4double a = pos.getX()*pos.getX()+pos.getY()*pos.getY()+(pos.getZ()-Tz)*(pos.getZ()-Tz);
  G4double b = 2*((pos.getZ()-Tz)*Tz);
  G4double c = Tz*Tz - R_G*R_G;
  G4double u1 = (-b + sqrt(b*b-4*a*c))/(2*a);
  G4double u2 = (-b - sqrt(b*b-4*a*c))/(2*a);
  G4ThreeVector pos1, pos2; // point of intersection gretina sphere and gamma trajectory

  pos1.set(u1*pos.getX(),u1*pos.getY(),Tz + u1*(pos.getZ()-Tz));
  pos2.set(u2*pos.getX(),u2*pos.getY(),Tz + u2*(pos.getZ()-Tz));

  if((pos-pos1).mag()<(pos-pos2).mag()){
    depth = (pos-pos1).mag();
  }else{
    depth = (pos-pos2).mag();
  }
  
  return depth;
}

// TM added for sphere addback in v4.2
void Analysis::AddbackMode3Simple(void){
  
  //ORG const float MaxDist=7.5*cm; //This sets the distance between the two points which will allow for summing.  
  MaxDist=Analysis::GetMaxDist(); 
  const  int MaxDet=DetM_G;
  int maxpos=highenmult;
  int numtested=0;
  vector<int> mult;
  vector<vector<int> > tested; //This is a vector which holds each set of points that are to be summed. 
  while(numtested!=MaxDet){
    int multinsum=0;
    vector<int> summed;
    summed.push_back(maxpos);
    multinsum++;
    numtested++;

    int j=0;

    while(j<MaxDet){

      if (j!=maxpos&& GetDistance(maxpos,j)<MaxDist){ //Simply tests to see if a point already in the region is within some distance from another point
	//	cout <<"Distance is " <<  GetDistance(highenmult,j) << endl;
	unsigned int k=0;
	bool untested=true;//
	while(k<tested.size()){
	  unsigned int m=0;
	  while(m<tested.at(k).size()){
	    if (j==tested.at(k).at(m)){untested=false;cout << "false"<< endl;}
	    m++;}//Here if the hit is already in the vector we skip the hit
	k++;}
		if(untested){//If the hit is not in the tested vector we add it in so that it will also be looped over		  	 
		  summed.push_back(j);multinsum++;numtested++;}
       }	
     j++;
     }
    //    cout << numtested << endl;
    //    cout << MaxDet;
    tested.push_back(summed);
    mult.push_back(multinsum);
    maxpos=GetLargestRemaining(tested);
    
  }

  for(unsigned int  it=0; it!=tested.size(); it++){
    G4double ensum=0;
    for(unsigned int it2=0; it2!=tested.at(it).size();it2++){
      ensum+=DetE_G[tested.at(it).at(it2)];}

    ensum=ensum* 1.0/(sqrt(1-pow(addbackbeta,2))) * (1-addbackbeta*CosSeg_G[tested.at(it).at(0)]);

    AddbackEdcModeSphere.push_back(ensum);
    AddbackAngleModeSphere.push_back(Theta_G[tested.at(it).at(0)]);
    AddbackMult.push_back(mult.at(it));//Number of hits in the sphere
  }
  //  cout << "the multiplicty is " << DetM_G<< "  ";
  //  cout << "The size of the addback is " << AddbackEdcModeSphere.size()<<endl;
 
  //  AddbackEdcMode3n  = AddbackEdcMode3n * 1.0/(sqrt(1-pow(addbackbeta,2))) * (1-addbackbeta*addbackcos);//Finally we use the position information of the highest multiplicty event and the beta setting to determine the doppler corrected energies
  //  Need to get beta and cosine for each position. Should be possible easily because that position will be at 0;  
  //  cout << "finish again" << endl;


}

//AR New for v4.3 -> Used for Rob's Cascade Analysis
void Analysis::AddbackMode3Simple(GretinaHitDetectorCollection* ghdc){

  int ipMult = ghdc->entries();
  
  const float MaxDist=8.0*cm;//This sets the distance between the two points which will allow for summing.  EDIT HERE
  const  int MaxDet=ipMult;
  int maxpos=highenmult;
  int numtested=0;
  vector<int> mult;
  vector<vector<int> > tested;//This is a vector which holds each set of points that are to be summed.

  while(numtested!=MaxDet){
    int multinsum=0;
    vector<int> summed;
    summed.push_back(maxpos);
    multinsum++;
    numtested++;
    int j=0;
    
    while(j<MaxDet){
      if (j!=maxpos && GetDistance((*ghdc)[maxpos]->GetPos(),(*ghdc)[j]->GetPos())<MaxDist){//Simply tests to see if a point already in the region is within some distance from another point
	unsigned int k=0;
	bool untested=true;//
	
	while(k<tested.size()){
	  unsigned int m=0;
	  while(m<tested.at(k).size()){
	    if (j==tested.at(k).at(m)){
	      untested=false;
	      //	      cout << "false"<< endl;
	    }m++;
	  }//Here if the hit is already in the vector we skip the hit
	  k++;
	}
	if(untested){//If the hit is not in the tested vector we add it in so that it will also be looped over
	  
	  summed.push_back(j);
	  multinsum++;
	  numtested++;
	}
      }
      j++;
    }
    tested.push_back(summed);
    mult.push_back(multinsum);
    maxpos=GetLargestRemaining(tested,ghdc);
  }

  for(unsigned int  it=0; it<tested.size(); it++){
    G4double ensum=0;
    for(unsigned int it2=0; it2!=tested.at(it).size();it2++){
      ensum+=(*ghdc)[tested.at(it).at(it2)]->GetEdep()/keV;
    }
    AddbackEModeSphere.push_back(ensum);
    //Rob Elder 2018-03-06.  Adding the FWHM response that is used for the "no AB" method above.
    if(FWHM_FLAG){
      AddbackEModeSphere[it] = detector->GetGretina()->FWHM_Response(AddbackEModeSphere[it]);
    }

    //Rob Elder 2018-03-06.  Adding pos res.
    AddbackPos.push_back((*ghdc)[tested.at(it).at(0)]->GetPos());

    if(DC_Gretina_Pos_Sigma >0){
      // Random X,Y,Z
      AddbackPos[it].setX(G4RandGauss::shoot(AddbackPos[it].getX()/mm, DC_Gretina_Pos_Sigma )*mm);
      AddbackPos[it].setY(G4RandGauss::shoot(AddbackPos[it].getY()/mm, DC_Gretina_Pos_Sigma )*mm);
      AddbackPos[it].setZ(G4RandGauss::shoot(AddbackPos[it].getZ()/mm, DC_Gretina_Pos_Sigma )*mm);
    }
    else{
      ;
    }

    ensum=ensum* 1.0/(sqrt(1-pow(addbackbeta,2))) * (1-addbackbeta*GetCosTheta((*ghdc)[tested.at(it).at(0)]->GetPos()));

    AddbackEdcModeSphere.push_back(ensum);
    AddbackAngleModeSphere.push_back(acos(GetCosTheta((*ghdc)[tested.at(it).at(0)]->GetPos())));
    AddbackMult.push_back(mult.at(it));//Number of hits in the sphere
    
  }

  //Rob Elder 2018-02-24.  Calculating the radius of the sphere containing a set of AB'ed points centered on the point
  //with the highest energy (also called the "span" in other codes).
  for(unsigned int i=0; i<tested.size(); i++){

    G4double span = CalcSpan(tested[i],ghdc);
    ABsphereRad.push_back(span);

  }

  //Rob Elder 2018-02-16.  Adding my Cascade DC method to the sphere AB method.
  vector<G4double> empty;
  
  //Rob Elder 2019-04-05.  I want to check how many gammas detected by GRETINA are consistent with "known gamma"
  //Then, if there are more than one, choose just one to use to find the decay position.
  int multABgoodgk = 0;
  vector<int> listGoodgk;
  int chosenGoodgk = -1;
  for(unsigned int i=0;i<AddbackEModeSphere.size();i++){
    double thisEnergy = AddbackEModeSphere[i];
    if(thisEnergy > GammaKnown_Low && thisEnergy < GammaKnown_High){
      multABgoodgk++;
      listGoodgk.push_back(i);
    }
  }
  if(multABgoodgk>1){
    int randInt = rand() % multABgoodgk;
    chosenGoodgk = listGoodgk[randInt];
  }
  else if(multABgoodgk==1){
    chosenGoodgk = listGoodgk[0];
  }
  else{
    ;
  }
    
  
  for(unsigned int i=0;i<AddbackEModeSphere.size();i++)
    {
  
      CorrectedE_AB.push_back(empty);
      
      if(AddbackEModeSphere[i]>GammaKnown_Low && AddbackEModeSphere[i]<GammaKnown_High)
	{
	  //Rob Elder 2019-04-05.  Now only find the decay position if "i" is the gamma that I chose to use earlier.  
	if(i==chosenGoodgk){
	  LabE_GammaKnown = AddbackEModeSphere[i];
	  Z_Origin_Tar = Target_Pos.getZ()/mm;
	  Z_Origin_GammaKnown = AddbackPos[i].z();
	  R_Origin_GammaKnown = AddbackPos[i].r();
	  Rho_Tar_GammaKnown = sqrt(R_Origin_GammaKnown*R_Origin_GammaKnown-Z_Origin_GammaKnown*Z_Origin_GammaKnown);

	  RelGamma = 1/sqrt(1-beta_r*beta_r);
	  Theta_GammaKnown = acos((1-(TrueE_GammaKnown)/(RelGamma*LabE_GammaKnown))/beta_r);
	  Z_Decay_GammaKnown = Rho_Tar_GammaKnown/tan(Theta_GammaKnown);
	  Z_Tar_Decay = Z_Origin_GammaKnown-Z_Decay_GammaKnown-Z_Origin_Tar;     
	  CascadeZ_AB.push_back(Z_Tar_Decay);

	  //Rob Elder 2018-12-17.  Moved the comptAngDiff calculation to occur outside of the cascade analysis section.
	  //the comptAngDiff calc now occurs just after the span calc.
	  //Rob Elder 2018-12-20.  Returned the comptAngDiff calculation to this location because it requires CascadeZ to
	  //have been calculated.  
	  //Rob Elder 2018-12-20.  Need to run the actual calculation for the truly simulated gammas, but just sample
	  //the 35P background distribution for the bgd gammas.  

	  if(i<tested.size()){ //then this is still a fully simulated gamma
	    if(tested[i].size()>1){
	      double comptAngle = acos(CosComptAngle((*ghdc)[tested[i][0]]->GetEdep()/keV,AddbackEModeSphere[i]-(*ghdc)[tested[i][0]]->GetEdep()/keV));
	      G4ThreeVector decay(0.0,0.0,CascadeZ_AB.back());
	      double geomAngle = CasZangle((*ghdc)[tested[i][0]]->GetPos(),(*ghdc)[tested[i][1]]->GetPos(),decay);
	      ABcomptAngDiff.push_back(abs(comptAngle-geomAngle));
	    }
	    else{
	      ABcomptAngDiff.push_back(0.0);
	    }
	  }
	  else{   
	  }
	  
	  for(unsigned int j=0;j<AddbackEModeSphere.size();j++)
	    {
	      if(j!=i)
		{
		  LabE_GammaOfInterest = AddbackEModeSphere[j];
		  Z_Origin_GammaOfInterest = AddbackPos[j].z();
		  R_Origin_GammaOfInterest = AddbackPos[j].r();
		  Rho_Tar_GammaOfInterest = sqrt(R_Origin_GammaOfInterest*R_Origin_GammaOfInterest-Z_Origin_GammaOfInterest*Z_Origin_GammaOfInterest);
		  Z_Decay_GammaOfInterest = Z_Origin_GammaOfInterest-Z_Origin_Tar-Z_Tar_Decay;
		  Theta_GammaOfInterest = atan(Rho_Tar_GammaOfInterest/Z_Decay_GammaOfInterest);
		  // atan is arbitrary for its return value.  
		  // If it returns a value between -pi and 0, 
		  // we change it into the corresponding 
		  // value between 0 and +pi.  
		  if(Theta_GammaOfInterest<0)
		    {
		      Theta_GammaOfInterest = Theta_GammaOfInterest+pi;
		    }
		  CorrectedE_AB[i].push_back(LabE_GammaOfInterest*RelGamma*(1-beta_r*cos(Theta_GammaOfInterest)));
		  
		}
	      else
		{
		  CorrectedE_AB[i].push_back(-1000.);
		}
	    }
	  }//close "if(i==chosenGoodgk)"
	  
	}
      else
	{
	  CascadeZ_AB.push_back(-1000);
	  ABcomptAngDiff.push_back(10.0);  
	}
    }

}

void Analysis::AddESegment(G4int r, G4int d, G4int s, G4int q, G4double E)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif  
  ESeg[r][d][s][q] += E;
}

void Analysis::Report(void)
{
#ifdef DEBUG
  G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif  

  G4cout << "-------------------------------------------" <<  endl; 
  G4cout << "  Reporting on Analyis : " <<  endl;
  G4cout << "-------------------------------------------" <<  endl; 
  G4cout << " S800_ACCEPTANCE : " << S800_ACCEPTANCE_FLAG  << endl;
  if(S800_ACCEPTANCE_FLAG){
	 G4cout << " --> Momentum : [" << AccDL*100. << ":" << AccDH*100. << "]%"<< endl;
	 G4cout << " --> Angle : [" << AccThL/mrad << ":" << AccThH /mrad<< "] rad"<< endl;
  }
  G4cout << " SeGA FWHM : " << FWHM_FLAG  << endl;
  G4cout << " Doppler Correction : " << endl;
  G4cout << " --> Beta : " << DC_Beta << " (0: User; 1: After Target; 2: True (decay); 3: After Degrader; 4: After Stopper;)"<< endl;
  G4cout << " --> Ion Axis : " << DC_IonAxis << " (0: Z; 1: Beam In ; 2: True (decay); 3: After Target; 4: After Degrader; 5: After Stopper;)" << endl;
  G4cout << " --> Ion Position: " << DC_IonPos<< " (0: Mid-Target; 1: True (decay); 2: User; 3: Mid-Degrader; 4: Mid-Stopper)" <<endl;
  G4cout << " --> Gamma Position : " <<DC_Gamma_Pos<< " (0: Center of Segments; 1: True ; )" << endl;

}
//---------------------------------------------------------------------
G4bool Analysis::S800Accept(ion_event_t evt)
{
#ifdef DEBUG
	G4cout << __PRETTY_FUNCTION__  << G4endl;
	G4cout << "S800 Acceptance Flag : " <<  S800_ACCEPTANCE_FLAG  << G4endl;

#endif  
	if(S800_ACCEPTANCE_FLAG) {	
	  // Momentum acceptance
	  G4ThreeVector cent,ion;
	  G4double     theta,a,b,d,x,y,z;
	  d=(evt.KE/BeamOut->GetURsetKE()-1.);
	  //d=CLHEP::RandGauss::shoot(d,sig_d);
	  //d/=100.;
	  if(d<AccDL||d>AccDH) return false;
	  
	  
	  // Angular acceptance
	  a=BeamIn->getAta0();
	  b=BeamIn->getBta0();
	  z=1./sqrt(1.+tan(a)*tan(a)+tan(b)*tan(b));
	  y=tan(b)*z;
	  x=tan(a)*z;
	  cent.setZ(z);
	  cent.setY(y);
	  cent.setX(x);
	  cent.setMag(1.);
	  ion.setZ(1.);
	  ion.setTheta(evt.theta);
	  ion.setPhi(evt.phi);
	  theta=acos(cent.dot(ion));
	  if(theta<AccThL||theta>AccThH) return false;
	}	  
	return true;

}

void Analysis::AddGretinaGeEnergy( G4double energy, G4int clunum, G4int detnum)
{
#ifdef DEBUG
	G4cout << __PRETTY_FUNCTION__  << G4endl;
#endif  
	ECr_G[clunum][detnum] += energy;
}


// TM added for sphere addback v4.2

float Analysis::GetDistance(int i, int j){//This simply gets the distance between two hit points in standard cartesian style
      float distance=(X_G[i]-X_G[j])*(X_G[i]-X_G[j]);
      distance=distance+((Y_G[i]-Y_G[j])*(Y_G[i]-Y_G[j]));
      distance=distance+(Z_G[i]-Z_G[j])*(Z_G[i]-Z_G[j]);
      return sqrt(distance);
    }

int Analysis::GetLargestRemaining(vector<vector<int> > tested){//Iterate over vector and compare to see what largest energy is
    int i=-1;
    const  int MaxDet=DetM_G;
    float high=0;
    bool untest=true;
    int highpoint=-1;
    while (i< MaxDet){
	untest=true;
	  if (DetE_G[i]>high){
	    unsigned int j=0;
	    while(j<tested.size()){
	      unsigned int k=0;
	      while(k<tested.at(j).size()){
		if (tested.at(j).at(k)==i){untest=false;}k++;}j++;}
	  if(untest){
	    high=DetE_G[i];highpoint=i;}}
	i++;}

    return highpoint;
}

//AR Added for v4.3 - Used for Rob's Cascade Analysis
//Making a new GetDistance method that takes two GretinaHitDetectors' positions
float Analysis::GetDistance(G4ThreeVector p1, G4ThreeVector p2){//This simply gets the distance between two hit points in standard cartesian style
  float distance=(p1.x()-p2.x())*(p1.x()-p2.x());
  distance=distance+((p1.y()-p2.y())*(p1.y()-p2.y()));
  distance=distance+((p1.z()-p2.z())*(p1.z()-p2.z()));
  return sqrt(distance);
}

//Making this function to work with the GretinaHitDetectorCollection instead of the list of DetE_G.  
int Analysis::GetLargestRemaining(vector<vector<int> > tested, GretinaHitDetectorCollection* ghdc){//Iterate over vector and compare to see what largest energy is
  int i=0;
  const  int MaxDet=ghdc->entries();
  float high=-1.0;
  bool untest=true;
  int highpoint=-1;
  while (i< MaxDet){
    untest=true;
    //    cout << "(*ghdc)["<<i<<"]->GetEdep(): " << (*ghdc)[i]->GetEdep() << endl;
    if ((*ghdc)[i]->GetEdep()>high){
      unsigned int j=0;
      while(j<tested.size()){
	unsigned int k=0;
	while(k<tested.at(j).size()){
	  //	  cout << "tested.at("<<j<<").at("<<k<<"): " << tested.at(j).at(k) << endl;
	  if (tested.at(j).at(k)==i){
	    untest=false;
	  }
	  k++;
	}
	j++;
      }
      if(untest){
	high=(*ghdc)[i]->GetEdep();
	highpoint=i;
      }
    }
    i++;
  }
  return highpoint;
}

//Function to find cos of angle between ion trajectory and gamma based on a GretinaHitDetector position
double Analysis::GetCosTheta(G4ThreeVector gpos){
  G4ThreeVector gammaV = gpos-Ion_Pos;
  gammaV.setMag(1.);
  G4ThreeVector ionV = Ion_V;
  if(!(ionV.x()==0. && ionV.y()==0. && ionV.z()==0.)){
    ionV.setMag(1.);
  }
  G4double costheta=Ion_V.dot(gammaV);
  
  return costheta;
}

//Function to calculate the radius of the sphere enclosing a group of AB points centered on the point with the highest energy (also called the span).
G4double Analysis::CalcSpan(vector<int> abPoint, GretinaHitDetectorCollection* ghdc){
  G4double span;

  if(abPoint.size()==1){
    span = 0.0;
  }
  else if(abPoint.size()>1){
    double maxDist = 0.0;
    int maxIndex = 0;
    for(unsigned int i=1;i<abPoint.size();i++){
      double thisDist = ((*ghdc)[abPoint[i]]->GetPos()-(*ghdc)[abPoint[0]]->GetPos()).mag();
      if(thisDist>maxDist){
	maxDist=thisDist;
	maxIndex=i;
      }
    }
    span=maxDist;
  }
  else{
    span = -1.0;
  }
  return span;
}

// This calculates the cosine of the angle between an incoming gamma-ray and the compton scattered
//gamma-ray.
//enDep is the energy deposited in the detector (ie, transferred to the electron)
//enRem is the energy remaining with the gamma-ray after the Compton Scatter.
//enDep + enRem is then the initial energy of the gamma-ray.  
double Analysis::CosComptAngle(double enDep, double enRem){
  double cosAngle = 1-511.0/enRem+511.0/(enDep+enRem);
  return cosAngle;
}

//This takes a gamma-ray interaction position and the next gamma-ray interaction position, as well as the
//decay position that originated the gamma-ray, and calculates the angle between the incoming and outgoing gamma-ray at the
//first interaction point.  
double Analysis::CasZangle(G4ThreeVector posCompt, G4ThreeVector posNext, const G4ThreeVector& decay){
  //the IPoint positions are initially w/rt the target.  transform them to be w/rt decay position.
  posCompt-=decay;
  posNext-=decay;
  double casZangle = posCompt.angle(posNext-posCompt);
  return casZangle;
}

/////////////////////   
