////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCutG.h"
#include "TEnv.h"
#include "TKey.h"
#include "TDirectory.h"

#include "CommandLineInterface.hh"
#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "Mode3Calc.hh"
#include "Scaler.hh"
#include "CalHistograms.hh"
using namespace TMath;
using namespace std;

//not nice but works at least
static TCutG* TimeCut;
static bool foundTCut = false;


//RE New in v4.1_LT -> Making function to return the cylindrical rho distance from the beam axis.
double CalHistograms::rho(TVector3 hit){
  return sqrt(hit.X()*hit.X()+hit.Y()*hit.Y());
}

//RE New in v4.1_LT -> Take the positions of the compton scatter interaction and the next interaction (either compt or photoelec) and the position
//of the decay.  calculate and return the angle between the incoming and outgoing compton scatter gammas.
double CalHistograms::CasZangle(TVector3 posCompt, TVector3 posNext,const TVector3& decay){
  posCompt-=decay;
  posNext-=decay;
  double casZangle = posCompt.Angle(posNext-posCompt);
  return casZangle;
}

//RE New in v4.1_LT
//This function is based on the Compton scatter formula. It takes the energy deposited "enDep" and the energy remaining "enRem" after a Compton scatter event and returns "Cos(theta)"
//The "Cos(theta)" from this function is based only on the energy information. It was compared with the angle determined only from gamma-ray position information to judge whether the gamma ray hits were properly reconstructed
double CosComptAngle(double enDep,double enRem){
  double cosAngle = 1-511.0/enRem+511.0/(enDep+enRem);
  return cosAngle;
}

//RE New in v4.1_LT
//This is the function that determines the decay position. It takes an observed gamma ray "gammaKnown", the expected ion-frame energy "eTrueGammaKnown", and the velocity of the ion "beta".
//It returns a vector corresponding to the position of the decay relative to the target
TVector3 CalHistograms::CascadeZdecay(HitCalc* gammaKnown, double eTrueGammaKnown, double beta){
  TVector3 decay(0.,0.,-10.);
  
  double relGamma = 1/sqrt(1-beta*beta);
  
  double eLabGammaKnown = gammaKnown->GetEnergy();

  TVector3 hitGammaKnown = gammaKnown->GetPosition();
  double zHitGammaKnown = hitGammaKnown.Z();
  double xHitGammaKnown = hitGammaKnown.X();
  double yHitGammaKnown = hitGammaKnown.Y();

  double thetaGammaKnown = acos((1-(eTrueGammaKnown)/(relGamma*eLabGammaKnown))/beta);
  double rhoHitGammaKnown = sqrt(xHitGammaKnown*xHitGammaKnown+yHitGammaKnown*yHitGammaKnown);

  double dDecayToHit;     //this is the distance along the z-axis from the decay to the gamma hit.  positive is downstream/  
  if(thetaGammaKnown==Pi()/2.){  //RE New in v4.1_LT -> Including ways to handle the error cases.  
    dDecayToHit = 0.;
  }
  else if(thetaGammaKnown==0){  //gamma was emitted straight ahead.  impossible to detect both ion and gamma
    return decay;                  //for this case in our geometry.  
  }                              
  else if(thetaGammaKnown==Pi()){  //gamma was emitted straight backward.  impossible to detect both ion and gamma
    return decay;                  //for this case in our geometry.  
  }
  else{
    dDecayToHit = rhoHitGammaKnown/tan(thetaGammaKnown);
  }
  
  decay.SetZ(zHitGammaKnown - dDecayToHit);
  
  return decay;
}

//RE New in v4.1_LT
//This function performs the Doppler-shift correction using the decay position from the CascadeZdecay function. It takes the other observed gamma ray "gammaOfInterest", the decay position "decay", and the velocity "beta".
//It returns the Doppler-shift corrected energy
double CalHistograms::CascadeCorrectedE(HitCalc* gammaOfInterest, const TVector3& decay, double beta){
  double relGamma = 1/sqrt(1-beta*beta);

  double eLabGammaOfInterest = gammaOfInterest->GetEnergy();
  TVector3 hitGammaOfInterest = gammaOfInterest->GetPosition();
  double zHitGammaOfInterest = hitGammaOfInterest.Z();
  double xHitGammaOfInterest = hitGammaOfInterest.X();
  double yHitGammaOfInterest = hitGammaOfInterest.Y();

  double rhoHitGammaOfInterest = sqrt(xHitGammaOfInterest*xHitGammaOfInterest+yHitGammaOfInterest*yHitGammaOfInterest);
  double dDecayToHit = zHitGammaOfInterest - decay.Z();
  
  double thetaGammaOfInterest = atan(rhoHitGammaOfInterest/dDecayToHit);

  if(thetaGammaOfInterest<0.){
    thetaGammaOfInterest = thetaGammaOfInterest+Pi();
  }
  else{   //theta is in the good range.  do nothing.  
    ;
  }
  
  double eCorrGammaOfInterest = eLabGammaOfInterest*relGamma*(1-beta*cos(thetaGammaOfInterest));
  
  return eCorrGammaOfInterest;
}

//RE New in v4.1_LT
//Function used for some gates. For one gamma ray, it calculates the distance from the first hit position to the furthest hit position. In the code I call this quantity the "span" but it also makes sense to call it the "maximum radius" of collection of hit positions for one gamma ray
double CalHistograms::CalcSpan(HitCalc* gamma){
  double span;
  
  if(gamma->GetIPMult()==1){
    span = 0.0;
  }
  else if(gamma->GetIPMult()>1){
    double maxDist = 0.0;
    int maxIndex = 0;
    for(int i=0;i<gamma->GetIPMult();i++){
      double thisDist = (gamma->GetIPoints()[i]->GetPosition()-gamma->GetPosition()).Mag();
      if(thisDist>maxDist){
	maxDist=thisDist;
	maxIndex=i;
      }
    }
    span=maxDist;
  }
  else{
    //ip mult is negative for some reason - error.
    span = -1.0;
  }
  //  cout << "span: " << span << endl << endl;
  return span;
}

void CalHistograms::FillHistograms(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c){

  //bool hass800 = s800->GetTimeS800()!=0;
  bool hasgret = gr->GetMult()>0;

  //Declare the histograms here.
  static vector<TCutG*> InPartCut;
  static vector<vector<TCutG*> > OutPartCut;

  static PAD* pad[2];
  static PPAC* ppac[2];
  static IC* ich;
  static STACKPIN* stackpin; //AR New in v4.1_LT
  static SCINT* scint;
  static TOF* tof;
  static TRACK* track;
  static IITRACK* iitrack;
  static HODO* hodo;

  static bool foundCuts = false;
  static bool NoInCuts = false;
  if (!foundCuts){
    //Read in the cuts file for incoming and outgoing particle ID
    if(fhasfile){
      cout << "Cuts were created for ";
      if (ftac&(1<<2))
	cout << "MTDC";
      else if (ftac&(1<<0))
	cout << "TDC";
      else
	cout << "TAC";
      cout << " data and are applied to the ";
      if (ftac&(1<<1))
	cout << "XFP";
      else
	cout << "OBJ";
      cout << " scintillator" << endl;


      //Remember the current directory so that we can revert to it later.
      TDirectory* outfile = gDirectory;

      char* Name = NULL;
      char* Name2 = NULL;
      TFile *cFile = new TFile(fcutfile);
      TIter next(cFile->GetListOfKeys());
      TKey *key;
      while((key=(TKey*)next())){
	if(strcmp(key->GetClassName(),"TCutG") == 0){
	  Name = (char*)key->GetName();
	  if(strstr(Name,"in") && !strstr(Name,"out")){
	    cout << "incoming cut found "<<Name << endl;
	    InPartCut.push_back((TCutG*)cFile->Get(Name));
	  }
	  if(strstr(Name,"tgam")){
	    cout << "timing cut found "<<Name << endl;
	    TimeCut = (TCutG*)cFile->Get(Name);
	    foundTCut = true;
	  }
	}
      }
      TIter next2(cFile->GetListOfKeys());
      if(InPartCut.size()>0){
	OutPartCut.resize(InPartCut.size());
	while((key=(TKey*)next2())){
	  if(strcmp(key->GetClassName(),"TCutG") == 0){
	    Name = (char*)key->GetName();
	    if(strstr(Name,"in") && strstr(Name,"out")){
	      for(unsigned short i=0;i<InPartCut.size();i++){
		Name2 = (char*)InPartCut[i]->GetName();
		if(strstr(Name,strstr(Name2,Name2))){
		  OutPartCut[i].push_back((TCutG*)cFile->Get(Name));
		  cout << "outgoing cut found "<<Name << endl;
		}
	      }
	    }
	  }
	}
      }
      else{
	OutPartCut.resize(1);
	while((key=(TKey*)next2())){
	  if(strcmp(key->GetClassName(),"TCutG") == 0){
	    Name = (char*)key->GetName();
	    if(!strstr(Name,"in") && strstr(Name,"out")){
	      OutPartCut[0].push_back((TCutG*)cFile->Get(Name));
	      cout << "outgoing cut found (no incoming!) "<<Name << endl;
	      NoInCuts = true;
	    }
	  }
	}
      }
      cFile->Close();
      outfile->cd();

      foundCuts = true;
      cout << "Cuts found" << endl;
    }
  }

  //-------------------------------------------------------------------------
  //*************************************************************************
  //Fill the histograms here.
  //*************************************************************************
  //-------------------------------------------------------------------------

  int hp = s800->GetRegistr();
  FillI("registr",20,0,20,hp);
  for(int j=0;j<16;j++){
    if(hp & (1<<j))
      FillI("trigbit",16,0,16,j);
  }

  stackpin = s800->GetSTACKPIN(); //AR New in v4.1_LT
  ich = s800->GetIC();
  tof = s800->GetTOF();
  track = s800->GetTRACK();
  iitrack = s800->GetIITRACK();
  hodo = s800->GetHODO();
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }

  FillHistogramsNoGate(gr,s800,m3c);
  FillHistogramsGateIn(gr,s800,m3c,"all");
  FillHistogramsGateOut(gr,s800,m3c,"all");
  for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
    if( (ftac==1 && InPartCut[in]->IsInside(tof->GetOBJ(),tof->GetXFP()))
	|| (ftac==3 && InPartCut[in]->IsInside(tof->GetOBJ(),tof->GetXFP())) 
	|| (ftac==0 && InPartCut[in]->IsInside(tof->GetTACOBJ(),tof->GetTACXFP())) 
	|| (ftac==2 && InPartCut[in]->IsInside(tof->GetTACOBJ(),tof->GetTACXFP())) 
	|| (ftac==5 && InPartCut[in]->IsInside(tof->GetMOBJ(),tof->GetMXFP()))
	|| (ftac==6 && InPartCut[in]->IsInside(tof->GetMOBJ(),tof->GetMXFP())) ){
      const char* inname = InPartCut[in]->GetName();
      FillHistogramsGateIn(gr,s800,m3c,inname);
      if(hasgret){
	FillHistogramsGateIn(gr,s800,m3c,Form("%s_coinc",inname));
      }
      FillHistogramsGateOut(gr,s800,m3c,inname);
      
      for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
	if((ftac == 1 && OutPartCut[in][ou]->IsInside(tof->GetOBJC(),ich->GetDE()))
	   ||(ftac == 0 && OutPartCut[in][ou]->IsInside(tof->GetTACOBJC(),ich->GetDE()))
	   ||(ftac == 3 && OutPartCut[in][ou]->IsInside(tof->GetXFPC(),ich->GetDE()))
	   ||(ftac == 2 && OutPartCut[in][ou]->IsInside(tof->GetTACXFPC(),ich->GetDE()))
	   ||(ftac == 5 && OutPartCut[in][ou]->IsInside(tof->GetMOBJC(),ich->GetDE()))
	   ||(ftac == 6 && OutPartCut[in][ou]->IsInside(tof->GetMXFPC(),ich->GetDE()))){
	  const char* outname = OutPartCut[in][ou]->GetName();

	  FillHistogramsGateOut(gr,s800,m3c,outname);
	  Fill(Form("txfp_%s",outname), 600,-300,300,track->GetXFP());
//	  if (hasgret){
//	    FillHistogramsGateOut(gr,s800,m3c,Form("%s_coinc",outname));
//	  }

	}

      }
    }
  }
  if(NoInCuts){
    // for pure beams, i.e. no incoming cut required or possible
    for(UShort_t ou=0;ou<OutPartCut[0].size();ou++){ // loop over PID cuts
      if((ftac == 1 && OutPartCut[0][ou]->IsInside(tof->GetOBJC(),ich->GetDE()))
	 ||(ftac == 0 && OutPartCut[0][ou]->IsInside(tof->GetTACOBJC(),ich->GetDE()))
	 ||(ftac == 3 && OutPartCut[0][ou]->IsInside(tof->GetXFPC(),ich->GetDE()))
	 ||(ftac == 2 && OutPartCut[0][ou]->IsInside(tof->GetTACXFPC(),ich->GetDE()))
	 ||(ftac == 5 && OutPartCut[0][ou]->IsInside(tof->GetMOBJC(),ich->GetDE()))
	 ||(ftac == 6 && OutPartCut[0][ou]->IsInside(tof->GetMXFPC(),ich->GetDE()))){
	const char* outname = OutPartCut[0][ou]->GetName();

	FillHistogramsGateOnlyOut(gr,s800,m3c,outname);

      }
    }
  }

}
void CalHistograms::FillHistogramsNoGate(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c){

  PAD* pad[2];
  PPAC* ppac[2];
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }
  IC* ich = s800->GetIC();
  STACKPIN* stackpin = s800->GetSTACKPIN(); //AR New in v4.1_TL
  MTDC* mtdc = s800->GetMTDC(); //AR New in v4.1_TL
  TOF* tof = s800->GetTOF();
  TRACK* track = s800->GetTRACK();
  IITRACK* iitrack = s800->GetIITRACK();
  HODO* hodo = s800->GetHODO();
  bool hasgret = gr->GetMult()>0;

  //Stack Pin
  for(UShort_t k=0 ; k<stackpin->GetCal().size() ; k++){
    Fill(Form("StackPinE_%d",k),fstackpin_range[0],fstackpin_range[1],fstackpin_range[2],stackpin->GetCal()[k]);//AR New in v4.1_TL
  }
  //Time diamond
  if(mtdc->Getch13().size()>0){
    Fill(Form("MTDC_%d",13),500,-500,0,mtdc->Getch13()[0]);//AR New in v4.1_TL
  }
  
  for(UShort_t p=0; p<2;p++){
    for(UShort_t i=0; i<ppac[p]->GetXStrips().size();i++){
      if(ppac[p]->GetXStrips()[i]>0)
	Fill(Form("ppac_%d_X",p),64,0,64,i,1000,0,100000,ppac[p]->GetXStrips()[i]);
    }
    for(UShort_t i=0; i<ppac[p]->GetYStrips().size();i++){
       if(ppac[p]->GetYStrips()[i]>0)
	 Fill(Form("ppac_%d_Y",p),64,0,64,i,1000,0,100000,ppac[p]->GetYStrips()[i]);
    }
  }
  Double_t sum =0;
  for(UShort_t c=0; c<hodo->GetEnergy()->size();c++){
    Short_t ch = hodo->GetChannel()->at(c);
    Fill(Form("hodo_%d",ch),1000,0,4000,hodo->GetEnergy()->at(c));
    Fill("hodo_vs_ch",32,0,32,ch,400,0,4000,hodo->GetEnergy()->at(c));
    Fill("hodo_all",1000,0,4000,hodo->GetEnergy()->at(c));
    Fill(Form("hodotime_%d",ch),4000,0,4000,hodo->GetTime());
    Fill("hodotime_vs_ch",32,0,32,ch,1000,0,4000,hodo->GetTime());
    Fill("hodotime_all",4000,0,4000,hodo->GetTime());
    if(hodo->GetEnergy()->at(c)>50)
      sum+=hodo->GetEnergy()->at(c);
  }// hodo energy
  Fill("hodo_sum",1000,0,4000,sum);

  if(fCal==0){
    Fill("obj_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fobj_range[0]/2,fobj_range[1],fobj_range[2],tof->GetOBJ());
    Fill("xfp_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fxfp_range[0]/2,fxfp_range[1],fxfp_range[2],tof->GetXFP());
    Fill("objtac_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 ftacobj_range[0]/2,ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ());
    Fill("xfptac_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 ftacxfp_range[0]/2,ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP());
    Fill("objm_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fmobj_range[0]/2,fmobj_range[1],fmobj_range[2],tof->GetMOBJ());
    Fill("xfpm_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fmxfp_range[0]/2,fmxfp_range[1],fmxfp_range[2],tof->GetMXFP());
    Fill("IC_vs_time",
	 fnentries/10000+1,0,(fnentries/10000+1)*10000,fentry,
	 fIC_range[0]/2,fIC_range[1],fIC_range[2],ich->GetDE());
  }

  

  // for(UShort_t p=0; p<2;p++){
  //   Fill(Form("ICde_vs_x_%d",p),
  // 	 600,-300,300,pad[p]->GetX(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  //   Fill(Form("ICde_vs_y_%d",p),
  // 	 200,-100,100,pad[p]->GetY(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
    // Fill(Form("ICsum_vs_x_%d",p),
    // 	 600,-300,300,pad[p]->GetX(),
    // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
    // Fill(Form("ICsum_vs_y_%d",p),
    // 	 200,-100,100,pad[p]->GetY(),
    // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // }

  if(hasgret){
    Fill("xfp_vs_obj_coinc",
	 fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
	 fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP());
  }


  Fill("xfp_vs_obj",
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP());
  Fill("x_vs_obj",
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_obj",
       fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
       100,-100,100,track->GetAFP());
  Fill("x_vs_objC",
       fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_objC",
       fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
       100,-50,50,track->GetAFP());

  Fill("x_vs_xfp",
       fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfp",
       fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
       100,-50,50,track->GetAFP());
  Fill("x_vs_xfpC",
       fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfpC",
       fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
       100,-50,50,track->GetAFP());

  Fill("xfp_vs_objtac",
       ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
       ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP());
  Fill("x_vs_objtac",
       ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_objtac",
       ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
       100,-50,50,track->GetAFP());
  Fill("x_vs_objtacC",
       ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_objtacC",
       ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
       100,-50,50,track->GetAFP());

  Fill("x_vs_xfptac",
       ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfptac",
       ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
       100,-50,50,track->GetAFP());
  Fill("x_vs_xfptacC",
       ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
       600,-300,300,pad[0]->GetX());
  Fill("afp_vs_xfptacC",
       ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
       100,-50,50,track->GetAFP());

  for(UShort_t i=0;i<tof->GetMOBJV()->size();i++){
    for(UShort_t j=0;j<tof->GetMXFPV()->size();j++){
      Fill("xfp_vs_objm",
	   fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJV()->at(i),
	   fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFPV()->at(j));
    }
  }
  for(UShort_t i=0;i<tof->GetMOBJCV()->size();i++){
    for(UShort_t j=0;j<tof->GetMXFPCV()->size();j++){
      Fill("xfpc_vs_objmc",
	   fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJCV()->at(i),
	   fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFPCV()->at(j));
    }
  }
  
  for(UShort_t i=0;i<tof->GetMOBJV()->size();i++){
    Fill("x_vs_objm",
	 fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJV()->at(i),
	 600,-300,300,pad[0]->GetX());
    Fill("afp_vs_objm",
	 fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJV()->at(i),
	 100,-50,50,track->GetAFP());
    Fill("x_vs_objmC",
	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJCV()->at(i),
	 600,-300,300,pad[0]->GetX());
    Fill("afp_vs_objmC",
	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJCV()->at(i),
	 100,-50,50,track->GetAFP());
  }
  for(UShort_t i=0;i<tof->GetMXFPV()->size();i++){
    Fill("x_vs_xfpm",
	 fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFPV()->at(i),
	 600,-300,300,pad[0]->GetX());
    Fill("afp_vs_xfpm",
	 fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFPV()->at(i),
	 100,-50,50,track->GetAFP());
    Fill("x_vs_xfpmC",
	 fmxfpC_range[0],fmxfpC_range[1],fmxfpC_range[2],tof->GetMXFPCV()->at(i),
	 600,-300,300,pad[0]->GetX());
    Fill("afp_vs_xfpmC",
	 fmxfpC_range[0],fmxfpC_range[1],fmxfpC_range[2],tof->GetMXFPCV()->at(i),
	 100,-50,50,track->GetAFP());
  }


  for (UShort_t g=0; g<gr->GetMult(); g++){
    HitCalc* hit = gr->GetHit(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    TVector3 pos = hit->GetPosition();//AR New in v4.1_TL
    Fill("egam",
	 8000,0,8000,energy);
    Fill(Form("egam_depth"),
	 8000,0,8000,energy,
	 2000,-50,150,hit->GetDepth());
    if(g==0 && gr->GetMult()==1){
      Fill(Form("egam_depth_MultOne"),
	   8000,0,8000,energy,
	   2000,-50,150,hit->GetDepth());
    }
    Fill("egam_tgam",
	 4000,0,4000,energy,
	 400,-200,200,hit->GetTime());
    Fill("egamdc",
	 8000,0,8000,energy_dc);
    Fill("egamdc_tgam",
	 4000,0,4000,energy_dc,
	 400,-200,200,hit->GetTime());
    Fill("egam_summary",
	 32,-0.5,31.5,4*fSett->Hole2Det(hit->GetHole())+hit->GetCrystal(),
	 4000,0,4000,energy);
   }
  for(UShort_t g=0;g<gr->GetMultAB();g++){
    HitCalc* hit = gr->GetHitAB(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    Fill("egamAB",
	 8000,0,8000,energy);
    Fill("egamAB_tgam",
	 4000,0,4000,energy,
	 400,-200,200,hit->GetTime());
    Fill("egamABdc",
	 8000,0,8000,energy_dc);
    Fill("egamABdc_tgam",
	 4000,0,4000,energy_dc,
	 400,-200,200,hit->GetTime());
    Fill("egamAB_summary",
	 32,-0.5,31.5,4*fSett->Hole2Det(hit->GetHole())+hit->GetCrystal(),
	 4000,0,4000,energy);
   }

}


void CalHistograms::FillHistogramsGateIn(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, const char* inname){
  PAD* pad[2];
  PPAC* ppac[2];
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }
  IC* ich = s800->GetIC();
  TOF* tof = s800->GetTOF();
  TRACK* track = s800->GetTRACK();
  IITRACK* iitrack = s800->GetIITRACK();
  HODO* hodo = s800->GetHODO();
  bool hasgret = gr->GetMult()>0;
  STACKPIN* stackpin = s800->GetSTACKPIN(); //AR New in v4.1_TL
  MTDC* mtdc = s800->GetMTDC(); //AR New in v4.1_TL

  Fill("hcrdcrawrange0xy",
       2400,-600,600, s800->GetPAD(0)->GetX(),
       4096,0,4096,s800->GetPAD(0)->GetY());
  Fill("hcrdcrawrange1xy",
       2400,-600,600, s800->GetPAD(1)->GetX(),
       4096,0,4096,s800->GetPAD(1)->GetY());
  Fill(Form("hcrdc0xy_%s",inname),
       600,-300,300,s800->GetPAD(0)->GetX(),
       600,-300,300,s800->GetPAD(0)->GetY());
  Fill(Form("hcrdc1xy_%s",inname),
       600,-300,300,s800->GetPAD(1)->GetX(),
       600,-300,300,s800->GetPAD(1)->GetY());

  //Stack Pin
  for(UShort_t k=0 ; k<stackpin->GetCal().size() ; k++){
    Fill(Form("StackPinE_%d_%s",k,inname),fstackpin_range[0],fstackpin_range[1],fstackpin_range[2],stackpin->GetCal()[k]);//AR New in v4.1_TL
  }
  //Time diamond
  if(mtdc->Getch13().size()>0){
    Fill(Form("MTDC_%d_%s",13,inname),500,-500,0,mtdc->Getch13()[0]);//AR New in v4.1_TL
  }
  
  // Fill(Form("ICde_vs_obj_%s",inname),
  //      fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // Fill(Form("ICde_vs_objc_%s",inname),
  //      fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // if (hasgret){
  //   Fill(Form("ICde_vs_objc_coinc_%s",inname),
  // 	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // }

  Double_t sum =0;
  for(UShort_t c=0; c<hodo->GetEnergy()->size();c++){
    if(hodo->GetEnergy()->at(c)>50)
      sum+=hodo->GetEnergy()->at(c);
  }// hodo energy
  Fill(Form("hodo_sum_%s",inname),1000,0,4000,sum);
  // if(sum>100){
  //   Fill(Form("ICde_vs_objtacc_hodo_%s",inname),
  // 	 ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  //   Fill(Form("ICde_vs_objc_hodo_%s",inname),
  // 	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  //   Fill(Form("ICde_vs_objmc_hodo_%s",inname),
  //      fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // }

  // Fill(Form("ICde_vs_objtacc_%s",inname),
  //      ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  Fill(Form("ICde_vs_objmc_%s",inname),
       fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
       fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // Fill(Form("ICde_vs_xfpc_%s",inname),
  //      fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // Fill(Form("ICde_vs_xfptacc_%s",inname),
  //      ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // Fill(Form("ICde_vs_xfpmc_%s",inname),
  //      fmxfpC_range[0],fmxfpC_range[1],fmxfpC_range[2],tof->GetMXFPC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());

  // Fill(Form("ICsum_vs_objc_%s",inname),
  //      fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // Fill(Form("ICsum_vs_objtacc_%s",inname),
  //      ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // Fill(Form("ICsum_vs_objmc_%s",inname),
  //      fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // Fill(Form("ICsum_vs_xfpc_%s",inname),
  //      fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // Fill(Form("ICsum_vs_xfptacc_%s",inname),
  //      ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // Fill(Form("ICsum_vs_xfpmc_%s",inname),
  //      fmxfpC_range[0],fmxfpC_range[1],fmxfpC_range[2],tof->GetMXFPC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());

  // for(UShort_t p=0; p<2;p++){
  //   Fill(Form("ICde_vs_x%d_%s",p,inname),
  //   	 600,-300,300,pad[p]->GetX(),
  //   	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  //   Fill(Form("ICde_vs_y%d_%s",p,inname),
  //   	 200,-100,100,pad[p]->GetY(),
  //   	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // Fill(Form("ICsum_vs_x%d_%s",p,inname),
  // 	 600,-300,300,pad[p]->GetX(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // Fill(Form("ICsum_vs_y%d_%s",p,inname),
  // 	 200,-100,100,pad[p]->GetY(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
  // }

  for(UShort_t i=0;i<tof->GetMOBJV()->size();i++){
    for(UShort_t j=0;j<tof->GetMXFPV()->size();j++){
      Fill(Form("xfp_vs_objm_%s",inname),
	   fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJV()->at(i),
	   fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFPV()->at(j));
    }
  }
}


void CalHistograms::FillHistogramsGateOut(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, const char* outname){
  TString oname = outname;

  PAD* pad[2];
  PPAC* ppac[2];
  for(UShort_t p=0; p<2;p++){
    pad[p] = s800->GetPAD(p);
    ppac[p] = s800->GetPPAC(p);
  }
  IC* ich = s800->GetIC();
  TOF* tof = s800->GetTOF();
  TRACK* track = s800->GetTRACK();
  IITRACK* iitrack = s800->GetIITRACK();
  HODO* hodo = s800->GetHODO();
  bool hasgret = gr->GetMult()>0;

  for(UShort_t c=0; c<hodo->GetEnergy()->size();c++){
    Short_t ch = hodo->GetChannel()->at(c);
    //Fill(Form("hodo_%d_%s",ch,oname.Data()),1000,0,4000,hodo->GetEnergy()->at(c));
    Fill(Form("hodo_vs_ch_%s",oname.Data()),32,0,32,ch,400,0,4000,hodo->GetEnergy()->at(c));
    Fill(Form("hodo_all_%s",oname.Data()),1000,0,4000,hodo->GetEnergy()->at(c));
    //Fill(Form("hodotime_%d_%s",ch,oname.Data()),4000,0,4000,hodo->GetTime());
    if(hodo->GetEnergy()->at(c)>100){
      Fill(Form("hodotime_vs_ch_%s",oname.Data()),32,0,32,ch,1000,0,4000,hodo->GetTime());
      Fill(Form("hodotime_vs_en_%s",oname.Data()),400,0,4000,hodo->GetEnergy()->at(c),1000,0,4000,hodo->GetTime());
      Fill(Form("hodotime_all_%s",oname.Data()),4000,0,4000,hodo->GetTime());
      for (UShort_t g=0; g<gr->GetMult(); g++){
	HitCalc* hit = gr->GetHit(g);
	float energy_dc = hit->GetDCEnergy();
	Fill(Form("egamdc_hodo_%s",oname.Data()),
	     1000,0,4000,energy_dc,1000,0,4000,hodo->GetEnergy()->at(c));
      }
      for(UShort_t g=0;g<gr->GetMultAB();g++){
	HitCalc* hit = gr->GetHitAB(g);
	float energy_dc = hit->GetDCEnergy();
	Fill(Form("egamABdc_hodo_%s",oname.Data()),
	     1000,0,4000,energy_dc,1000,0,4000,hodo->GetEnergy()->at(c));
      }
      // Fill(Form("hodo_ppar_%s",oname.Data()),
      // 	   1000,0,4000,hodo->GetEnergy()->at(c),
      // 	   fPP_range[0],fPP_range[1],fPP_range[2],track->GetPpar());
      // Fill(Form("hodo_pparc_%s",oname.Data()),
      // 	   1000,0,4000,hodo->GetEnergy()->at(c),
      // 	   fPP_range[0],fPP_range[1],fPP_range[2],track->GetPparC());
    }//good hodo energy
  }// hodo energy

  // Fill(Form("ICde_vs_objmc_%s",oname.Data()),
  //      fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
  //      fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());

  //registr gates
  // if(s800->GetRegistr()==1){
  //   Fill(Form("ICde_vs_objmc_%s_registr_1",oname.Data()),
  // 	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // }
  // if(s800->GetRegistr()==3){
  //   Fill(Form("ICde_vs_objmc_%s_registr_3",oname.Data()),
  // 	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // }
  // if(s800->GetRegistr()==1 || s800->GetRegistr()==3){
  //   Fill(Form("ICde_vs_objmc_%s_registr_1or3",oname.Data()),
  // 	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
  // 	 fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
  // }
  
  for (UShort_t g=0; g<gr->GetMult(); g++){
    HitCalc* hit = gr->GetHit(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    TVector3 pos = hit->GetPosition();//AR New in v4.1_TL
    Fill(Form("egam_%s",oname.Data()),
	 8000,0,8000,energy);
    Fill(Form("egam_scatter_%s",oname.Data()),
	 500,0,0.5,track->GetTheta(),
	 4000,0,4000,energy);
    Fill(Form("egam_tgam_%s",oname.Data()),
	 4000,0,4000,energy,
	 400,-200,200,hit->GetTime());
    Fill(Form("egamdc_tgam_%s",oname.Data()),
	 4000,0,4000,energy_dc,
	 400,-200,200,hit->GetTime());
    if(hit->GetDepth()>=15. && hit->GetDepth()<120.){
      Fill(Form("egamdc_tgam_%s_dcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   400,-200,200,hit->GetTime());
    }
    Fill(Form("tgam_egam_%s",oname.Data()),
	 400,-200,200,hit->GetTime(),
	 4000,0,4000,energy);
    Fill(Form("egam_summary_%s",oname.Data()),
	 32,-0.5,31.5,4*fSett->Hole2Det(hit->GetHole())+hit->GetCrystal(),
	 4000,0,4000,energy);
    Fill(Form("egamdc_%s",oname.Data()),
	 4000,0,4000,energy_dc);
    Fill(Form("egamdc_scatter_%s",oname.Data()),
	 500,0,0.5,track->GetTheta(),
	 4000,0,4000,energy_dc);
    Fill(Form("egamdc_grettheta_%s",oname.Data()),
	 4000,0,4000,energy_dc,
	 180,0,3.14159,hit->GetPosition().Theta());
    // if(hit->GetDepth()>=-10. && hit->GetDepth()<20.){
    //   Fill(Form("egamdc_grettheta_Dcut1_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   180,0,3.14159,hit->GetPosition().Theta());
    // }
    // if(hit->GetDepth()>=20. && hit->GetDepth()<40.){
    //   Fill(Form("egamdc_grettheta_Dcut2_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   180,0,3.14159,hit->GetPosition().Theta());
    // }
    // if(hit->GetDepth()>=40. && hit->GetDepth()<60.){
    //   Fill(Form("egamdc_grettheta_Dcut3_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   180,0,3.14159,hit->GetPosition().Theta());
    // }
    // if(hit->GetDepth()>=60. && hit->GetDepth()<80.){
    //   Fill(Form("egamdc_grettheta_Dcut4_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   180,0,3.14159,hit->GetPosition().Theta());
    // }
    // if(hit->GetDepth()>=80. && hit->GetDepth()<100.){
    //   Fill(Form("egamdc_grettheta_Dcut5_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   180,0,3.14159,hit->GetPosition().Theta());
    // }
    // if(hit->GetDepth()>=100. && hit->GetDepth()<120.){
    //   Fill(Form("egamdc_grettheta_Dcut6_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   180,0,3.14159,hit->GetPosition().Theta());
    // }
    // Fill(Form("egamdc_ppar_%s",oname.Data()),
    // 	 4000,0,4000,energy_dc,
    // 	 fPP_range[0],fPP_range[1],fPP_range[2],track->GetPpar());
    // Fill(Form("egamdc_pparc_%s",oname.Data()),
    // 	 4000,0,4000,energy_dc,
    // 	 fPP_range[0],fPP_range[1],fPP_range[2],track->GetPparC());
    Fill(Form("egamdc_depth_%s",oname.Data()),
	 4000,0,4000,energy_dc,
	 2000,-50,150,hit->GetDepth());
    if(g==0 && gr->GetMult()==1){ //AR New in v4.1_TL
      Fill(Form("egamdc_tgam_MultOne_%s",oname.Data()),
	   4000,0,4000,energy_dc,
	   400,-200,200,hit->GetTime());
      Fill(Form("egamdc_MultOne_%s",oname.Data()),
	   4000,0,4000,energy_dc);	
      Fill(Form("egamdc_depth_MultOne_%s",oname.Data()),
	   4000,0,4000,energy_dc,
	   2000,-50,150,hit->GetDepth());
      if(hit->GetDepth()>=20. && hit->GetDepth()<80.){
	Fill(Form("egamdc_tgam_MultOne_%s_dcut",oname.Data()),
	     4000,0,4000,energy_dc,
	     400,-200,200,hit->GetTime());
      }
    }
    //if((g==0 || g==1) && gr->GetMult()==2){ //aldric 1/9/2020 Multiplicity 2 only
    //   Fill(Form("egamdc_tgam_MultTwo_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   400,-200,200,hit->GetTime());
    //   Fill(Form("egamdc_MultTwo_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc);	
    // }
    // if((g==0 || g==1 || g==2) && gr->GetMult()==3){ //aldric 1/9/2020 Multiplicity 3 only
    //   Fill(Form("egamdc_tgam_MultThree_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc,
    // 	   400,-200,200,hit->GetTime());
    //   Fill(Form("egamdc_MultThree_%s",oname.Data()),
    // 	   4000,0,4000,energy_dc);	
    // }
    if(foundTCut && TimeCut->IsInside(hit->GetTime(),energy)){
      Fill(Form("egam_%s_tcut",oname.Data()),
	   4000,0,4000,energy);
      Fill(Form("egam_scatter_%s_tcut",oname.Data()),
	   500,0,0.5,track->GetTheta(),
	   4000,0,4000,energy);
      Fill(Form("egam_summary_%s_tcut",oname.Data()),
	   28,-0.5,27.5,4*fSett->Hole2Det(hit->GetHole())+hit->GetCrystal(),
	   4000,0,4000,energy);
      Fill(Form("egamdc_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc);
      Fill(Form("egamdc_scatter_%s_tcut",oname.Data()),
	   500,0,0.5,track->GetTheta(),
	   4000,0,4000,energy_dc);
      if(hit->GetDepth()>=15. && hit->GetDepth()<120.){
	Fill(Form("egamdc_grettheta_%s_tcut_dcut",oname.Data()),
	     4000,0,4000,energy_dc,
	     180,0,3.14159,hit->GetPosition().Theta());
      }
      Fill(Form("egamdc_grettheta_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   180,0,3.14159,hit->GetPosition().Theta());
      // Fill(Form("egamdc_ppar_%s_tcut",oname.Data()),
      // 	   4000,0,4000,energy_dc,
      // 	   fPP_range[0],fPP_range[1],fPP_range[2],track->GetPpar());
      // if(hit->GetDepth()>=-10. && hit->GetDepth()<20.){
      // 	Fill(Form("egamdc_grettheta_Dcut1_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc,
      // 	     180,0,3.14159,hit->GetPosition().Theta());
      // }
      // if(hit->GetDepth()>=20. && hit->GetDepth()<40.){
      // 	Fill(Form("egamdc_grettheta_Dcut2_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc,
      // 	     180,0,3.14159,hit->GetPosition().Theta());
      // }
      // if(hit->GetDepth()>=40. && hit->GetDepth()<60.){
      // 	Fill(Form("egamdc_grettheta_Dcut3_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc,
      // 	     180,0,3.14159,hit->GetPosition().Theta());
      // }
      // if(hit->GetDepth()>=60. && hit->GetDepth()<80.){
      // 	Fill(Form("egamdc_grettheta_Dcut4_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc,
      // 	     180,0,3.14159,hit->GetPosition().Theta());
      // }
      // if(hit->GetDepth()>=80. && hit->GetDepth()<100.){
      // 	Fill(Form("egamdc_grettheta_Dcut5_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc,
      // 	     180,0,3.14159,hit->GetPosition().Theta());
      // }
      // if(hit->GetDepth()>=100. && hit->GetDepth()<120.){
      // 	Fill(Form("egamdc_grettheta_Dcut6_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc,
      // 	     180,0,3.14159,hit->GetPosition().Theta());
      // }
      Fill(Form("egamdc_depth_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   2000,-50,150,hit->GetDepth());
      if(g==0 && gr->GetMult()==1){ //AR New in v4.1_TL -> Multiplicity 1 only
	// if(hit->GetDepth()>=-10. && hit->GetDepth()<20.){
	//   Fill(Form("egamdc_grettheta_Dcut1_MultOne_%s_tcut",oname.Data()),
	//        4000,0,4000,energy_dc,
	//        180,0,3.14159,hit->GetPosition().Theta());
	// }
	// if(hit->GetDepth()>=20. && hit->GetDepth()<40.){
	//   Fill(Form("egamdc_grettheta_Dcut2_MultOne_%s_tcut",oname.Data()),
	//        4000,0,4000,energy_dc,
	//        180,0,3.14159,hit->GetPosition().Theta());
	// }
	// if(hit->GetDepth()>=40. && hit->GetDepth()<60.){
	//   Fill(Form("egamdc_grettheta_Dcut3_MultOne_%s_tcut",oname.Data()),
	//        4000,0,4000,energy_dc,
	//        180,0,3.14159,hit->GetPosition().Theta());
	// }
	// if(hit->GetDepth()>=60. && hit->GetDepth()<80.){
	//   Fill(Form("egamdc_grettheta_Dcut4_MultOne_%s_tcut",oname.Data()),
	//        4000,0,4000,energy_dc,
	//        180,0,3.14159,hit->GetPosition().Theta());
	// }
	// if(hit->GetDepth()>=80. && hit->GetDepth()<100.){
	//   Fill(Form("egamdc_grettheta_Dcut5_MultOne_%s_tcut",oname.Data()),
	//        4000,0,4000,energy_dc,
	//        180,0,3.14159,hit->GetPosition().Theta());
	// }
	// if(hit->GetDepth()>=100. && hit->GetDepth()<120.){
	//   Fill(Form("egamdc_grettheta_Dcut6_MultOne_%s_tcut",oname.Data()),
	//        4000,0,4000,energy_dc,
	//        180,0,3.14159,hit->GetPosition().Theta());
	// }
	Fill(Form("egamdc_MultOne_%s_tcut",oname.Data()),
	     4000,0,4000,energy_dc);
	Fill(Form("egamdc_depth_MultOne_%s_tcut",oname.Data()),
	     4000,0,4000,energy_dc,
	     2000,-50,150,hit->GetDepth());
	if(hit->GetDepth()>=15. && hit->GetDepth()<120.){
	  Fill(Form("egamdc_grettheta_MultOne_%s_tcut_dcut",oname.Data()),
	       4000,0,4000,energy_dc,
	       180,0,3.14159,hit->GetPosition().Theta());
	}
	Fill(Form("egamdc_grettheta_MultOne_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   180,0,3.14159,hit->GetPosition().Theta());
      }
      // if((g==0 || g==1) && gr->GetMult()==2){ //aldric 1/9/2020 Multiplicity 2 only
      // 	Fill(Form("egamdc_MultTwo_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc);
      // 	Fill(Form("egamdc_grettheta_MultTwo_%s_tcut",oname.Data()),
      // 	   4000,0,4000,energy_dc,
      // 	   180,0,3.14159,hit->GetPosition().Theta());
      // }
      // if((g==0 || g==1 || g==2)  && gr->GetMult()==3){ //aldric 1/9/2020 Multiplicity 3 only
      // 	Fill(Form("egamdc_MultThree_%s_tcut",oname.Data()),
      // 	     4000,0,4000,energy_dc);
      // 	Fill(Form("egamdc_grettheta_MultThree_%s_tcut",oname.Data()),
      // 	   4000,0,4000,energy_dc,
      // 	   180,0,3.14159,hit->GetPosition().Theta());
      // }
    }
    int q = fSett->Hole2Det(hit->GetHole())+1;
    if((q==1 || q==7) && hit->GetCrystal()==1){
      Fill(Form("egamdc_BW_%s",oname.Data()),
	   4000,0,4000,hit->GetDCEnergy(0.3));
    }
    else if((q==2 || q==6 || q==4 || q==5 ) && hit->GetCrystal()==2){
      Fill(Form("egamdc_FW_%s",oname.Data()),
	   4000,0,4000,hit->GetDCEnergy(0.3));
    }
    Fill(Form("egamdc_s800phi_%s",oname.Data()),
	 4000,0,4000,energy_dc,
	 180,-2*3.14159,0,s800->GetTRACK()->GetPhi());
    Fill(Form("egamdc_gretphi_%s",oname.Data()),
	 4000,0,4000,energy_dc,
	 180,-3.14159,3.14159,hit->GetPosition().Phi());
  }

  for(UShort_t g=0;g<gr->GetMultAB();g++){
    HitCalc* hit = gr->GetHitAB(g);
    float energy = hit->GetEnergy();
    float energy_dc = hit->GetDCEnergy();
    Fill(Form("egamAB_%s",oname.Data()),
	 8000,0,8000,energy);
    Fill(Form("egamABdc_%s",oname.Data()),
	 8000,0,8000,energy_dc);
    Fill(Form("egamABdc_tgam_%s",oname.Data()),
	 8000,0,8000,energy_dc,
	 400,-200,200,hit->GetTime());
    Fill(Form("egamABdc_grettheta_%s",oname.Data()),
	 4000,0,4000,energy_dc,
	 180,0,3.14159,hit->GetPosition().Theta());
    // Fill(Form("egamABdc_ppar_%s",oname.Data()),
    // 	 4000,0,4000,energy_dc,
    // 	 fPP_range[0],fPP_range[1],fPP_range[2],track->GetPpar());
    // Fill(Form("egamABdc_pparc_%s",oname.Data()),
    // 	 4000,0,4000,energy_dc,
    // 	 fPP_range[0],fPP_range[1],fPP_range[2],track->GetPparC());
    if(foundTCut && TimeCut->IsInside(hit->GetTime(),energy)){
      Fill(Form("egamAB_%s_tcut",oname.Data()),
	   4000,0,4000,energy);
      Fill(Form("egamABdc_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc);
      Fill(Form("egamABdc_grettheta_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   180,0,3.14159,hit->GetPosition().Theta());
      Fill(Form("egamABdc_depth_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   2000,-50,150,hit->GetDepth());
      Fill(Form("egamABdc_span_%s_tcut",oname.Data()),
	   4000,0,4000,energy_dc,
	   2000,-50,350,hit->GetSpan());
      if(hit->GetSpan()<20.){
	Fill(Form("egamABdc_%s_ABVeto_tcut",oname.Data()),
	     4000,0,4000,energy_dc);
	Fill(Form("egamABdc_grettheta_%s_ABVeto_tcut",oname.Data()),
	     4000,0,4000,energy_dc,
	     180,0,3.14159,hit->GetPosition().Theta());
      }
      if(g==0 && gr->GetMultAB()==1){
	Fill(Form("egamABdc_grettheta_MultOne_%s_tcut",oname.Data()),
	     4000,0,4000,energy_dc,
	     180,0,3.14159,hit->GetPosition().Theta());
	Fill(Form("egamABdc_depth_MultOne_%s_tcut",oname.Data()),
	     4000,0,4000,energy_dc,
	     2000,-50,150,hit->GetDepth());
      }
      // Fill(Form("egamABdc_ppar_%s_tcut",oname.Data()),
      // 	   4000,0,4000,energy_dc,
      // 	   fPP_range[0],fPP_range[1],fPP_range[2],track->GetPpar());
   }
  }

  // Fill(Form("IC_vs_trackxfp_%s",oname.Data()),
  //      600,-300,300,track->GetXFP(),1000,0,1000,ich->GetDE());

  Fill(Form("ata_%s",oname.Data()),
       200,-100,100,track->GetATA());
  Fill(Form("bta_%s",oname.Data()),
       200,-100,100,track->GetBTA());
  Fill(Form("ata_bta_%s",oname.Data()),
       200,-100,100,track->GetATA(),
       200,-100,100,track->GetBTA());
  Fill(Form("scatter_%s",oname.Data()),
       500,0,0.5,track->GetTheta());

  // Various gamma-gamma spectra with no additional gates.
  {
    int highestg=-1;
    double highesten =0;
    for(UShort_t g=0;g<gr->GetMult();g++){
      if(gr->GetHit(g)->GetEnergy()<1 || gr->GetHit(g)->GetEnergy()>6000)
	continue;
      if(gr->GetHit(g)->GetDCEnergy()>highesten){
	highesten = gr->GetHit(g)->GetDCEnergy();
	highestg =g;
      }
    }
    for(UShort_t g=0;g<gr->GetMult();g++){
      if(gr->GetHit(g)->GetEnergy()<1 || gr->GetHit(g)->GetEnergy()>6000)
	continue;
      if(highestg>-1 && g!=highestg){
	Fill(Form("egamegamdc_fold_%s",oname.Data()),
	     4000,0,4000,highesten,
	     4000,0,4000,gr->GetHit(g)->GetDCEnergy());
      }
    }
    highestg=-1;
    highesten =0;
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      if(gr->GetHitAB(g)->GetEnergy()<1 || gr->GetHitAB(g)->GetEnergy()>6000)
	continue;
      if(gr->GetHitAB(g)->GetDCEnergy()>highesten){
	highesten = gr->GetHitAB(g)->GetDCEnergy();
	highestg =g;
      }
    }
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      if(gr->GetHitAB(g)->GetEnergy()<1 || gr->GetHitAB(g)->GetEnergy()>6000)
	continue;
      if(highestg>-1 && g!=highestg){
	Fill(Form("egamegamABdc_fold_%s",oname.Data()),
	     4000,0,4000,highesten,
	     4000,0,4000,gr->GetHitAB(g)->GetDCEnergy());
      }
    }
  }
  for(UShort_t g1=0;g1<gr->GetMult();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMult();g2++){
      Fill(Form("egamegamdc_%s",oname.Data()),
	   1000,0,4000,gr->GetHit(g1)->GetDCEnergy(),
	   1000,0,4000,gr->GetHit(g2)->GetDCEnergy());
      Fill(Form("egamegamdc_sym_%s",oname.Data()),
	   1000,0,4000,gr->GetHit(g1)->GetDCEnergy(),
	   1000,0,4000,gr->GetHit(g2)->GetDCEnergy());
      Fill(Form("egamegamdc_sym_%s",oname.Data()),
	   1000,0,4000,gr->GetHit(g2)->GetDCEnergy(),
	   1000,0,4000,gr->GetHit(g1)->GetDCEnergy());
    }
  }
  
  //AR New in v4.1.LT -> Necessary Histograms for Rob's Cascade Analysis
  if(fUpstream>0){
    
    Fill(Form("mult_global_%s",oname.Data()),       
	 20,0,20,gr->GetMult());
    
    Fill(Form("multAB_global_%s",oname.Data()),       
	 20,0,20,gr->GetMultAB());
    
    
    True_GammaKnown = fSett->GetTrueE_gk(); // in keV
    double True_GoI = fSett->GetTrueE_goi(); // in keV 
    beta_r = fSett->TargetBeta();
    
    //RE New in v4.1_LT -> Use the inverse map dta info to correct beta.
    RelGamma = 1/sqrt(1-beta_r*beta_r);
 
    //RE New in v4.1_LT ->  checking effect of dta correction
    if(!isnan(track->GetDTA())){                                    //ie, dta is valid, use it to correct beta
      RelGamma = RelGamma + (RelGamma - 1.0)*track->GetDTA()/100.;
      beta_r = sqrt(1.0-1.0/(RelGamma*RelGamma));
    }
    else{                                                           //dta is "not a number", just use the central beta.  
      ;                                                             //do nothing.  perhaps these events should just 
    } 

    Z_Origin_Tar = fSett->GetZ_Origin_Tar();  //in mm
    Z_Origin_FocalPoint = -721.265;  //in mm.  Not currently used.  
    
    double LabFrameLoE_gk = fSett->GetLabFrameLowE_gk();
    double LabFrameHiE_gk = fSett->GetLabFrameHighE_gk();
    double LabFrameLoE_goi = fSett->GetLabFrameLowE_goi();
    double LabFrameHiE_goi = fSett->GetLabFrameHighE_goi();
    double comptThetaThresh = 0.5;

    //RE New in v4.1_LT ->  I want to keep track of which gamma's are consistent with the "known gamma"
    vector<int> listGoodgk;
    int chosenGoodgk = -1;   //This refers to that gamma which is consistent with the "known gamma" that I choose to calculate the decay position from.  default value -1 signifies that NO gammas are consistent with "known gamma"
    for(int i=0; i<gr->GetMultAB(); i++){
      int thisEnergy = gr->GetHitAB(i)->GetEnergy();
      if(thisEnergy > LabFrameLoE_gk && thisEnergy < LabFrameHiE_gk){
	listGoodgk.push_back(i);
      }
      
    }

    
    // We do the Cascade Doppler-correction assuming a thin pencil beam.  (ie V1.0)  No addback.
    if(gr->GetMult()>=2){       
      for(UShort_t g1=0;g1<gr->GetMult();g1++){
	if((gr->GetHit(g1)->GetEnergy()>LabFrameLoE_gk && gr->GetHit(g1)->GetEnergy()<LabFrameHiE_gk)){
	  Shifted_GammaKnown = gr->GetHit(g1)->GetEnergy(); 
	  TVector3 decay = CascadeZdecay(gr->GetHit(g1),True_GammaKnown,beta_r);
	  CascadeZ = decay.Z();
	  
	  Rho_Tar_GammaKnown = rho(gr->GetHit(g1)->GetPosition());

	  //RE New in v4.1_LT ->  Calculate the difference between compton and geometric angles for THIS ordering.
	  double gkSegAngDiff = 10;
	  double gkIPAngDiff = 10;

	  if(gr->GetHit(g1)->GetSegMult()>1){
	    gkSegAngDiff = abs(CasZangle(gr->GetHit(g1)->GetSegPos(0),
					 gr->GetHit(g1)->GetSegPos(1),
					 decay)
			       -acos(CosComptAngle(gr->GetHit(g1)->GetSegEn(0),
						   gr->GetHit(g1)->GetEnergy()-gr->GetHit(g1)->GetSegEn(0))));
	  }
	  else{
	    gkSegAngDiff=0;     //no possible compton scatter to campare this with.
	  }
	  
	  Fill(Form("CascadeZold_%s",oname.Data()),
	       400,0,2000,CascadeZ);
	  
	  Fill(Form("CascadeZold_thisegam_%s",oname.Data()),
	       2000,0,8000,Shifted_GammaKnown,
	       400,0,2000,CascadeZ);
	  Fill(Form("CascadeZ_thisegamdcOld_%s",oname.Data()),
	       2000,0,8000,gr->GetHit(g1)->GetDCEnergy(),
	       400,0,2000,CascadeZ);

	  for(UShort_t g2=0;g2<gr->GetMult();g2++){
	    if(g2!=g1){
	      Shifted_GammaOfInterest = gr->GetHit(g2)->GetEnergy(); 
	      CorrectedE = CascadeCorrectedE(gr->GetHit(g2),decay,beta_r);
	      
	      Rho_Tar_GammaOfInterest = rho(gr->GetHit(g2)->GetPosition());
	      
	      if(gr->GetHit(g1)->GetSegMult()>1){                                  
		
		Fill(Form("gkSegAngDiff_vs_CorrectedEold_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     1000,0,3.2,gkSegAngDiff);
		
		if(gkSegAngDiff<0.5){
		  Fill(Form("CascadeZ_vs_CorrectedEold_gkSegAngDiff0p5_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		  
		  if(gr->GetHit(g2)->GetSegMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedEold_gkSegAngDiff0p5_GoIseg1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    
		    if(gr->GetHit(g1)->GetSegMult()==2){
		      Fill(Form("CascadeZ_vs_CorrectedEold_gkSegAngDiff0p5_GoIseg1_GKseg2_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		    if(gr->GetHit(g1)->GetSegMult()==3){
		      Fill(Form("CascadeZ_vs_CorrectedEold_gkSegAngDiff0p5_GoIseg1_GKseg3_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }
		}
	      }
	      
	      Fill(Form("CorrectedEold_%s",oname.Data()),
		   2000,0,4000,CorrectedE);
	      
	      //RE New in v4.1_LT -> Making spectrum to see how making a timing cut may affect peak/bgd.
	      Fill(Form("tgam_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   400,-200,200,gr->GetHit(g1)->GetTime());
	      
	      //RE New in v4.1_LT -> Making spectrum to see how making gating on events within a certain cylindrical radius
	      //of the "low energy" hit may affect peak/bgd.
	      Fill(Form("RhoGoI_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   500,0,500,Rho_Tar_GammaOfInterest);
	      //RE New in v4.1_LT -> sim for the "high energy" gamma
	      Fill(Form("RhoGK_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   500,0,500,Rho_Tar_GammaKnown);

	      Fill(Form("CascadeZ_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   600,-1000,2000,CascadeZ);
	      
	      Fill(Form("thisegam_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   500,0,2000,Shifted_GammaKnown);
	      
	      Fill(Form("otheregam_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   500,0,2000,Shifted_GammaOfInterest);
	      
	      Fill(Form("thisegamdc_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   1000,0,4000,gr->GetHit(g1)->GetDCEnergy());
	      
	      Fill(Form("othergamdc_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   1000,0,4000,gr->GetHit(g2)->GetDCEnergy());

	      Fill(Form("mult_vs_CorrectedEold_%s",oname.Data()),   
		   1000,0,4000,CorrectedE,                          
		   10,0,10,gr->GetMult());
	     
	      Fill(Form("GoISegmult_vs_CorrectedEold_%s",oname.Data()),
		   1000,0,4000,CorrectedE,
		   20,0,20,gr->GetHit(g2)->GetSegMult());
	      
	      Fill(Form("GKSegmult_vs_CorrectedEold_%s",oname.Data()), 
		   1000,0,4000,CorrectedE,
		   20,0,20,gr->GetHit(g1)->GetSegMult());
	      
	      Fill(Form("GoIipmult_vs_CorrectedEold_%s",oname.Data()), 
		   1000,0,4000,CorrectedE,
		   20,0,20,gr->GetHit(g2)->GetIPMult());
	      
	      Fill(Form("GKipmult_vs_CorrectedEold_%s",oname.Data()),  
		   1000,0,4000,CorrectedE,
		   20,0,20,gr->GetHit(g1)->GetIPMult());
	      
	      Fill(Form("RhoGoI_vs_ZGoI_%s",oname.Data()),       
		   1000,0,2000,Z_Tar_GammaOfInterest,            
		   500,0,500,Rho_Tar_GammaOfInterest);           
	      
	      if(gr->GetHit(g2)->GetIPMult()==1){
		Fill(Form("RhoGoI_vs_CorrectedEold_GoIipmult1_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaOfInterest);
	      }
	      
	      else if(gr->GetHit(g2)->GetIPMult()==2){
		Fill(Form("RhoGoI_vs_CorrectedEold_GoIipmult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaOfInterest);
	      }
	      
	      else if(gr->GetHit(g2)->GetIPMult()==3){
		Fill(Form("RhoGoI_vs_CorrectedEold_GoIipmult3_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaOfInterest);
	      }
	      
	      
	      if(gr->GetHit(g2)->GetIPMult()==1){
		Fill(Form("GKipmult_vs_CorrectedEold_GoIipmult1_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g1)->GetIPMult());
	      }
	      else if(gr->GetHit(g2)->GetIPMult()==2){
		Fill(Form("GKipmult_vs_CorrectedEold_GoIipmult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g1)->GetIPMult());
	      }
	      
	      else if(gr->GetHit(g2)->GetIPMult()==3){
		Fill(Form("GKipmult_vs_CorrectedEold_GoIipmult3_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g1)->GetIPMult());
	      }
	      
	      //RE New in v4.1_LT ->  Want to see cases where IPMult=1 for both gammas (g1 and g2).
	      //In these cases, our choice of the "first interaction point" does not matter.
	      
	      if(gr->GetHit(g1)->GetIPMult()==1){
		if(gr->GetHit(g2)->GetIPMult()==1){
		  Fill(Form("CascadeZ_vs_CorrectedEold_goiIPMult1_gkIPMult1_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
	      }
	      if(gr->GetHit(g1)->GetSegMult()==1){
		if(gr->GetHit(g2)->GetSegMult()==1){
		  Fill(Form("CascadeZ_vs_CorrectedEold_goiSegMult1_gkSegMult1_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
	      }
	      
	      //RE New in v4.1_LT ->  Adding more histograms to look at energy through potentially useful gates.  
	      if(gr->GetMult()==2){


		//RE New in v4.1_LT ->  Want to see cases where IPMult=1 for both gammas (g1 and g2).
		//In these cases, our choice of the "first interaction point" does not matter.
		if(gr->GetHit(g1)->GetIPMult()==1){
		  if(gr->GetHit(g2)->GetIPMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_goiIPMult1_gkIPMult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		}
		if(gr->GetHit(g1)->GetSegMult()==1){
		  if(gr->GetHit(g2)->GetSegMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_goiSegMult1_gkSegMult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		}
		

		//RE New in v4.1_LT ->  looking at rho gate without any segment multiplicity gate.
		//RE New in v4.1_LT -> Adding gates on my new "span" quantity.  
		if(Rho_Tar_GammaOfInterest<230.){
		  Fill(Form("CascadeZ_vs_CorrectedEold_mult2_goiRho230_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		  Fill(Form("CorrectedEold_mult2_goiRho230_%s",oname.Data()),
		       2000,0,4000,CorrectedE);
		}
		
		Fill(Form("dta_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     200,-10,10,track->GetDTA());

		if(CascadeZ < 500.0){             //RE New in v4.1_LT ->  With this gate I can look at only forward emitted
		                                  //gammas and judge if CorrectedE has dependence on uncorrected beta.  
		  Fill(Form("dta_vs_CorrectedEold_mult2_nearDecays_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       200,-10,10,track->GetDTA());
		}

		Fill(Form("gkEn_vs_CorrectedEold_mult2_%s",oname.Data()),  //RE New in v4.1_LT ->  Adding histos to check
		     1000,0,4000,CorrectedE,                              //whether CorrectedE has problematic dependencies.  
		     1000,0,2000,gr->GetHit(g1)->GetEnergy());

		Fill(Form("goiEn_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     1000,0,2000,gr->GetHit(g2)->GetEnergy());

		Fill(Form("gkZ_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,1000,gr->GetHit(g1)->GetPosition().Z());

		Fill(Form("gkRho_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaKnown);

		Fill(Form("goiZ_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,1000,gr->GetHit(g2)->GetPosition().Z());
		
		Fill(Form("goiRho_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaOfInterest);
		
		if(gr->GetHit(g1)->GetSegMult()>1){                                  
		  
		  Fill(Form("gkSegAngDiff_vs_CorrectedEold_mult2_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       1000,0,3.2,gkSegAngDiff);
		  
		  if(gkSegAngDiff<0.5){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_gkSegAngDiff0p5_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    
		    if(gr->GetHit(g2)->GetSegMult()==1){
		      Fill(Form("CascadeZ_vs_CorrectedEold_mult2_gkSegAngDiff0p5_GoIseg1_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		      
		      if(gr->GetHit(g1)->GetSegMult()==2){
			Fill(Form("CascadeZ_vs_CorrectedEold_mult2_gkSegAngDiff0p5_GoIseg1_GKseg2_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		      if(gr->GetHit(g1)->GetSegMult()==3){
			Fill(Form("CascadeZ_vs_CorrectedEold_mult2_gkSegAngDiff0p5_GoIseg1_GKseg3_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		  }
		}
		
		
		Fill(Form("CorrectedEold_mult2_%s",oname.Data()),
		     2000,0,4000,CorrectedE);

		//RE New in v4.1_LT -> Changing the range of CascadeZ axis so that the 30P case (target is at GRETINA center can be plotted usefully.  Change back for 32Mg Upstream case.  
		Fill(Form("CascadeZ_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     600,-1000,2000,CascadeZ);
		
		Fill(Form("GoISegmult_vs_CorrectedEold_mult2_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g2)->GetSegMult());
		
		Fill(Form("GKSegmult_vs_CorrectedEold_mult2_%s",oname.Data()),    
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g1)->GetSegMult());
		
		
		Fill(Form("RhoGoI_vs_CorrectedEold_mult2_%s",oname.Data()),      
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaOfInterest);
		
		
		Fill(Form("GoIipmult_vs_CorrectedEold_mult2_%s",oname.Data()),  
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g2)->GetIPMult());
		
		
		Fill(Form("GKipmult_vs_CorrectedEold_mult2_%s",oname.Data()),   
		     1000,0,4000,CorrectedE,
		     20,0,20,gr->GetHit(g1)->GetIPMult());
		
		
		if(Rho_Tar_GammaKnown<=200.){
		  Fill(Form("RhoGoI_vs_CorrectedEold_mult2_rhoGK200_%s",oname.Data()),    
		       1000,0,4000,CorrectedE,
		       500,0,500,Rho_Tar_GammaOfInterest);
		}
		
		
		Fill(Form("RhoGK_vs_CorrectedEold_mult2_%s",oname.Data()),      
		     1000,0,4000,CorrectedE,
		     500,0,500,Rho_Tar_GammaKnown);
		
		
		if(Rho_Tar_GammaOfInterest<=180.){
		  Fill(Form("RhoGK_vs_CorrectedEold_mult2_rhoGoI180_%s",oname.Data()),   
		       1000,0,4000,CorrectedE,
		       500,0,500,Rho_Tar_GammaKnown);
		}
		
		//"low energy" seg multiplicity = 1 (since the physical events we care about should mostly do photoelec.
		if(gr->GetHit(g2)->GetSegMult()==1){
		  Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoIsegmult1_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		  
		  Fill(Form("GKSegmult_vs_CorrectedEold_mult2_GoIsegmult1_%s",oname.Data()),        
		       1000,0,4000,CorrectedE,
		       20,0,20,gr->GetHit(g1)->GetSegMult());
		  
		  //further, "high energy" seg multiplicity = 2,3,4 (since physically it should compton scatter some
		  if(gr->GetHit(g1)->GetSegMult()==2){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);

		    Fill(Form("dta_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 200,-10,10,track->GetDTA());
		    
		    //RE New in v4.1_LT -> Making spectrum to see how making a timing cut may affect peak/bgd.
		    Fill(Form("tgam_vs_CorrectedEold_mult2_GoIsegmult1_GKsegmult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,-200,200,gr->GetHit(g1)->GetTime());
		    
		    //RE New in v4.1_LT ->  Making spectrum to see how making gating on events within a certain cylindrical radius
		    //of the "low energy" hit may affect peak/bgd.
		    Fill(Form("RhoGoI_vs_CorrectedEold_mult2_GoIsegmult1_GKsegmult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaOfInterest);
		    //RE New in v4.1_LT ->  sim for the "high energy" gamma
		    Fill(Form("RhoGK_vs_CorrectedEold_mult2_GoIsegmult1_GKsegmult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaKnown);
		    
		    //RE New in v4.1_LT ->  Making a very selective spectra based on the useful ranges of the above specra.
		    if(Rho_Tar_GammaOfInterest<180.0 && Rho_Tar_GammaKnown<200.0){
		      Fill(Form("CorrectedEold_mult2_GoIsegmult1_GKsegmult2_rhoGoI180_rhoGK200_%s",oname.Data()),
			   2000,0,2000,CorrectedE);
		      Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_rhoGoI180_rhoGK200_%s",oname.Data()),
			   2000,0,2000,CorrectedE,
			   400,0,2000,CascadeZ);
		      if(gr->GetHit(g1)->GetTime()>-60.0 && gr->GetHit(g1)->GetTime()<-10.0){
			Fill(Form("CorrectedEold_mult2_GoIsegmult1_GKsegmult2_rhoGoI180_rhoGK200_tgamn60n10_%s",oname.Data()),
			     2000,0,2000,CorrectedE);
			Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_rhoGoI180_rhoGK200_tgamn60n10_%s",oname.Data()),
			     2000,0,2000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		    if(Rho_Tar_GammaOfInterest<230.0){
		      Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_rhoGoI230_%s",oname.Data()),
			   2000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		    
		    //RE New in v4.1_LT ->   Using both rho conditions together resulted in no peak remaining.
		    //trying just one rho condition and the timing condition.  
		    if(Rho_Tar_GammaOfInterest<180.0){
		      if(gr->GetHit(g1)->GetTime()>-60.0 && gr->GetHit(g1)->GetTime()<-10.0){
			Fill(Form("CorrectedEold_mult2_GoIsegmult1_GKsegmult2_rhoGoI180_tgamn60n10_%s",oname.Data()),
			     2000,0,2000,CorrectedE);
			Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_rhoGoI180_tgamn60n10_%s",oname.Data()),
			     2000,0,2000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		    //RE New in v4.1_LT ->   Using both rho conditions together resulted in no peak remaining.
		    //trying just one rho condition and the timing condition.  
		    if(Rho_Tar_GammaKnown<200.0){
		      if(gr->GetHit(g1)->GetTime()>-60.0 && gr->GetHit(g1)->GetTime()<-10.0){
			Fill(Form("CorrectedEold_mult2_GoIsegmult1_GKsegmult2_rhoGK200_tgamn60n10_%s",oname.Data()),
			     2000,0,2000,CorrectedE);
			Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult2_rhoGK200_tgamn60n10_%s",oname.Data()),
			     2000,0,2000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		  }
		  
		  if(gr->GetHit(g1)->GetSegMult()==3){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult3_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    if(Rho_Tar_GammaOfInterest<230.0){
		      Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult3_rhoGoI230_%s",oname.Data()),
			   2000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }
		  
		  if(gr->GetHit(g1)->GetSegMult()==4){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult4_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		  //yet, "high energy" seg multiplicty = 1 is not unheard of - photo is about 1/10 compt scatt prob.
		  //so keep track of that too.  
		  if(gr->GetHit(g1)->GetSegMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    //RE New in v4.1_LT ->  Want to make same gate for 885.3 seg mult = 1 case
		    //because there ought to be many true events here as well.  
		    if(Rho_Tar_GammaOfInterest<230.0){     
		      Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult1_GKSegmult1_rhoGoI230_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }
		}
		
		//"low energy" seg multiplicity = 2 (since some decent number of physical events compton scatter.)
		if(gr->GetHit(g2)->GetSegMult()==2){
		  Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoIsegmult2_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		  
		  Fill(Form("GKSegmult_vs_CorrectedEold_mult2_GoIsegmult2_%s",oname.Data()),        
		       1000,0,4000,CorrectedE,
		       20,0,20,gr->GetHit(g1)->GetSegMult());
		  
		  //further, "high energy" seg multiplicity = 2,3,4 (since physically it should compton scatter some
		  if(gr->GetHit(g1)->GetSegMult()==2){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult2_GKSegmult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		  
		  if(gr->GetHit(g1)->GetSegMult()==3){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult2_GKSegmult3_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		  
		  if(gr->GetHit(g1)->GetSegMult()==4){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult2_GKSegmult4_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		  //yet, "high energy" seg multiplicty = 1 is not unheard of - photo is about 1/10 compt scatt prob.
		  //so keep track of that too.  
		  if(gr->GetHit(g1)->GetSegMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedEold_mult2_GoISegmult2_GKSegmult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    //RE New in v4.1_LT ->  adding a couple most histos here since this is likely the best condition for
		    //using 172 to doppler correct the 885.  
		    Fill(Form("RhoGoI_vs_CorrectedEold_mult2_GoIsegmult2_GKsegmult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaOfInterest);
		    Fill(Form("RhoGK_vs_CorrectedEold_mult2_GoIsegmult2_GKsegmult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaKnown);
		  }
		}
	      } //close "if(mult == 2)"
	    } //close "if(g2!=g1)"
	  } //close "for all g2"    
	} //close "if g1 w/in labframe range"
      } //close "for all g1"
    } //closing the "if crystal mult >= 2" condition.  

    //RE New in v4.1_LT ->  Making some counters to keep track of the number of times
    //some histograms are filled per event.
    int counter_ABold_mult2_Zgate=0,counter_ABold_mult2_goi20_Zgate=0,counter_ABold_mult2_gkIPAngDiff0p7_goi20_gk60_Zgate=0,counter_ABold_mult2_gkIPAngDiff0p7_goi20_gk20t60_Zgate=0,counter_ABold_mult2_goi20_gk20t60_Zgate=0,counter_ABold_mult2_goi20_gk20t60=0,counter_ABold_mult2_gk20t60_Zgate=0;
      
    //RE New in v4.1_LT ->  Adding the AB version of requiring at least two crystal hits.
    //RE New in v4.1_LT ->   Only run of events that have exactly 1 gamma consistent with 885.  This gaurentees that decay
    //position is only calculated once.  So one event, one fill of decay pos.
    //    if(multABgoodgk==1){
    if(gr->GetMultAB()>=2){

      //RE New in v4.1_LT ->  I only want to determine the decay position once per event based on the chosenGood885 gamma.
      //So instead of looping over all gammas, I just use that one instead.  
      for(UShort_t g1=0;g1<gr->GetMultAB();g1++){

	if((gr->GetHitAB(g1)->GetEnergy()>LabFrameLoE_gk && gr->GetHitAB(g1)->GetEnergy()<LabFrameHiE_gk)){
	  //RE New in v4.1_LT -> Here, just use my one chosenGood885.  
	  if(g1==chosenGoodgk){
	    //RE New in v4.1_LT -> Using my new function to calculate the decay position.  
	    TVector3 decay = CascadeZdecay(gr->GetHitAB(g1),True_GammaKnown,beta_r);
	    CascadeZ = decay.Z();

	    //RE New in v4.1_LT -> Currently, this is only valid for the situation based on the not-rearranged hit order.
	    //RE New in v4.1_LT ->  Getting the hit position from the segment with the highest energy deposited
	    //(correct for segment-level analysis) instead of from the default position with is from the highest-energy IP
	    //(not necessarily in the highest-energy segment).
	  
	    Rho_Tar_GammaKnown = rho(gr->GetHitAB(g1)->GetPosition());

	    //RE New in v4.1_LT ->  Using span gate without AddBack makes for strange cases where the sphere intersects
	    //the edge of the crystal and is effectively cut off (since all the hits in a crystal are summed and called
	    //"one gamma ray" in this method.  Instead, try this span idea with Charlie's Sphere AddBack.
	    //The span is essentially the minimum radius that would keep all these hits together with Charlie's AddBack.  
	    double gkSpan = CalcSpan(gr->GetHitAB(g1));
	  
	    //RE New in v4.1_LT ->  Caluclate the difference between compton and geometric angles for THIS ordering.
	    double gkSegAngDiff = 10;
	    double gkIPAngDiff = 10;
	  
	    if(gr->GetHitAB(g1)->GetSegMult()>1){
	      gkSegAngDiff = abs(CasZangle(gr->GetHitAB(g1)->GetSegPos(0),
					   gr->GetHitAB(g1)->GetSegPos(1),
					   decay)
				 -acos(CosComptAngle(gr->GetHitAB(g1)->GetSegEn(0),
						     gr->GetHitAB(g1)->GetEnergy()-gr->GetHitAB(g1)->GetSegEn(0))));
	    }
	    else{
	      gkSegAngDiff=0;     //no possible compton scatter to campare this with.
	    }
	  
	    if(gr->GetHitAB(g1)->GetIPMult()>1){
	      gkIPAngDiff = abs(CasZangle(gr->GetHitAB(g1)->GetIPoints()[0]->GetPosition(),
					  gr->GetHitAB(g1)->GetIPoints()[1]->GetPosition(),
					  decay)
				-acos(CosComptAngle(gr->GetHitAB(g1)->GetIPoints()[0]->GetEnergy(),
						    gr->GetHitAB(g1)->GetEnergy()-gr->GetHitAB(g1)->GetIPoints()[0]->GetEnergy())));
	    }
	    else{
	      gkIPAngDiff=0;     //no possible compton scatter to campare this with.
	    }
	  
	    Fill(Form("CascadeZ_ABold_%s",oname.Data()),
		 400,0,2000,CascadeZ);
	  
	    Fill(Form("CascadeZ_thisegam_ABold_%s",oname.Data()),
		 2000,0,8000,Shifted_GammaKnown,
		 400,0,2000,CascadeZ);
	    Fill(Form("CascadeZ_thisegamdc_ABold_%s",oname.Data()),
		 2000,0,8000,gr->GetHitAB(g1)->GetDCEnergy(),
		 400,0,2000,CascadeZ);
	  
	    for(UShort_t g2=0;g2<gr->GetMultAB();g2++){

	      if(g2!=g1){

		//RE New in v4.1_LT ->  Using my new function to calculate the corrected energy. 
		CorrectedE = CascadeCorrectedE(gr->GetHitAB(g2),decay,beta_r);

		Rho_Tar_GammaOfInterest = rho(gr->GetHitAB(g2)->GetPosition());
		double goiSpan = CalcSpan(gr->GetHitAB(g2));
	       
		Fill(Form("CorrectedE_ABold_%s",oname.Data()),
		     2000,0,4000,CorrectedE);
	    
		//RE New in v4.1_LT ->  Changing the range of CascadeZ axis so that the 30P case (target is at GRETINA center can be plotted usefully.  Change back for 32Mg Upstream case.  
		Fill(Form("CascadeZ_vs_CorrectedE_ABold_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     600,-1000,2000,CascadeZ);

		//RE New in v4.1_LT ->  Trying the new tgam gate: -60 < tgam "885" < -40 && -60 < tgam "164" < 0.
		if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
		  if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		}
		  
		Fill(Form("thisegam_vs_CorrectedE_ABold_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,2000,Shifted_GammaKnown);
		  
		Fill(Form("otheregam_vs_CorrectedE_ABold_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     500,0,2000,Shifted_GammaOfInterest);
		  
		Fill(Form("thisegamdc_vs_CorrectedE_ABold_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy());
		  
		Fill(Form("othergamdc_vs_CorrectedE_ABold_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy());

		Fill(Form("mult_vs_CorrectedE_ABold_%s",oname.Data()), 
		     1000,0,4000,CorrectedE,
		     10,0,10,gr->GetMultAB());

		Fill(Form("CascadeZ_vs_CorrectedE_ABold_gkSpan20t60_%s",oname.Data()),
		     1000,0,4000,CorrectedE,
		     400,0,2000,CascadeZ);

		//RE New in v4.1_LT -> Checking best set of gates w/o mult=2 gate.  
		if(gkIPAngDiff<0.7){
		  if(goiSpan<20.0){
		    if(gkSpan>20.0 && gkSpan<60.0){
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);

		      //RE New in v4.1_LT -> Trying the new tgam gate: -60 < tgam "885" < -40 && -60 < tgam "164" < 0.
		      if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
			if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
			  Fill(Form("CascadeZ_vs_CorrectedE_ABold_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			       1000,0,4000,CorrectedE,
			       400,0,2000,CascadeZ);
			}
		      }
		    }
		  }
		}

		//RE New in v4.1_LT -> Checking spectrum without mult=2 and w/o AngDiff gate
		if(goiSpan<20.0){
		  if(gkSpan>20.0 && gkSpan<60.0){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_goiSpan20_gkSpan20t60_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);

		    //RE New in v4.1_LT -> See corrE spectrum with varying upper limit for "both" timing gate.
		    //3d histo  - tgamGK_tgamGOI_corrE
		    if(CascadeZ>522.0 && CascadeZ<922.0){
		      Fill(Form("CorrectedE_vs_tgamgk_vs_tgamgoi_ABold_goiSpan20_gkSpan20t60_Zgate_%s",oname.Data()),
			   500,0,4000,CorrectedE,
			   80,-200,200,gr->GetHitAB(g1)->GetTime(),
			   80,-200,200,gr->GetHitAB(g2)->GetTime());
		    }
		  
		    //RE New in v4.1_LT -> Trying the new tgam gate: -60 < tgam "885" < -40 && -60 < tgam "164" < 0.
		    if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
		      if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_goiSpan20_gkSpan20t60_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		  }
		  if(gkSpan<60.0){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_goiSpan20_gkSpan60_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		}
		if(goiSpan<20.0){
		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_goiSpan20_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
		if(gkSpan>20.0 && gkSpan<60.0){
		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_gkSpan20t60_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);

		  //RE New in v4.1_LT -> Trying the new tgam gate: -60 < tgam "885" < -40 && -60 < tgam "164" < 0.
		  if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
		    if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_gkSpan20t60_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }
		}

		//RE New in v4.1_LT -> Studying how peak changes as mult gate is gradually removed - without other gates
	      
		if(gr->GetMultAB()<=2){
		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2orLess_%s",oname.Data()),
		       4000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
		if(gr->GetMultAB()<=3){
		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult3orLess_%s",oname.Data()),
		       4000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
		if(gr->GetMultAB()<=4){
		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult4orLess_%s",oname.Data()),
		       4000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
		if(gr->GetMultAB()<=5){
		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult5orLess_%s",oname.Data()),
		       4000,0,4000,CorrectedE,
		       400,0,2000,CascadeZ);
		}
	      
		//RE New in v4.1_LT -> Studying how peak changes as mult gate is gradually removed.
		if(gkIPAngDiff<0.7){
		  if(goiSpan<20.0){
		    if(gkSpan>20.0 && gkSpan<60.0){
		      if(gr->GetMultAB()<=2){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2orLess_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_%s",oname.Data()),
			     4000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		      if(gr->GetMultAB()<=3){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult3orLess_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_%s",oname.Data()),
			     4000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		      if(gr->GetMultAB()<=4){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult4orLess_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_%s",oname.Data()),
			     4000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		      if(gr->GetMultAB()<=5){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult5orLess_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_%s",oname.Data()),
			     4000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		      if(CascadeZ>522.0 && CascadeZ<922.0){
			Fill(Form("multAB_vs_CorrectedE_ABold_gkIPAngDiff0p7_goiSpan20_gkSpan20t60_Zgate_%s",oname.Data()),
			     4000,0,4000,CorrectedE,
			     20,0,20,gr->GetMultAB());
		      }
		    }
		  }
		}
	      
		//RE New in v4.1_LT -> Want to see cases where IPMult=1 for both gammas (g1 and g2).
		//In these cases, our choice of the "first interaction point" does not matter.
		if(gr->GetHitAB(g1)->GetIPMult()==1){
		  if(gr->GetHitAB(g2)->GetIPMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_goiIPMult1_gkIPMult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		}
		if(gr->GetHitAB(g1)->GetSegMult()==1){
		  if(gr->GetHitAB(g2)->GetSegMult()==1){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_goiSegMult1_gkSegMult1_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		  }
		}
	      
		//RE New in v4.1_LT -> Now properly uses GetMultAB, not GetMult.  
		if(gr->GetMultAB()==2){


		  //RE New in v4.1_LT -> Want to see cases where IPMult=1 for both gammas (g1 and g2).
		  //In these cases, our choice of the "first interaction point" does not matter.
		  if(gr->GetHitAB(g1)->GetIPMult()==1){
		    if(gr->GetHitAB(g2)->GetIPMult()==1){
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiIPMult1_gkIPMult1_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }
		  if(gr->GetHitAB(g1)->GetSegMult()==1){
		    if(gr->GetHitAB(g2)->GetSegMult()==1){
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiSegMult1_gkSegMult1_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }

		
		  Fill(Form("dta_vs_CorrectedE_ABold_mult2_%s",oname.Data()),    //RE New in v4.1_LT -> Now properly fills AB
		       1000,0,4000,CorrectedE,                                   //spectrum instead of the corresponding
		       200,-10,10,track->GetDTA());                              //no-AB spectrum.  

		  //RE New in v4.1_LT -> Adding histos for span and rho gates.
		  if(Rho_Tar_GammaOfInterest<230.){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiRho230_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    Fill(Form("CorrectedE_ABold_mult2_goiRho230_%s",oname.Data()),
			 2000,0,4000,CorrectedE);
		  }
		  //RE New in v4.1_LT 
		  if(CascadeZ>522.0 && CascadeZ<922.0){
		    Fill(Form("gkSpan_vs_CorrectedE_ABold_mult2_p522p922_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 200,0,200,gkSpan);
		  }
		
		  if(goiSpan<20.0){  //in mm.  
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiSpan20_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);
		    Fill(Form("CorrectedE_ABold_mult2_goiSpan20_%s",oname.Data()),
			 2000,0,4000,CorrectedE);
		    //RE New in v4.1_LT -> energy-timing spectrum with goi gate.  
		    Fill(Form("CorrectedE_vs_gkTime_ABold_mult2_goiSpan20_%s",oname.Data()),
			 400,-200,200,gr->GetHitAB(g1)->GetTime(),
			 4000,0,4000,CorrectedE);
		    Fill(Form("CorrectedE_vs_goiTime_ABold_mult2_goiSpan20_%s",oname.Data()),
			 400,-200,200,gr->GetHitAB(g2)->GetTime(),
			 4000,0,4000,CorrectedE);
		    if(Rho_Tar_GammaOfInterest<230.){
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiRho230_goiSpan20_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);		    
		    }
		    //RE New in v4.1_LT -> Want to test a range of limits on the gk sphere radius.
		    if(CascadeZ>522.0 && CascadeZ<922.0){
		      Fill(Form("gkSpan_vs_CorrectedE_ABold_mult2_goiSpan20_p522p922_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   200,0,200,gkSpan);

		      //RE New in v4.1_LT -> Counter for this condition.
		      counter_ABold_mult2_goi20_Zgate++;
		      Fill(Form("counter_ABold_mult2_goi20_Zgate_%s",oname.Data()),
			   5,0,5,counter_ABold_mult2_goi20_Zgate);
		    }

		    //RE New in v4.1_LT -> timing gates but no gkSpan gate, no compt gate
		    if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
		      if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiSpan20_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		  
		  }
		  if(gkSpan>20.0 && gkSpan<60.0){
		    Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_gkSpan20t60_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 400,0,2000,CascadeZ);

		    //RE New in v4.1_LT ->Trying the new tgam gate: -60 < tgam "885" < -40 && -60 < tgam "164" < 0.
		    if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
		      if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_gkSpan20t60_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		  
		  
		    Fill(Form("CorrectedE_ABold_mult2_gkSpan20t60_%s",oname.Data()),
			 2000,0,4000,CorrectedE);
		    //RE New in v4.1_LT -> Want to test a range of limits on the goi sphere radius.
		    if(CascadeZ>522.0 && CascadeZ<922.0){
		      Fill(Form("goiSpan_vs_CorrectedE_ABold_mult2_gkSpan20t60_p522p922_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   200,0,200,goiSpan);

		      //RE New in v4.1_LT -> counter here
		      counter_ABold_mult2_gk20t60_Zgate++;
		      Fill(Form("counter_ABold_mult2_gk20t60_Zgate_%s",oname.Data()),
			   5,0,5,counter_ABold_mult2_gk20t60_Zgate);
		    }
		    if(goiSpan<20.0){  //in mm.  
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiSpan20_gkSpan20t60_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);

		      //RE New in v4.1_LT -> Trying the new tgam gate: -60 < tgam "885" < -40 && -60 < tgam "164" < 0.
		      if(gr->GetHitAB(g1)->GetTime() > -60.0 && gr->GetHitAB(g1)->GetTime() < -40.0){
			if(gr->GetHitAB(g2)->GetTime() > -60.0 && gr->GetHitAB(g2)->GetTime() < 0.0){
			  Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiSpan20_gkSpan20t60_tgamgkn60n40_tgamgoin60n00_%s",oname.Data()),
			       1000,0,4000,CorrectedE,
			       400,0,2000,CascadeZ);
			}
		      }

		      //RE New in v4.1_LT -> counter here
		      if(CascadeZ>522.0 && CascadeZ<922.0){
			counter_ABold_mult2_goi20_gk20t60_Zgate++;
			Fill(Form("counter_ABold_mult2_goi20_gk20t60_Zgate_%s",oname.Data()),
			     5,0,5,counter_ABold_mult2_goi20_gk20t60_Zgate);
		      }

		      //RE New in v4.1_LT -> counter here
		      counter_ABold_mult2_goi20_gk20t60++;
		      Fill(Form("counter_ABold_mult2_goi20_gk20t60_%s",oname.Data()),
			   5,0,5,counter_ABold_mult2_goi20_gk20t60);
		    
		      if(Rho_Tar_GammaOfInterest<230.){
			Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiRho230_goiSpan20_gkSpan20t60_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     400,0,2000,CascadeZ);
		      }
		    }
		    if(Rho_Tar_GammaOfInterest<230.){
		      Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_goiRho230_gkSpan20t60_%s",oname.Data()),
			   1000,0,4000,CorrectedE,
			   400,0,2000,CascadeZ);
		    }
		  }
		  		
		  if(gr->GetHitAB(g1)->GetIPMult()>1){
		    
		    Fill(Form("gkIPAngDiff_vs_CorrectedE_ABold_mult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 1000,0,3.2,gkIPAngDiff);
		    
		  }
		
		  if(gr->GetHitAB(g1)->GetSegMult()>1){                                  
			  
		    Fill(Form("gkSegAngDiff_vs_CorrectedE_ABold_mult2_%s",oname.Data()),
			 1000,0,4000,CorrectedE,
			 1000,0,3.2,gkSegAngDiff);
		  
		    //RE New in v4.1_LT -> adding another look at compt gate and span gate.
		    if(CascadeZ>522.0 && CascadeZ<922.0){
		      if(goiSpan<20.0){
			if(gkSpan>20.0 && gkSpan<60.0){
			  Fill(Form("gkSegAngDiff_vs_CorrectedE_ABold_mult2_goiSpan20_gkSpan20t60_p522p922_%s",oname.Data()),
			       1000,0,4000,CorrectedE,
			       1000,0,3.2,gkSegAngDiff);
			}
		      }
		      if(gkSpan>20.0 && gkSpan<60.0){
			Fill(Form("gkSegAngDiff_vs_CorrectedE_ABold_mult2_gkSpan20t60_p522p922_%s",oname.Data()),
			     1000,0,4000,CorrectedE,
			     1000,0,3.2,gkSegAngDiff);
		      }
		    }
		  }
		    
		  Fill(Form("CorrectedE_ABold_mult2_%s",oname.Data()),
		       2000,0,4000,CorrectedE);

		  Fill(Form("CascadeZ_vs_CorrectedE_ABold_mult2_%s",oname.Data()),
		       1000,0,4000,CorrectedE,
		       600,-1000,2000,CascadeZ);

		  Fill(Form("GoISegmult_vs_CorrectedE_ABold_mult2_%s",oname.Data()),  
		       1000,0,4000,CorrectedE,
		       20,0,20,gr->GetHitAB(g1)->GetSegMult());
				    
		  Fill(Form("GKSegmult_vs_CorrectedE_ABold_mult2_%s",oname.Data()),   
		       1000,0,4000,CorrectedE,
		       20,0,20,gr->GetHitAB(g1)->GetSegMult());

		  Fill(Form("RhoGoI_vs_CorrectedE_ABold_mult2_%s",oname.Data()),      
		       1000,0,4000,CorrectedE,
		       500,0,500,Rho_Tar_GammaOfInterest);

		    
		  if(Rho_Tar_GammaKnown<=200.){
		    Fill(Form("RhoGoI_vs_CorrectedE_ABold_mult2_rhoGK200%s",oname.Data()),   
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaOfInterest);
		  }

		  Fill(Form("RhoGK_vs_CorrectedE_ABold_mult2_%s",oname.Data()),     
		       1000,0,4000,CorrectedE,
		       500,0,500,Rho_Tar_GammaKnown);

		    
		  if(Rho_Tar_GammaOfInterest<=180.){
		    Fill(Form("RhoGK_vs_CorrectedE_ABold_mult2_rhoGoI180%s",oname.Data()),  
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaKnown);
		  }

		  if(gr->GetHitAB(g2)->GetIPMult()==1){
		    Fill(Form("RhoGoI_vs_CorrectedE_ABold_mult2_GoIipmult1_%s",oname.Data()),  
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaOfInterest);
		  }
		  else if(gr->GetHitAB(g2)->GetIPMult()==2){
		    Fill(Form("RhoGoI_vs_CorrectedE_ABold_mult2_GoIipmult2_%s",oname.Data()),  
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaOfInterest);
		  }
		    
		  else if(gr->GetHitAB(g2)->GetIPMult()==3){
		    Fill(Form("RhoGoI_vs_CorrectedE_ABold_mult2_GoIipmult3_%s",oname.Data()),  
			 1000,0,4000,CorrectedE,
			 500,0,500,Rho_Tar_GammaOfInterest);
		  }

		} //closing "if(mult==2)
	      } //closing "if g2 != g1"
	      //	    } //closing "if g2 is w/in good timing region".  
	    } //closing "for all g2"
	  } //closing "if g1==chosenGood885"
	} //closing "if g1 w/in labframe range"
	//	} //closing "if g1 is w/in good tgam region".  
      } //closing "for all g1"
    }  //closing the "if AB crystal mult >= 2" condition.
    
  }//closing "upstream"
  
  for(UShort_t g1=0;g1<gr->GetMultAB();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMultAB();g2++){
      Fill(Form("egamegamABdc_%s",oname.Data()),
	   1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy(),
	   1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy());
      Fill(Form("egamegamABdc_sym_%s",oname.Data()),
	   1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy(),
	   1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy());
      Fill(Form("egamegamABdc_sym_%s",oname.Data()),
	   1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy(),
	   1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy());
    }
  }


  // Various gamma-gamma spectra with an additional time gate.
  {
    int highestg=-1;
    double highesten =0;
    for(UShort_t g=0;g<gr->GetMult();g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(hit->GetDCEnergy()>highesten){
	  highesten = hit->GetDCEnergy();
	  highestg =g;
	}
      }
    }
    for(UShort_t g=0;g<gr->GetMult();g++){
      HitCalc* hit = gr->GetHit(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(highestg>-1 && g!=highestg){
	  Fill(Form("egamegamdc_fold_%s_tcut",oname.Data()),
	       4000,0,4000,highesten,
	       4000,0,4000,hit->GetDCEnergy());
	}
      }
    }
    highestg=-1;
    highesten =0;
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      HitCalc* hit = gr->GetHitAB(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(hit->GetDCEnergy()>highesten){
	  highesten = hit->GetDCEnergy();
	  highestg =g;
	}
      }
    }
    for(UShort_t g=0;g<gr->GetMultAB();g++){
      HitCalc* hit = gr->GetHitAB(g);
      if(hit->GetEnergy()<1 || hit->GetEnergy()>6000)
	continue;
      if(foundTCut && TimeCut->IsInside(hit->GetTime(),hit->GetEnergy())){
	if(highestg>-1 && g!=highestg){
	  Fill(Form("egamegamABdc_fold_%s_tcut",oname.Data()),
	       4000,0,4000,highesten,
	       4000,0,4000,hit->GetDCEnergy());
	}
      }
    }
  }
  for(UShort_t g1=0;g1<gr->GetMult();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMult();g2++){
      if(foundTCut && TimeCut->IsInside(gr->GetHit(g1)->GetTime(),gr->GetHit(g1)->GetEnergy()) &&
	 TimeCut->IsInside(gr->GetHit(g2)->GetTime(),gr->GetHit(g2)->GetEnergy())){
	Fill(Form("egamegamdc_%s_tcut",oname.Data()),
	     1000,0,4000,gr->GetHit(g1)->GetDCEnergy(),
	     1000,0,4000,gr->GetHit(g2)->GetDCEnergy());
	if(gr->GetHit(g1)->GetPosition().Theta()>1.03 && gr->GetHit(g2)->GetPosition().Theta()>1.03){
	  Fill(Form("egamegamdc_%s_tcut_thetacut",oname.Data()),
	       1000,0,4000,gr->GetHit(g1)->GetDCEnergy(),
	       1000,0,4000,gr->GetHit(g2)->GetDCEnergy());
	}
	Fill(Form("egamegamdc_sym_%s_tcut",oname.Data()),
	     1000,0,4000,gr->GetHit(g1)->GetDCEnergy(),
	     1000,0,4000,gr->GetHit(g2)->GetDCEnergy());
	Fill(Form("egamegamdc_sym_%s_tcut",oname.Data()),
	     1000,0,4000,gr->GetHit(g2)->GetDCEnergy(),
	     1000,0,4000,gr->GetHit(g1)->GetDCEnergy());
      }
    }
  }
  for(UShort_t g1=0;g1<gr->GetMultAB();g1++){
    for(UShort_t g2=g1+1;g2<gr->GetMultAB();g2++){
      if(foundTCut && TimeCut->IsInside(gr->GetHitAB(g1)->GetTime(),gr->GetHitAB(g1)->GetEnergy()) &&
	 TimeCut->IsInside(gr->GetHitAB(g2)->GetTime(),gr->GetHitAB(g2)->GetEnergy())){
	Fill(Form("egamegamABdc_%s_tcut",oname.Data()),
	     1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy(),
	     1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy());
	if(gr->GetHitAB(g1)->GetPosition().Theta()>1.03 && gr->GetHitAB(g2)->GetPosition().Theta()>1.03){
	  Fill(Form("egamegamABdc_%s_tcut_thetacut",oname.Data()),
	       1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy(),
	       1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy());
	}
	Fill(Form("egamegamABdc_sym_%s_tcut",oname.Data()),
	     1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy(),
	     1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy());
	Fill(Form("egamegamABdc_sym_%s_tcut",oname.Data()),
	     1000,0,4000,gr->GetHitAB(g2)->GetDCEnergy(),
	     1000,0,4000,gr->GetHitAB(g1)->GetDCEnergy());
      }
    }
  }

  int hp = s800->GetRegistr();
  Fill(Form("registr_%s",oname.Data()),
       20,-0.5,19.5,hp);
  for(int j=0;j<16;j++){
    if(hp & (1<<j)){
      Fill(Form("trigbit_%s",oname.Data()),
	   16,-0.5,15.5,j);
    }
  }

  if(fCal&(1<<2)){
    for(UShort_t p=0; p<2;p++){
      Fill(Form("ICde_vs_x%d_%s",p,oname.Data()),
    	   600,-300,300,pad[p]->GetX(),
    	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
      Fill(Form("ICde_vs_y%d_%s",p,oname.Data()),
    	   600,-100,100,pad[p]->GetY(),
    	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetDE());
      Fill(Form("ICsum_vs_x%d_%s",p,oname.Data()),
      	   600,-300,300,pad[p]->GetX(),
      	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
      Fill(Form("ICsum_vs_y%d_%s",p,oname.Data()),
      	   600,-100,100,pad[p]->GetY(),
      	   fIC_range[0],fIC_range[1],fIC_range[2],ich->GetSum());
    }
    Fill(Form("x_vs_obj_%s",oname.Data()),
	 fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_obj_%s",oname.Data()),
	 fobj_range[0],fobj_range[1],fobj_range[2],tof->GetOBJ(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_objC_%s",oname.Data()),
	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objC_%s",oname.Data()),
	 fobjC_range[0],fobjC_range[1],fobjC_range[2],tof->GetOBJC(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfp_%s",oname.Data()),
	 fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfp_%s",oname.Data()),
	 fxfp_range[0],fxfp_range[1],fxfp_range[2],tof->GetXFP(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfpC_%s",oname.Data()),
	 fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfpC_%s",oname.Data()),
	 fxfpC_range[0],fxfpC_range[1],fxfpC_range[2],tof->GetXFPC(),
	 100/2,-100,100,track->GetAFP());

    Fill(Form("x_vs_objtac_%s",oname.Data()),
	 ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objtac_%s",oname.Data()),
	 ftacobj_range[0],ftacobj_range[1],ftacobj_range[2],tof->GetTACOBJ(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_objtacC_%s",oname.Data()),
	 ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objtacC_%s",oname.Data()),
	 ftacobjC_range[0],ftacobjC_range[1],ftacobjC_range[2],tof->GetTACOBJC(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfptac_%s",oname.Data()),
	 ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfptac_%s",oname.Data()),
	 ftacxfp_range[0],ftacxfp_range[1],ftacxfp_range[2],tof->GetTACXFP(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfptacC_%s",oname.Data()),
	 ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfptacC_%s",oname.Data()),
	 ftacxfpC_range[0],ftacxfpC_range[1],ftacxfpC_range[2],tof->GetTACXFPC(),
	 100/2,-100,100,track->GetAFP());

    Fill(Form("x_vs_objm_%s",oname.Data()),
	 fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJ(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objm_%s",oname.Data()),
	 fmobj_range[0],fmobj_range[1],fmobj_range[2],tof->GetMOBJ(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_objmC_%s",oname.Data()),
	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_objmC_%s",oname.Data()),
	 fmobjC_range[0],fmobjC_range[1],fmobjC_range[2],tof->GetMOBJC(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfpm_%s",oname.Data()),
	 fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFP(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfpm_%s",oname.Data()),
	 fmxfp_range[0],fmxfp_range[1],fmxfp_range[2],tof->GetMXFP(),
	 100/2,-100,100,track->GetAFP());
    Fill(Form("x_vs_xfpmC_%s",oname.Data()),
	 fmxfpC_range[0],fmxfpC_range[1],fmxfpC_range[2],tof->GetMXFPC(),
	 600/2,-300,300,pad[0]->GetX());
    Fill(Form("afp_vs_xfpmC_%s",oname.Data()),
	 fmxfpC_range[0],fmxfpC_range[1],fmxfpC_range[2],tof->GetMXFPC(),
	 100/2,-100,100,track->GetAFP());
  }

  Fill(Form("xfp_%s",oname.Data()),
       600,-300,300,track->GetXFP());
  Fill(Form("xfpf_%s",oname.Data()),
       600,-0.3,0.3,track->GetXFP()/1000.);
  Fill(Form("xfpafp_%s",oname.Data()),
       200,-100,100,track->GetAFP(),
       600,-300,300,track->GetXFP());
  Fill(Form("xfpazita_%s",oname.Data()),
       360,0,360,track->GetPhi(),
       600,-300,300,track->GetXFP());
  Fill(Form("afp_%s",oname.Data()),
       400,-200,200,track->GetAFP());
  Fill(Form("yfp_%s",oname.Data()),
       200,-100,100,track->GetYFP());
  Fill(Form("bfp_%s",oname.Data()),
       100,-50,50,track->GetBFP());
  Fill(Form("yta_%s",oname.Data()),
       100,-50,50,track->GetYTA());
  Fill(Form("dta_%s",oname.Data()),
       200,-10,10,track->GetDTA());
  // Fill(Form("dtac_%s",oname.Data()),
  //      200,-10,10,track->GetDTAC());
  Fill(Form("azita_%s",oname.Data()),
       360,0,360,track->GetPhi());
  Fill(Form("ptot_%s",oname.Data()),
       500,8,20,track->GetPtot());
  Fill(Form("ppar_%s",oname.Data()),
       fPP_range[0],fPP_range[1],fPP_range[2],track->GetPpar());
  Fill(Form("pparc_%s",oname.Data()),
       fPP_range[0],fPP_range[1],fPP_range[2],track->GetPparC());
  Fill(Form("ptra_%s",oname.Data()),
       200,0,2,track->GetPtra());
  Fill(Form("etot_%s",oname.Data()),
       1000,2000,3000,track->GetEtot());
  Fill(Form("dta_vs_ata_%s",oname.Data()),
       120,-60,60,track->GetATA(),
       400,-10,10,track->GetDTA());
  // Fill(Form("dta_vs_dtac_%s",oname.Data()),
  //      400,-10,10,track->GetDTAC(),
  //      400,-10,10,track->GetDTA());


#ifdef S800_DETAILEDTREE
  if(fCal&(1<<1)){
    for(UShort_t c=0;c<ich->GetChan().size();c++){
      Fill(Form("ic_ch%d_%s",ich->GetChan()[c],oname.Data()),
	   4000,0,4000,ich->GetCal()[c]);
      Fill(Form("ic_%s",oname.Data()),
	   16,0,16,ich->GetChan()[c],
	   2500,0,2500,ich->GetCal()[c]);
    }
  }
  if(fCal&(1<<0)) {
    for(UShort_t p=0; p<2;p++){
      if(pad[p]->GetMaxPad()>0&&pad[p]->GetMaxPad()<224){
	Fill(Form("crdcmaxpad_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[1]);

	Fill(Form("crdcpad_%d_%d_%s",p,pad[p]->GetMaxPad(),oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[1]);

	//max and neighbours

	Fill(Form("crdcmaxpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[1]);
	Fill(Form("crdcmaxpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad()-1,
	     1024,0,1024,pad[p]->GetPadMax()[0]);
	Fill(Form("crdcmaxpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad()+1,
	     1024,0,1024,pad[p]->GetPadMax()[2]);

	Fill(Form("crdcpad0_%d_%d_%s",p,pad[p]->GetMaxPad(),oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[1]);
	Fill(Form("crdcpad0_%d_%d_%s",p,pad[p]->GetMaxPad()-1,oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[0]);
	Fill(Form("crdcpad0_%d_%d_%s",p,pad[p]->GetMaxPad()+1,oname.Data()),
	     1024,0,1024,pad[p]->GetPadMax()[2]);


	Fill(Form("xpad_%d_%s",p,oname.Data()),
	     200,-2,2,pad[p]->GetXGravity()-pad[p]->GetXFit(),
	     224,0,224,pad[p]->GetMaxPad());

	Fill(Form("xpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[0]);
	Fill(Form("xpad0_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1024,0,1024,pad[p]->GetPadMax()[2]);

	Fill(Form("xpadl_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1000,-500,500,pad[p]->GetPadMax()[1]-pad[p]->GetPadMax()[0]);
	Fill(Form("xpadr_%d_%s",p,oname.Data()),
	     224,0,224,pad[p]->GetMaxPad(),
	     1000,-500,500,pad[p]->GetPadMax()[1]-pad[p]->GetPadMax()[2]);

	Fill(Form("xgra_%d_%s",p,oname.Data()),
	     200,-2,2,pad[p]->GetXGravity()-pad[p]->GetMaxPad(),
	     224,0,224,pad[p]->GetMaxPad());

	Fill(Form("xfit_%d_%s",p,oname.Data()),
	     200,-2,2,pad[p]->GetXFit()-pad[p]->GetMaxPad(),
	     224,0,224,pad[p]->GetMaxPad());
      }//pads
    }//crdc

  }
#endif

}
void CalHistograms::FillHistogramsGateOnlyOut(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, const char* outname){
  //in case you have pure beam, fill here
}
