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
#include <string>
#include <vector>
#include <sys/time.h> 
#include <signal.h>
#include <omp.h>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TKey.h"

#include "Gretinadefs.h"

#include "CommandLineInterface.hh"
#include "Settings.hh"
#include "RunInfo.hh"
#include "Calibration.hh"
#include "Gretina.hh"
#include "GretinaCalc.hh"
#include "GretinaTrack.hh"
#include "Mode3Calc.hh"
#include "Trace.hh"
#include "S800.hh"
#include "S800Calc.hh"

using namespace TMath;
using namespace std;

bool signal_received = false;
void signalhandler(int sig){
  if (sig == SIGINT){
    signal_received = true;
  }
}
double get_time(){  
    struct timeval t;  
    gettimeofday(&t, NULL);  
    double d = t.tv_sec + (double) t.tv_usec/1000000;  
    return d;  
}  
int main(int argc, char* argv[]){
  double time_start = get_time();  
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);

  char *InputFile = NULL;
  char *RootFile = NULL;
  char* SettingFile;
  char* CutFile = NULL;
  int LastEvent = -1;
  int addbackType = -1;
  bool trackMe = false;
  int tac = 0;
  CommandLineInterface* interface = new CommandLineInterface();
  interface->Add("-i", "input file", &InputFile);
  interface->Add("-o", "output file", &RootFile);
  interface->Add("-s", "settingsfile", &SettingFile);
  interface->Add("-c", "cutfile", &CutFile);
  interface->Add("-t", "0 for OBJTAC, 1 for OBJTDC, 2 for XFPTAC, 1 for XFPTDC (default = 0)", &tac);
  interface->Add("-n", "last event", &LastEvent);
  interface->Add("-ab","type of addback to use (overrides settings file)", &addbackType);
  interface->Add("-tr","track (overrides settings file)", &trackMe);
  interface->CheckFlags(argc, argv);


  if(InputFile == NULL){
    cout << "No input file given " << endl;
    return 1;
  }
  if(RootFile == NULL){
    cout << "No output ROOT file given " << endl;
    return 2;
  }
  if(SettingFile == NULL){
    cout << "No settings file given " << endl;
    return 3;
  }




  TFile* infile = new TFile(InputFile);
  TTree* gtr = (TTree*) infile->Get("gtr");
  if(gtr == NULL){
    cout << "could not find tree tr in file " << infile->GetName() << endl;
    return 3;
  }

  Mode3Event* m3e = new Mode3Event;
  gtr->SetBranchAddress("mode3Event",&m3e);
  S800 *s800 = new S800;
  gtr->SetBranchAddress("s800",&s800);
  Gretina* gretina = new Gretina;
  gtr->SetBranchAddress("gretina",&gretina);



  TFile *ofile = new TFile(RootFile,"RECREATE");
  cout<<"input file: "<<InputFile<< endl;
  cout<<"writing to output file: "<<RootFile<< endl;
  cout << "------------------------------------" << endl;
  Settings* set = new Settings(SettingFile);
  if(addbackType!=-1){
    cout << "Modified the addbackType to " << addbackType << endl;
    set->SetAddBackType(addbackType);
  }
  if(trackMe)
    cout << "Tracking is ON " << endl;
  else
    cout << "Tracking is OFF " << endl;
    
  set->SetTracking(trackMe);
  set->Write("settings",TObject::kOverwrite);

  RunInfo *info = (RunInfo*)infile->Get("runinfo");
  if(info==NULL){
    cout << InputFile << "does not contain run information" <<endl;
    info = new RunInfo;
  }
  TEnv *evtnumbers = new TEnv(set->EvtNrFile());
  //int startevent = evtnumbers->GetValue(Form("Start.Event.Number.%d",info->GetRunNr()),0);
  //int RunNumber = info->GetRunNr(); //AR New in v4.1_LT
  //get the run number from the filename  //AR New in v4.1_LT
  int run;
  TString ifname(InputFile);
  ifname.Remove(0,ifname.Length()-12); //cut to keep only RunXXXX.root
  sscanf(ifname.Data(),"Run%04d.root",&run);
  int RunNumber = run;
  int startevent = evtnumbers->GetValue(Form("Start.Event.Number.%d",RunNumber),0);
  cout << "run number " << RunNumber << " starts with event nr " << startevent << endl;   

  Calibration *cal = new Calibration(set,startevent,RunNumber); //AR New in v4.1_TL
  cout << "Brho " << cal->GetBrho() << " mass " <<  cal->GetMass() << " Z " << cal->GetZ() << endl;
  info->SetBeam(cal->GetBrho(), cal->GetMass(), cal->GetZ());

  TTree* ctr = new TTree("ctr","Gretina/S800 calibrated and builtevents");
  S800Calc* s800Calc = new S800Calc;
  GretinaCalc* gretinaCalc = new GretinaCalc;
  GretinaEvent* gretinaEvent = new GretinaEvent;
  Mode3Calc* mode3Calc = new Mode3Calc;
  ctr->Branch("s800calc",&s800Calc, 320000);
  ctr->Branch("gretinacalc",&gretinaCalc, 320000);
  if(trackMe)
    ctr->Branch("gretinaevent",&gretinaEvent, 320000);
  ctr->Branch("mode3calc",&mode3Calc, 320000);
  ctr->BranchRef();

  vector<TCutG*> InPartCut;
  vector<vector<TCutG*> > OutPartCut;
  vector<vector<TTree*> > splittree;


  if(CutFile!=NULL){
    cout << "Cuts were created for ";
    if(tac%2==0)
      cout << "TAC";
    else if (tac%2==1)
      cout << "TDC";
    cout << " data and are applied to the ";
    if(tac<2)
      cout << "OBJ";
    else
      cout << "XFP";
    cout << " scintillator" << endl;
    
    //Read in the cuts file for incoming and outgoing particle ID
    char* Name = NULL;
    char* Name2 = NULL;
    TFile *cFile = new TFile(CutFile);
    TIter next(cFile->GetListOfKeys());
    TKey *key;
    while((key=(TKey*)next())){
      if(strcmp(key->GetClassName(),"TCutG") == 0){
	Name = (char*)key->GetName();
	if(strstr(Name,"in") && !strstr(Name,"out")){
	  cout << "incoming cut found "<<Name << endl;
	  InPartCut.push_back((TCutG*)cFile->Get(Name));
	}
      }      
    }
    TIter next2(cFile->GetListOfKeys());
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
    cFile->Close();
    ofile->cd();
    splittree.resize(InPartCut.size());
    for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
      splittree[in].resize(OutPartCut[in].size());
      for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
	splittree[in][ou] = new TTree(Form("ctr_%s",OutPartCut[in][ou]->GetName()),"Gretina/S800 calibrated and builtevents");
	splittree[in][ou]->Branch("s800calc",&s800Calc, 320000);
	splittree[in][ou]->Branch("gretinacalc",&gretinaCalc, 320000);
	if(trackMe)
	  splittree[in][ou]->Branch("gretinaevent",&gretinaEvent, 320000);
	splittree[in][ou]->Branch("mode3calc",&mode3Calc, 320000);
	splittree[in][ou]->BranchRef();
      }
    }
  }
  Double_t nentries = gtr->GetEntries();
  Int_t nbytes = 0;
  Int_t status;

  cout << "Nr of Events " << nentries << endl;
  if (LastEvent!=-1){
    nentries = LastEvent;
  } else if(set->LastEvent()>0){
    nentries = set->LastEvent();
  }
  if(nentries != gtr->GetEntries()){
    cout << "reading only first " << nentries << " events " << endl;
  }


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This is a for loop that goes over all the enties. THIS WOULD NEED TO BE PARALLELIZED

  int max_threads = omp_get_max_threads();
  printf("MAX NUMBER OF THREADS: %d\n",max_threads);
  for(int i=0; i<nentries; i++){
    if(signal_received)
      break;

    if(set->VLevel()>2){
      cout << "-----------------------------------------"<< endl;
      cout << "processing event " << i << endl;
    }
    if(i%10000 == 0){
      //cal->PrintCtrs();
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries<<" % done\t" << 
      (Float_t)i/(time_end - time_start) << " events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
    if(i%1000000 == 0 && CutFile==NULL)
      ctr->AutoSave();

    //It clears all of the data entries in that object
    s800Calc->Clear();
    gretinaCalc->Clear();
    gretinaEvent->Clear();
    mode3Calc->Clear();
    s800->Clear();
    gretina->Clear();
    m3e->Clear();

    status = gtr->GetEvent(i);

    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<gtr->GetName()<<" in file "<<infile->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<gtr->GetName()<<" in file "<<infile->GetName()<<" doesn't exist"<<endl;
      return 6;
    }
    nbytes += status;



    //Build all of the calibrated objects, using the calibration in cal.
    //Use the data from the first three parameters, output into the last three parameters.

    cal->BuildAllCalc(s800,gretina,m3e,
      s800Calc,gretinaCalc,mode3Calc);

    if(trackMe)
      cal->GammaTrack(gretinaCalc,gretinaEvent);

    if(CutFile==NULL)
      ctr->Fill();
    else{
      for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
        if(tac%2 == 1 && !InPartCut[in]->IsInside(s800Calc->GetTOF()->GetOBJ(),s800Calc->GetTOF()->GetXFP()))
          continue;
        if(tac%2 == 0 && !InPartCut[in]->IsInside(s800Calc->GetTOF()->GetTACOBJ(),s800Calc->GetTOF()->GetTACXFP()))
          continue;
        for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
          if(tac == 1 && !OutPartCut[in][ou]->IsInside(s800Calc->GetTOF()->GetOBJC(),s800Calc->GetIC()->GetDE()))
            continue;
          if(tac == 0 && !OutPartCut[in][ou]->IsInside(s800Calc->GetTOF()->GetTACOBJC(),s800Calc->GetIC()->GetDE()))
            continue;
          if(tac == 3 && !OutPartCut[in][ou]->IsInside(s800Calc->GetTOF()->GetXFPC(),s800Calc->GetIC()->GetDE()))
            continue;
          if(tac == 2 && !OutPartCut[in][ou]->IsInside(s800Calc->GetTOF()->GetTACXFPC(),s800Calc->GetIC()->GetDE()))
            continue;
          splittree[in][ou]->Fill();
          }	
        }
    }//split tree
  }


  info->SetEntries(nentries);
  info->SetCounters(cal->GetICHCtr(),cal->GetHodoCtr(),cal->GetCARD29Ctr(),cal->GetGRETACtr(),cal->GetSCINTCtr());
  info->SetEff(cal->GetICHCtr(),cal->GetOBJCtr(),cal->GetXFPCtr(),cal->GetTOFCtr(),cal->GetPADCtr(),
	cal->GetTRACKCtr(),cal->GetPPACCtr(),cal->GetIITRACKCtr(),cal->GetCARD29Ctr(),cal->GetSCINTCtr());
  cout << endl;
  
  Long64_t filesize =0;
  if(CutFile==NULL){
    ctr->Write("",TObject::kOverwrite);
    filesize = ctr->GetZipBytes();
  }
  else{
    for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
      for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
        splittree[in][ou]->Write("",TObject::kOverwrite);
        filesize += splittree[in][ou]->GetZipBytes();
      }
    }
  }


  info->Write("runinfo",TObject::kOverwrite);
  cout <<  endl;
  cout << "closing file ..." << endl;

  cout<<"Size of input tree  "<<setw(7)<<gtr->GetZipBytes()/(1024*1024)<<" MB"<<endl
      <<"size of calibrated tree(s) "<<setw(7)<<filesize/(1024*1024)<<" MB"<<endl
      <<"=> size of calibrated tree(s) is "<<setw(5)<<setiosflags(ios::fixed)<<setprecision(1)<<(100.*filesize)/gtr->GetZipBytes()<<" % of size of input tree"<<endl;

  infile->Close();
  ofile->Close();
  double time_end = get_time();
  cout << "Program Run time " << time_end - time_start << " s." << endl;
  cout << "Calculated " << nentries/(time_end - time_start) << " events/s." << endl;
  timer.Stop();
  cout << "\n CPU time: " << timer.CpuTime() << "\tReal time: " << timer.RealTime() << endl;

  return 0;
}
