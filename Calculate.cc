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




  //gtr is the input TTree
  TFile* infile = new TFile(InputFile);
  TTree* gtr = (TTree*) infile->Get("gtr");
  if(gtr == NULL){
    cout << "could not find tree tr in file " << infile->GetName() << endl;
    return 3;
  }

  /*
  // The serial version
  Mode3Event* m3e = new Mode3Event;
  gtr->SetBranchAddress("mode3Event",&m3e);
  S800 *s800 = new S800;
  gtr->SetBranchAddress("s800",&s800);
  Gretina* gretina = new Gretina;
  gtr->SetBranchAddress("gretina",&gretina);
  */  


  // The parallel version

  std::vector<TTree*> gtr_objs; 
  std::vector<Gretina*> gretina_objs;
  std::vector<S800*> s800_objs;
  std::vector<Mode3Event*> m3e_objs;
  int thd_cnt =  omp_get_max_threads();

  //Clone the input TTree into each one
  for(int i=0;i<thd_cnt;i++){
    gtr_objs.emplace_back((TTree*)gtr->CloneTree());
  }


  //Generatate all the data objects that are going to be branches in each TTree
  /*
  Note:
  The reason we created one input TTree for each thread is due to the "SetBranchAddress" that occurs.
  If we only have one input TTree, everytime we set the BranchAddress, it is going to get overriden. 
  To overcome this, I made multiple copies of the input TTree and set the branch address of the ith 
  thread to the address of the ith input TTree. THis way, each thread has its own input TTree. 
  */
  for(int i=0;i<thd_cnt;i++){
    s800_objs.emplace_back(new S800);
    gretina_objs.emplace_back(new Gretina);
    m3e_objs.emplace_back(new Mode3Event);

    Gretina *gretina = gretina_objs.at(i);
    gtr_objs.at(i)->SetBranchAddress("gretina",&gretina);

    Mode3Event *m3e = m3e_objs.at(i);
    gtr_objs.at(i)->SetBranchAddress("mode3Event",&m3e);

    S800 *s800 = s800_objs.at(i);
    gtr_objs.at(i)->SetBranchAddress("s800",&s800);


  }
  



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

  /*
  Calibration *cal = new Calibration(set,startevent,RunNumber); //AR New in v4.1_TL
  cout << "Brho " << cal->GetBrho() << " mass " <<  cal->GetMass() << " Z " << cal->GetZ() << endl;
  info->SetBeam(cal->GetBrho(), cal->GetMass(), cal->GetZ());
  */

  // The parallel Version
  vector<Calibration*> cal_objs;
  for(int i =0;i< thd_cnt;i++)
    cal_objs.emplace_back(new Calibration(set,startevent,RunNumber));

  cout << "Brho " << cal_objs.at(0)->GetBrho() << " mass " <<  cal_objs.at(0)->GetMass() << " Z " << cal_objs.at(0)->GetZ() << endl;
  info->SetBeam(cal_objs.at(0)->GetBrho(), cal_objs.at(0)->GetMass(), cal_objs.at(0)->GetZ());


  /*
  // The serial version
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
  */

  
  // ============ Section generating separate TTrees for Parallelization Purposes ============ //
  std::vector<TTree*> ctr_objs;
  std::vector<S800Calc*> S800Calc_objs;
  std::vector<GretinaCalc*> gretinaCalc_objs;
  std::vector<GretinaEvent*> gretinaEvent_objs;
  std::vector<Mode3Calc*> mode3Calc_objs;
  //int thd_cnt =  omp_get_max_threads();

  //Generatate all the data objects that are going to be branches in each TTree
  for(int i=0;i<thd_cnt;i++){
    S800Calc_objs.emplace_back(new S800Calc);
    gretinaCalc_objs.emplace_back(new GretinaCalc);
    gretinaEvent_objs.emplace_back(new GretinaEvent);
    mode3Calc_objs.emplace_back(new Mode3Calc);
  }
  

    
  //Create vector containing TTrees. We make a TTree for each thread.
  for(int i =0;i< thd_cnt;i++){
    //hists[i] = new TH1D(Form("Hist%d",i),"Histogram",1000,0,1000);
    ctr_objs.emplace_back(new TTree(Form("ctr%d",i),"Gretina/S800 calibrated and builtevents"));
    S800Calc* s800Calc = S800Calc_objs.at(i);
    GretinaCalc* gretinaCalc = gretinaCalc_objs.at(i);
    GretinaEvent* gretinaEvent = gretinaEvent_objs.at(i);
    Mode3Calc* mode3Calc = mode3Calc_objs.at(i);

    //Set different memory addresses for each branch in each different TTree
    ctr_objs.at(i)->Branch("s800calc",&s800Calc, 320000);
    ctr_objs.at(i)->Branch("gretinacalc",&gretinaCalc, 320000);
    if(trackMe)
    ctr_objs.at(i)->Branch("gretinaevent",&gretinaEvent, 320000);
    ctr_objs.at(i)->Branch("mode3calc",&mode3Calc, 320000);
    ctr_objs.at(i)->BranchRef();
  }
  



  vector<TCutG*> InPartCut;
  vector<vector<TCutG*> > OutPartCut;
  //vector<vector<TTree*> > splittree;
  vector<vector<vector<TTree*> > > splittree;
  


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


    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    /*
    NOTE: Could time this and see how long this takes. If further improvements in timing are necessary, 
    this could also be parallelized. In the mean time, it is going to be run serially.
    */
    for(int i =0;i< thd_cnt;i++){
      splittree_i = splittree.at(i);
      S800Calc* s800Calc = S800Calc_objs.at(i);
      GretinaCalc* gretinaCalc = gretinaCalc_objs.at(i);
      GretinaEvent* gretinaEvent = gretinaEvent_objs.at(i);
      Mode3Calc* mode3Calc = mode3Calc_objs.at(i);
      splittree_i.resize(InPartCut.size());
      for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
        splittree_i[in].resize(OutPartCut[in].size());
        for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
        	splittree_i[in][ou] = new TTree(Form("ctr_%s",OutPartCut[in][ou]->GetName()),"Gretina/S800 calibrated and builtevents");
        	splittree_i[in][ou]->Branch("s800calc",&s800Calc, 320000);
        	splittree_i[in][ou]->Branch("gretinacalc",&gretinaCalc, 320000);
        	if(trackMe)
        	  splittree_i[in][ou]->Branch("gretinaevent",&gretinaEvent, 320000);
        	splittree_i[in][ou]->Branch("mode3calc",&mode3Calc, 320000);
  	      splittree_i[in][ou]->BranchRef();
        }
      }
    }
  }
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  /*
  //After creating the input the Trees, make a copy for each thread by pushing them into the vector
  // NOTE: The current issue with doing this is that all of the cloned TTrees would have 
  // the same branch address and would all write to the same memory location in the foor loop
  // creating a data race
  std::vector<TTree*> ctr_objs;
  int thd_cnt =  omp_get_max_threads();

  for(int i=0;i<thd_cnt;i++){
    ctr_objs.emplace_back((TTree*)ctr->CloneTree());
  }

*/

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// This is a for loop that goes over all the enties. THIS WOULD NEED TO BE PARALLELIZED

  int max_threads = omp_get_max_threads();
  printf("MAX NUMBER OF THREADS: %d\n",max_threads);

  #pragma omp parallel for num_threads(max_threads) private(status) reduction(+:nbytes)
  for(int i=0; i<nentries; i++){
    if(signal_received)
      break;

    if(set->VLevel()>2){
      cout << "-----------------------------------------"<< endl;
      cout << "processing event " << i << endl;
    }
    int ID = omp_get_thread_num();
    if(i%10000 == 0 && ID == 0){
      //cal->PrintCtrs();
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries<<" % done\t" << 
      (Float_t)i/(time_end - time_start) << " events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
    }
    if(i%1000000 == 0 && CutFile==NULL){
      //ctr->AutoSave();
      ctr_objs.at(ID)->AutoSave();
    }


    /*
    //It clears all of the data entries in that object
    s800Calc->Clear();
    gretinaCalc->Clear();
    gretinaEvent->Clear();
    mode3Calc->Clear();
    s800->Clear();
    gretina->Clear();
    m3e->Clear();
    */

    S800Calc_objs.at(ID)->Clear();
    gretinaCalc_objs.at(ID)->Clear();
    gretinaEvent_objs.at(ID)->Clear();
    mode3Calc_objs.at(ID)->Clear();
    s800_objs.at(ID)->Clear();
    gretina_objs.at(ID)->Clear();
    m3e_objs.at(ID)->Clear();



    //For the entry "i", grab the event that occured from the TTree gtr (input tree)
    /*
    
    In using omp parallel for, the threads get the "ith" iteration (i.e. each loop doesn't start from 0, it starts on the true value)

    */
    status = gtr->GetEvent(i);

    if(status == -1){
      cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<gtr->GetName()<<" in file "<<infile->GetName()<<endl;
      return 5;
    }
    else if(status == 0){
      cerr<<"Error occured, entry "<<i<<" in tree "<<gtr->GetName()<<" in file "<<infile->GetName()<<" doesn't exist"<<endl;
      return 6;
    }

    //There will be a data race with this variable
    nbytes += status;



    //Build all of the calibrated objects, using the calibration in cal.
    //Use the data from the first three parameters, output into the last three parameters.

    /*
    cal->BuildAllCalc(s800,gretina,m3e,
      s800Calc,gretinaCalc,mode3Calc);
    */
    cal_objs.at(ID)->BuildAllCalc(s800_objs.at(ID),gretina_objs.at(ID),m3e_objs.at(ID),
      s800Calc_objs.at(ID),gretinaCalc_objs.at(ID),mode3Calc_objs.at(ID))

    if(trackMe){
      // cal->GammaTrack(gretinaCalc,gretinaEvent)
      cal_objs.at(ID)->GammaTrack(gretinaCalc_objs.at(ID),gretinaEvent_objs.at(ID));
    }

    // ctr is a TTree that we are filling with the calibrated data. The TTree::Fill() fills all the branches 
    /*
    If one were to create a parallel loop as is, the if-statement that follows would result in a seg fault due to the data 
    race that would take place (i.e. all of the threads would be writing to the same memory location, screwing things up, creating the 
    data race). 

    To overcome the data race, a good idea to try is to make individual TTree objects for each thread, and then at the end of the for loop,
    merge all of the TTrees together. 

    */
    if(CutFile==NULL){
      //ctr->Fill();
      ctr_objs.at(ID)->Fill();
    }
    /*
    //Serial version
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
    */

    //Parallel Version
    else{
      for(UShort_t in=0;in<InPartCut.size();in++){ // loop over incoming cuts
        if(tac%2 == 1 && !InPartCut[in]->IsInside(s800Calc_objs.at(ID)->GetTOF()->GetOBJ(),s800Calc_objs.at(ID)->GetTOF()->GetXFP()))
          continue;
        if(tac%2 == 0 && !InPartCut[in]->IsInside(s800Calc_objs.at(ID)->GetTOF()->GetTACOBJ(),s800Calc_objs.at(ID)->GetTOF()->GetTACXFP()))
          continue;
        for(UShort_t ou=0;ou<OutPartCut[in].size();ou++){ // loop over PID cuts
          if(tac == 1 && !OutPartCut[in][ou]->IsInside(s800Calc_objs.at(ID)->GetTOF()->GetOBJC(),s800Calc_objs.at(ID)->GetIC()->GetDE()))
            continue;
          if(tac == 0 && !OutPartCut[in][ou]->IsInside(s800Calc_objs.at(ID)->GetTOF()->GetTACOBJC(),s800Calc_objs.at(ID)->GetIC()->GetDE()))
            continue;
          if(tac == 3 && !OutPartCut[in][ou]->IsInside(s800Calc_objs.at(ID)->GetTOF()->GetXFPC(),s800Calc_objs.at(ID)->GetIC()->GetDE()))
            continue;
          if(tac == 2 && !OutPartCut[in][ou]->IsInside(s800Calc_objs.at(ID)->GetTOF()->GetTACXFPC(),s800Calc_objs.at(ID)->GetIC()->GetDE()))
            continue;
          splittree.at(ID)[in][ou]->Fill();
          } 
        }
    }
  }

  // To merge the individual TTrees into one 
  TList *cal_list = new TList;
  TList *ctr_list = new TList;
  int size = cal_objs.size();
  for(int i =0;i<size;i++){
    cal_list->Add(cal_objs.at(i));
    ctr_list->Add(ctr_objs.at(i)); 
  }

  TTree *ctr =  TTree::MergeTrees(ctr_list);

  info->SetEntries(nentries);
  print("PAST THE PARALLELIZATION AND SET NENTRIES")
  /*
  info->SetCounters(cal->GetICHCtr(),cal->GetHodoCtr(),cal->GetCARD29Ctr(),cal->GetGRETACtr(),cal->GetSCINTCtr());
  info->SetEff(cal->GetICHCtr(),cal->GetOBJCtr(),cal->GetXFPCtr(),cal->GetTOFCtr(),cal->GetPADCtr(),
	cal->GetTRACKCtr(),cal->GetPPACCtr(),cal->GetIITRACKCtr(),cal->GetCARD29Ctr(),cal->GetSCINTCtr());
  cout << endl;
  */
  
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
