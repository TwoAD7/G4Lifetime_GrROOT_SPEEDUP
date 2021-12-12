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
#include<thread>
#include<mutex>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TROOT.h"
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

//Create our shared lock
mutex mtx;
S800Calc* s800Calc; 
GretinaCalc* gretinaCalc; 
GretinaEvent* gretinaEvent; 
Mode3Calc* mode3Calc;
double time_start=0;
char* CutFile = NULL;
vector<TCutG*> InPartCut;
vector<vector<TCutG*> > OutPartCut;
vector<vector<TTree*>> splittree;
int tac =0;


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


//Function to specifically fill the Tree to assure there is no data race among the threads
void clearing(Mode3Event *m3e, S800 *s800, Gretina *gretina){
  mtx.lock();
  /*
  s800.Clear();
  gretina.Clear();
  m3e.Clear(); 
  */
  s800->Clear();
  gretina->Clear();
  m3e->Clear();
  mtx.unlock();
}

// void fill_tree(TTree *ctr, S800Calc *s800Calc, GretinaCalc *gretinaCalc, 
// GretinaEvent *gretinaEvent, Mode3Calc *mode3Calc, vector<vector<TTree*>> &splittree_i,
// TTree *gtr, Mode3Event *m3e, S800 *s800, Gretina *gretina, 
//  int n_itr_per_thrd, int ID, Int_t &nbytes){

void fill_tree(TTree *ctr, S800Calc *s800Calc, GretinaCalc *gretinaCalc, 
GretinaEvent *gretinaEvent, Mode3Calc *mode3Calc,
TTree *gtr, Mode3Event *m3e, S800 *s800, Gretina *gretina, 
 int n_itr_per_thrd, int ID, Int_t &nbytes){

  //printf("INSIDE THREAD FUNCTION, THD %d\n", ID);

  bool id_0 = false;
  if(ID == 0)
    id_0 = true;

  for(int i = ID*n_itr_per_thrd ;i< n_itr_per_thrd + n_itr_per_thrd*ID;i++){

    int status = 0;
    //nentries = nent;


    //Commented this out because can't have break statements in OMP
    /*
    if(signal_received)
      break;
    */

    /*
    if(set->VLevel()>2){
      cout << "-----------------------------------------"<< endl;
      cout << "processing event " << i << endl;
    }
    */

    if(id_0){
      //printf("INSIDE THE FIRST IF STATEMENT\n");
      cout<<"On entry: "<<i<<endl;
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/n_itr_per_thrd<<" % done\t" << 
      (Float_t)i/(time_end - time_start) << " events/s " << (n_itr_per_thrd-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
      cout<<endl;
      
      //cal->PrintCtrs();
      /*
      double time_end = get_time();
      cout << setw(5) << setiosflags(ios::fixed) << setprecision(1) << (100.*i)/nentries<<" % done\t" << 
      (Float_t)i/(time_end - time_start) << " events/s " << (nentries-i)*(time_end - time_start)/(Float_t)i << "s to go \r" << flush;
      */
    }

      
    if(i%1000000 == 0 && CutFile==NULL){
      mtx.lock();
      //printf("ATTEMPTING TO AUTOSAVE. THREAD: %d with run %d\n",ID,i);
      //ctr.AutoSave();
      ctr->AutoSave();
      mtx.unlock();
      //ofile->cd();
      //ctr_objs.at(ID)->AutoSave();
    }
    

    /*
    s800Calc.Clear();
    gretinaCalc.Clear();
    gretinaEvent.Clear();
    mode3Calc.Clear();    
    */
    
    //It clears all of the data entries in that object
    //s800Calc.Clear();
    //gretinaCalc.Clear();
    //gretinaEvent.Clear();
    //mode3Calc.Clear();
      
    //mtx.lock();
    //printf("CLEARING THE CALC\n");
    s800Calc->Clear();
    gretinaCalc->Clear();
    gretinaEvent->Clear();
    mode3Calc->Clear();
    //mtx.unlock();



    mtx.lock();
    //printf("CLEARING\n");
    s800->Clear();
    gretina->Clear();
    m3e->Clear(); 
    status = gtr->GetEntry(i,1);
    mtx.unlock();
    

    //This section requires a shared lock
    //clearing(m3e,s800,gretina);
    //mtx.lock();
    //printf("ATTEMPTING TO GET ENTRY, THREAD %d\n",ID);
    
    /*
    if(gtr == NULL)
      printf("THREAD %d, GTR IS NULL\n",ID);
    if(gtr!= NULL){
          TObjArray *cp = (TObjArray*)gtr->GetListOfBranches();
        for(int j =0;j<cp->GetEntries();++j)
          cout<<cp->At(j)->GetName()<<endl;
      }
    cout<<"VAL: "<<i<<" "<<ID<<endl;
    */

    //status = gtr->GetEntry(i,1);
    //status = (int)gtr->GetEvent(i);
    //printf("THE STATUS IS: %d\n", status);
    //status = (int)gtr.GetEntry(i);
    //mtx.unlock();

    if(status == -1){
      //cerr<<"Error occured, couldn't read entry "<<i<<" from tree "<<gtr->GetName()<<" in file "<<infile->GetName()<<endl;
      continue;
      //return 5;
    }
    else if(status == 0){
      //cerr<<"Error occured, entry "<<i<<" in tree "<<gtr->GetName()<<" in file "<<infile->GetName()<<" doesn't exist"<<endl;
      continue;
      //return 6;
    }

    nbytes += status;

    /*
    if(trackMe){
      cal->GammaTrack(gretinaCalc,gretinaEvent)
    }
    */

    if(CutFile==NULL){
      //ctr.Fill();
      //mtx.lock();
      //printf("CUTFILE IS NULL\n");
      ctr->Fill();
      //mtx.unlock();
    }

    
    else{
      cout<<"IN THE ELSE FOR CUTFILE\n";
      /*
      mtx.lock();
      printf("IN THE ELSE STATEMENT\n");
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
          splittree_i[in][ou]->Fill();
          } 
        }
      mtx.unlock();
      */
    }//split tree
  


  }
}



int main(int argc, char* argv[]){
  //double time_start = get_time();  
  time_start = get_time();
  TStopwatch timer;
  timer.Start();
  signal(SIGINT,signalhandler);

  char *InputFile = NULL;
  char *RootFile = NULL;
  char* SettingFile;
  //char* CutFile = NULL;
  int LastEvent = -1;
  int addbackType = -1;
  bool trackMe = false;
  //int tac = 0;

  vector<thread> thds;

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

  cout<<"OPENING OUTPUT FILE\n";
  TFile *ofile = new TFile(RootFile,"RECREATE");
  cout<<"OPENED OUTPUT FILE\n";



  //gtr is the input TTree
  TFile* infile = new TFile(InputFile);
  TTree* gtr = (TTree*) infile->Get("gtr");
  if(gtr == NULL){
    cout << "could not find tree tr in file " << infile->GetName() << endl;
    return 3;
  }

  int nentries = (int)gtr->GetEntries();
  //int max_threads = omp_get_max_threads();

  int max_threads = 8;
  //int itr_per_thread =  nentries/max_threads + (nentries % max_threads >0);
  int itr_per_thread =  nentries/max_threads + (nentries % max_threads >0);
  printf("ITERATION PER THREAD %d",itr_per_thread);

  // The serial version
  Mode3Event* m3e = new Mode3Event;
  gtr->SetBranchAddress("mode3Event",&m3e);
  S800 *s800 = new S800;
  gtr->SetBranchAddress("s800",&s800);
  Gretina* gretina = new Gretina;
  gtr->SetBranchAddress("gretina",&gretina);

  if(gtr!= NULL){
          TObjArray *cp = (TObjArray*)gtr->GetListOfBranches();
        for(int j =0;j<cp->GetEntries();++j)
          cout<<cp->At(j)->GetName()<<endl;
      }
  

  //TFile *ofile = new TFile(RootFile,"RECREATE");
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
  
  printf("SETTING FILE: %s\n",SettingFile);
  set->SetTracking(trackMe);
  //printf("IS IT HERE\n");
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

  //vector<TCutG*> InPartCut;
  //vector<vector<TCutG*> > OutPartCut;
  //vector<vector<TTree*>> splittree;
  vector<vector<vector<TTree*> > > splittree_vec;

  // ============= CREATING THE TTREE FOR EACH THREAD ============= //
  std::vector<TTree*> ctr_objs;
  std::vector<S800Calc*> S800Calc_objs;
  std::vector<GretinaCalc*> gretinaCalc_objs;
  std::vector<GretinaEvent*> gretinaEvent_objs;
  std::vector<Mode3Calc*> mode3Calc_objs;

  //Generatate all the data objects that are going to be branches in each TTree
    for(int i=0;i<max_threads;i++){
      S800Calc_objs.emplace_back(new S800Calc);
      gretinaCalc_objs.emplace_back(new GretinaCalc);
      gretinaEvent_objs.emplace_back(new GretinaEvent);
      mode3Calc_objs.emplace_back(new Mode3Calc);
    }
        
    // Create vector containing TTrees. We make a TTree for each thread with their 
    // own branches.
    for(int i =0;i< max_threads;i++){
      printf("CREATING THE TTREES %d\n",i);
      TTree *ctr = new TTree("ctr","Gretina/S800 calibrated and builtevents");
      //ctr_objs.emplace_back(new TTree(Form("ctr%d",i),"Gretina/S800 calibrated and builtevents"));
      //TTree * ctr =  ctr_objs.at(i);

      /*
      S800Calc* s800Calc = S800Calc_objs.at(i);
      GretinaCalc* gretinaCalc = gretinaCalc_objs.at(i);
      GretinaEvent* gretinaEvent = gretinaEvent_objs.at(i);
      Mode3Calc* mode3Calc = mode3Calc_objs.at(i);
      */

      //Set different memory addresses for each branch in each different TTree
      /*
      ctr_objs.at(i)->Branch("s800calc",&s800Calc, 320000);
      ctr_objs.at(i)->Branch("gretinacalc",&gretinaCalc, 320000);
      if(trackMe)
        ctr_objs.at(i)->Branch("gretinaevent",&gretinaEvent, 320000);
      ctr_objs.at(i)->Branch("mode3calc",&mode3Calc, 320000);
      ctr_objs.at(i)->BranchRef();
      */


      ctr->Branch("s800calc",S800Calc_objs.at(i), 320000);
      ctr->Branch("gretinacalc",gretinaCalc_objs.at(i), 320000);
      if(trackMe)
        ctr->Branch("gretinaevent",gretinaEvent_objs.at(i), 320000);
      ctr->Branch("mode3calc",mode3Calc_objs.at(i), 320000);
      ctr->BranchRef();

      ctr_objs.emplace_back(ctr);
    }



  if(CutFile!=NULL){
    cout<<"HEY CUTFILE IS NULL\n";
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

  //   //vector<TTree*> trees;// = new vector<TTree*>;
  // vector<thread> thds;





    /*
    // Serial version 
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
    */

    for(int i =0;i< max_threads;i++){
      vector<vector<TTree*> > splittree_i = splittree_vec[i];
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
  

  //Double_t nentries = gtr->GetEntries();
  int nbytes=0;
  //Int_t nbytes = 0;
  //Int_t status;

  cout << "Nr of Events " << nentries << endl;
  if (LastEvent!=-1){
    nentries = LastEvent;
  }
   else if(set->LastEvent()>0){
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


  printf("MAX NUMBER OF THREADS: %d\n",max_threads);

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  // FUNCTION PLACE HOLDER FOR THE 'FOR LOOP'

/*

  TTree* ctr_test = new TTree("ctr_test","Gretina/S800 calibrated and builtevents");
  S800Calc* s800Calc = new S800Calc;
  GretinaCalc* gretinaCalc = new GretinaCalc;
  GretinaEvent* gretinaEvent = new GretinaEvent;
  Mode3Calc* mode3Calc = new Mode3Calc;
  ctr_test->Branch("s800calc",&s800Calc, 320000);
  ctr_test->Branch("gretinacalc",&gretinaCalc, 320000);
  if(trackMe)
    ctr_test->Branch("gretinaevent",&gretinaEvent, 320000);
  ctr_test->Branch("mode3calc",&mode3Calc, 320000);
  ctr_test->BranchRef();

*/

  vector<Int_t> nbytes_vec;
  for(int i =0;i<max_threads;i++)
    nbytes_vec.emplace_back(0);

  double time_s = omp_get_wtime();
  // FOR LOOP FOR ALL THE THREADS
  for(int i=0;i<max_threads;i++){
    // pass in the function to fill the tree followed by the input parameters of 
    // the TTree at index i (for the ith thread), the input tree (shared among all trees),
    // and the number of iterations per thread. 
    // thds.emplace_back(fill_tree,std::ref(ctr_objs[i]), std::ref(S800Calc_objs[i]),std::ref(gretinaCalc_objs[i]),std,
    //   std::ref(gretinaEvent_objs[i]), std::ref(mode3Calc_objs[i]),
    //   std::ref(gtr),std::ref(s800), std::ref(gretina),std::ref(m3e),
    //   itr_per_thread,i);

      TTree *ctr =  ctr_objs[i];
      S800Calc *s800Calc = S800Calc_objs[i];
      GretinaCalc *gretinaCalc = gretinaCalc_objs[i];
      GretinaEvent *gretinaEvent = gretinaEvent_objs[i];
      Mode3Calc *mode3Calc = mode3Calc_objs[i];
      Int_t nbytes = (Int_t)nbytes_vec[i];
      
      //vector<vector<TTree*>> st;
      //if(CutFile == NULL)
      //  vector<vector<TTree*>> st = splittree_vec[i];

      /*

      // This is for testing purposes
      if(ctr == NULL)
        printf("THREAD %d, CTR IS NULL\n",i);
      if(ctr!= NULL){
          printf("CTR IS NOT NULL. PRINTING CONTENTS.\n");
            TObjArray *cp = (TObjArray*)ctr->GetListOfBranches();
          for(int j =0;j<cp->GetEntries();++j)
            cout<<cp->At(j)->GetName()<<endl;
        }
      */
    


      //cout<<"STARTING A THREAD\n";
      /*
      thds.emplace_back(fill_tree,ctr_objs[i], S800Calc_objs[i],gretinaCalc_objs[i],
      gretinaEvent_objs[i], mode3Calc_objs[i], std::ref(splittree[i]),
      gtr,m3e,s800,gretina,
      itr_per_thread,i, std::ref(nbytes_vec[i]));
      */

       // thds.emplace_back(fill_tree,ctr, s800Calc,gretinaCalc,
       // gretinaEvent, mode3Calc, std::ref(st),
       // gtr,m3e,s800,gretina,
       // itr_per_thread,i, std::ref(nbytes));

        thds.emplace_back(fill_tree,ctr, s800Calc,gretinaCalc,
        gretinaEvent, mode3Calc,
        gtr,m3e,s800,gretina,
        itr_per_thread,i, std::ref(nbytes));

  }


    // Joining all threads
    int id=0;
    for(auto &t : thds){
      if( t.joinable() ){
        cout<<"JOINING THREADS\n";
        t.join();
      }
    }
    cout<<"DONE JOINING THREADS\n";
    double time_e = omp_get_wtime();
    cout<<"THREADING TOOK "<<time_e-time_s<<" seconds."<<endl;


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  // To merge the individual TTrees into one 
 // TList *cal_list = new TList;
  TList *ctr_list = new TList;
  //int size = cal_objs.size();
  int size = ctr_objs.size();
  for(int i =0;i<size;i++){
    //cal_list->Add(cal_objs.at(i));
    ctr_list->Add(ctr_objs.at(i)); 
  }

  TTree *ctr =  TTree::MergeTrees(ctr_list);

  info->SetEntries(nentries);
  printf("PAST THE PARALLELIZATION AND SET NENTRIES");
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
        //splittree[in][ou]->Write("",TObject::kOverwrite);
        //filesize += splittree[in][ou]->GetZipBytes();
        splittree_vec.at(0)[in][ou]->Write("",TObject::kOverwrite);
        filesize += splittree_vec.at(0)[in][ou]->GetZipBytes();
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
