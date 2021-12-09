//Testing parallelization of filling a histogram 
#include<omp.h>
#include<iostream>
#include "TRandom.h"
#include "TH1.h"
#include<vector>
#include "TFile.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TRootCanvas.h"
#include<thread>
//using std::thread;

/*
Helpful link:

https://root.cern/manual/creating_a_user_application/
https://root-forum.cern.ch/t/compiling-a-c-macro-using-root-classes-with-g/7515
http://wwwinfo.jinr.ru/~tanyusha/OpenMp_en.html

To run this on one of the gateways on the nscl/frib server (after loading the proper containers 
and modules):

g++ -std=c++11 -fopenmp histo_test.cpp `root-config --cflags --glibs`

*/


// Function to merge all histograms. This is used by 
// the parallel implementations 
TH1D *merge_histograms(TString name){
//TH1D *merge_histograms(std::vector<TH1D*> vec){
	int l = omp_get_max_threads();
	TList *list =  new TList();
	TFile *f = new TFile(name,"read");
	TH1D *tot = (TH1D*)f->Get("Hist0");
	for(int i =1;i<l;i++)
		tot->Add((TH1D*)f->Get(Form("Hist%d",i)));

	return tot; 
}


//Function to use c++ threads
void fill_hist(std::vector<TH1D*> &hists, int thd_id){
	TRandom *rng = new TRandom();
	double myRand;
	for (int j = 0; j < 10000000; j++) {
		myRand = rng->Gaus(350+thd_id*8,20+2*thd_id);
		hists.at(thd_id)->Fill(myRand);
	}
}

int main(int argc,char **argv){


	TH1D* hist_serial = new TH1D("hist_serial","Histogram",1000,0,1000);

	//Uncomment if using the omp method
	//TRandom *rng = new TRandom();
	//double myRand;

	//omp_set_num_threads(4);
	int nhist = omp_get_max_threads(); //one histogram per thread
	//int thds = omp_get_max_threads();
	printf("Using %d threads\n",nhist);

	TString ofname = "omp_histo_test.root";
	TFile *f = new TFile(ofname,"recreate");

	std::vector<TH1D*> hists;
	//TH1D *hists[nhist];

	//Create vector containing histograms
	for(int i =0;i< nhist;i++){
		//hists[i] = new TH1D(Form("Hist%d",i),"Histogram",1000,0,1000);
		hists.emplace_back(new TH1D(Form("Hist%d",i),"Histogram",1000,0,1000));
	}

	double t1 = omp_get_wtime();

	
	//TH1D* hist_test = new TH1D("hist_test","Histogram",1000,800,600);
	std::vector<std::thread> thrds;

	for(int i=0;i<nhist;i++){
		thrds.emplace_back(fill_hist, std::ref(hists), i);
	}

	//Join all of the threads
	for(auto &t : thrds){
		if( t.joinable() )
			t.join();
	}
	

	/*
	//int i,j;
	//# pragma omp parallel for num_threads(nhist) shared(hists) shared(hists) schedule(guided) collapse(2)
	//for (int i = 0; i < nhist; i++) { schedule(dynamic,500) shared(hists) collapse(2) reduction(+:hist[:1000])
	//nhist = 6;
	//double hist[1000];
	//double hist[nhist];
	//nhist = 6;
	# pragma omp parallel for num_threads(nhist) shared(hists) schedule(dynamic,1000)
	for(int i=0; i<nhist;i++){
		//printf("On thread %d with loop iteration %d\n",omp_get_thread_num(),i);
	for (int j = 0; j < 10000000; j++) {
		//printf("Thread ID: %d\n",i);
		//int i = omp_get_thread_num();
			 //myRand = rng->Gaus(i,i);
			 
			 myRand = rng->Gaus(350+i*8,20+2*i);
			 //hist[i]+=rng->Gaus(350+i*8,20+2*i);

			 //myRand = rng->Gaus(j,j);
			 //#pragma omp critical
			 //hist_test->Fill(myRand);
			 
			 hists.at(i)->Fill(myRand);
			 
			 //hist_test->Fill(myRand);
			 //hists[i]->Fill(myRand);
			 //std::cout<<bin<<std::endl;
	}
}
	
	//TH1D *total = new TH1D("total","Total",1000,0.0,1000.);
	//total->FillN(1000,hist,1);
	*/
	



	double t2 = omp_get_wtime();
	std::cout << "Elapsed: " << t2-t1 << std::endl;

	


	f->Write();
	f->Close();
	
	
	
	/*
	// Serial version of filling a histogram
	for (int i = 0; i < nhist; i++) {
		for (int j = 0; j < 10000000; j++) {
			 //myRand = rng->Gaus(i,i);
			 myRand = rng->Gaus(350+i*8,20+2*i);
			 hist_serial->Fill(myRand);
		}
	}
	
	
	
	double t2 = omp_get_wtime();
	std::cout << "Elapsed: " << t2-t1 << std::endl;
	
	*/

	
	TH1D *total = merge_histograms(ofname);
	
	//TH1D *total = merge_histograms(hists);

	
	TApplication app("app",&argc,argv);
	TCanvas* c = new TCanvas("c", "Something", 0, 0, 800, 600);
	//hist->GetXaxis()->SetRangeUser(0,800);
	//printf("size of vector: %d\n",hists.size());
	//TH1D *t = hists[0];
	//t->Draw();
	//hist_test->Draw();
	total->Draw();
	c->Modified(); 
	c->Update();
	TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
	
	

	
	/*
	TApplication app("app",&argc,argv);
	TCanvas* c = new TCanvas("c", "Something", 0, 0, 800, 600);
	hist_serial->Draw();
	c->Modified(); 
	c->Update();
	TRootCanvas *rc = (TRootCanvas *)c->GetCanvasImp();
	rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");
	app.Run();
   
	*/

	return 0;

}