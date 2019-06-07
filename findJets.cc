///
/// modified from: 01-basic.cc
#include <fstream>
#include <vector>
#include <iterator>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstdio>   // needed for io
#include <ctime>
#include "fastjet/ClusterSequence.hh"
//--- ROOT includes
#include <TApplication.h>
#include <TTree.h>
#include <TChain.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>

using namespace std;

void outputJetTotalsWfracs(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets);
void outputFullJets(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets);
void outputFullJets_WTA(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets, double R);

void processEvents(TString dirname, int color);

TCanvas *cpt, *cRho;

int main(int argc, char* argv[]) {
 TApplication theApp("App", &argc, argv);
 cpt = new TCanvas("pT","pT") ;
 cRho = new TCanvas("jet structure","jet structure") ;
 time_t time0 = time(nullptr);
 processEvents("/dlocal/eposvhlle.out/b0b34_trig50_vac7/events/", kBlack);
 processEvents("/dlocal/eposvhlle.out/b0b34_trig50_hybM_NOrec/events/", kGreen);
 processEvents("/dlocal/eposvhlle.out/b0b34_trig50_hybM_rec/events/", kRed);
 time_t time1 = time(nullptr);
 cout << "walltime: " << difftime(time1, time0) << " s.\n";
 //----------------------------------------------------------
 //for(int iEvent=0; iEvent<nEvents+1; iEvent++) { // find jets in each event
  //vector<fastjet::PseudoJet> input_particles;
  //input_particles.clear();

    //// input_particles.push_back(fastjet::PseudoJet(px,py,pz,E));
  //for(int i=0; i<EIn[iEvent].size(); i++) {
   ////double pt = sqrt(pow(pxIn[iEvent].at(i),2) + pow(pyIn[iEvent].at(i),2));
   ////if(pt > 0.5)  // no difference for vacuum jets!
   ////double pmod = sqrt(pow(pxIn[iEvent].at(i),2) + pow(pyIn[iEvent].at(i),2) + 
   ////  pow(pzIn[iEvent].at(i),2));
   ////double rap = 0.5*log((pmod+pzIn[iEvent].at(i))/(pmod-pzIn[iEvent].at(i)));
   ////if(fabs(rap)<1.5)  // no unexpected difference for vacuum jets!
   //fastjet::PseudoJet particle (pxIn[iEvent].at(i),
     //pyIn[iEvent].at(i), pzIn[iEvent].at(i), EIn[iEvent].at(i));
   //particle.set_user_index(jetIdIn[iEvent].at(i));
   //input_particles.push_back(particle);
  //}
  //// create a jet definition:
  //// a jet algorithm with a given radius parameter
  //double R = 0.4;
  //fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  //// run the jet clustering with the above jet definition
  //fastjet::ClusterSequence clust_seq(input_particles, jet_def);

  //// get the resulting jets ordered in pt
  //double ptmin = 3.0;
  //vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
 
  //outputJetTotalsWfracs(fout, iEvent, inclusive_jets);
  ////outputFullJets(fout, iEvent, inclusive_jets);
  ////outputFullJets_WTA(fout, iEvent, inclusive_jets, R);

 //} // end event loop
  theApp.Run();
  return 0;
}


void processEvents(TString dirname, int color)
{
 static int iPlotSeq = 0;
 // loadnig the tree
 TChain *tree = new TChain("treefin;1");
 TSystemDirectory dir ("rootfiles",dirname);
 TList *files = dir.GetListOfFiles();
 TSystemFile *file;
 TIter next(files);
 while ((file=(TSystemFile*)next())) {
  TString fname = file->GetName();
  cout << "[" << fname.Data() << "]";
  if(strstr(fname,".root"))
   tree->Add(dirname+fname);
 }
 cout << endl;
 // accessing the events
 const int NP = 50000 ;
 Int_t id[NP], status[NP] ;
 Float_t px [NP], py[NP], pz[NP], E[NP] ;
 Short_t ele[NP] ;
 Int_t npart ;
 int nevents = tree->GetEntries() ;
 tree->SetBranchAddress("px",&px[0]) ;
 tree->SetBranchAddress("py",&py[0]) ;
 tree->SetBranchAddress("pz",&pz[0]) ;
 tree->SetBranchAddress("E",&E[0]) ;
 tree->SetBranchAddress("id",&id[0]) ;
 tree->SetBranchAddress("status",&status[0]) ;
 tree->SetBranchAddress("ele",&ele[0]) ;
 tree->SetBranchAddress("npart",&npart) ;
 cout<<"processing, events = "<<nevents<<endl ;
 TH1F *hpt = new TH1F(("hpt"+to_string(iPlotSeq)).c_str(),
    ("hpt"+to_string(iPlotSeq)).c_str(), 20, 0., 200.);
 TH1F *hptJet = new TH1F(("hptJet"+to_string(iPlotSeq)).c_str(),
    ("hptJet"+to_string(iPlotSeq)).c_str(), 20, 0., 200.);
 TH1F *hRho = new TH1F(("rhoJet"+to_string(iPlotSeq)).c_str(),
    ("rhoJet"+to_string(iPlotSeq)).c_str(), 20, 0., 0.6);
  TH1F *hRhoMed = new TH1F(("rhoMed"+to_string(iPlotSeq)).c_str(),
    ("rhoMed"+to_string(iPlotSeq)).c_str(), 20, 0., 0.6);
 for(int iev=0; iev<nevents; iev++) { // event loop
  tree->GetEntry(iev) ;
  vector<fastjet::PseudoJet> input_particles;
  input_particles.clear();
  for(int ip=0; ip<npart; ip++) { // particle loop
   double rap = 0.5*log((E[ip]+pz[ip])/(E[ip]-pz[ip]));
   double pt = sqrt(px[ip]*px[ip]+py[ip]*py[ip]);
   if(status[ip]==10 && fabs(rap)<1.0) {
    hpt->Fill(pt);
   }
   // constructing Fastjet input
   if(status[ip]==1 && pt>0.5 && fabs(rap)<1.0) {
   fastjet::PseudoJet particle (px[ip], py[ip], pz[ip], E[ip]);
   particle.set_user_index(status[ip]);
   input_particles.push_back(particle);
   }
  } // end of particle loop
  // jet finding
  double R = 0.3;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);
  // get the resulting jets ordered in pt
  double ptmin = 3.0;
  vector<fastjet::PseudoJet> jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  for(unsigned int i = 0; i < jets.size(); i++) {
   hptJet->Fill(jets[i].perp());
   // jet structure
   if(jets[i].perp()>40.0) {
    vector<fastjet::PseudoJet> constituents = jets[i].constituents();
    for(fastjet::PseudoJet ptl : constituents) {
     double phi = ptl.delta_R(jets[i]);
     hRho->Fill(phi, ptl.perp());
     if(ptl.user_index()==0) hRhoMed->Fill(phi, ptl.perp());
    }
   }
  }
 } // end of event loop
 // normalizing histos
 const char* drawOpt = iPlotSeq==0? "" : "same";
 cpt->cd();
 hpt->Scale(1./(hpt->GetBinWidth(1)*2. * nevents));
 hpt->Draw(drawOpt);
 hptJet->Scale(1./(hptJet->GetBinWidth(1)*2. * nevents));
 hptJet->Draw("same");
 cRho->cd();
 hRho->Scale(1./(nevents));
 hRho->SetLineColor(color);
 hRho->Draw(drawOpt);
 hRhoMed->Scale(1./(nevents));
 hRhoMed->SetLineColor(color);
 hRhoMed->SetMarkerStyle(22);
 hRhoMed->Draw("same");
 iPlotSeq++;
}


void outputJetTotalsWfracs(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets)
// prints global jet properties (energy-momentum and fractions) to a file
{
 for (unsigned int i = 0; i < jets.size(); i++) {
  vector<fastjet::PseudoJet> constituents = jets[i].constituents();
  map<int,double> origins; // map containing the indexes of jets the partons are coming from
  bool leadingTrigg = false;
  for(fastjet::PseudoJet ptl : constituents) {
   if(ptl.perp() > 5.0) leadingTrigg = true;
   int origJet = ptl.user_index(); // get back the index of original jet
   if(origins.count(origJet)==0)
    origins[origJet] = sqrt(ptl.px()*ptl.px()+ptl.py()*ptl.py());
   else
    origins[origJet] += sqrt(ptl.px()*ptl.px()+ptl.py()*ptl.py());//ptl.perp();
  }
  vector <double> listFractions;
  listFractions.reserve(10);
  for(auto it = origins.begin(); it != origins.end(); ++it) {
   listFractions.push_back(it->second);
  }
  while(listFractions.size()<4)
   listFractions.push_back(0);
  sort(listFractions.begin(), listFractions.end(), std::greater<double>());
  //if(leadingTrigg)
  fout << setw(8) << iEvent << setw(14) << jets[i].E()
     << setw(14) << jets[i].px() << setw(14) << jets[i].py()
     << setw(14) << jets[i].pz()
     << setw(14) << listFractions[0] << setw(14) << listFractions[1]
     << setw(14) << listFractions[2] << setw(14) << listFractions[3] << endl;
 }
}


void outputFullJets(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets)
// output for jet structure plots
{
 for (unsigned int i = 0; i < jets.size(); i++) {
  vector<fastjet::PseudoJet> constituents = jets[i].constituents();
  bool leadingTrigg = false;
  for(fastjet::PseudoJet ptl : constituents)
   if(ptl.perp() > 5.0) leadingTrigg = true;
  //if(leadingTrigg)
  for(fastjet::PseudoJet ptl : constituents) {
   double phi = ptl.delta_R(jets[i]);
   fout << setw(8) << iEvent << setw(8) << i
     << setw(14) << jets[i].perp() << setw(14) << jets[i].rap()
     << setw(14) << ptl.perp() << setw(14) << phi << endl;
   }
 }
}


void outputFullJets_WTA(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets, double R)
// output for jet structure plots
{
 for (unsigned int i = 0; i < jets.size(); i++) {
  vector<fastjet::PseudoJet> constituents = jets[i].constituents();
  bool leadingTrigg = false;
  vector<fastjet::PseudoJet> input_particles; // to be filled with constituents for the 2nd jet finding
  input_particles.clear();
  for(fastjet::PseudoJet ptl : constituents) {
   if(ptl.perp() > 5.0) leadingTrigg = true; // hard z trigger
   input_particles.push_back(ptl);
  }
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R, fastjet::WTA_pt_scheme);
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);
  double ptmin = 3.0;
  vector<fastjet::PseudoJet> foundJets2 = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  //cout << "-------\n";
  //for(fastjet::PseudoJet jets2 : foundJets2) {
   //cout << "FoundJets2: " << jets2.perp() << endl;
  //}
  //if(leadingTrigg)
  for(fastjet::PseudoJet ptl : constituents) {
   double phi = ptl.delta_R(foundJets2[0]);
   fout << setw(8) << iEvent << setw(8) << i
     << setw(14) << jets[i].perp() << setw(14) << jets[i].rap()
     << setw(14) << ptl.perp() << setw(14) << phi << endl;
   }
 }
}
