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
#include "fastjet/ClusterSequence.hh"

using namespace std;


int main(int argc, char* argv[]) {

 if(argc != 3) {
  cout << "usage: ./findJets  in_file  out_file\n";
  return 1;
 }
 ifstream fin(argv[1]);
 if(!fin) {
  cout << "cannot open input file: " << argv[1] << endl;
  return 1;
 }
 ofstream fout(argv[2]);
 if(!fout) {
  cout << "cannot open output file: " << argv[2] << endl;
  return 1;
 }
 // read in the events
 int evIdCurr = 2147483646;  // a huge number not equal to any real event ID
 int nEvents = -1;
 vector<vector<int> > typeIn, colIn, acolIn, jetIdIn;
 vector<vector<double> > EIn, pxIn, pyIn, pzIn, timeIn;
 const int __Nevents = 200;
 jetIdIn.resize(__Nevents);
 typeIn.resize(__Nevents);
 colIn.resize(__Nevents);
 acolIn.resize(__Nevents);
 EIn.resize(__Nevents);
 pxIn.resize(__Nevents);
 pyIn.resize(__Nevents);
 pzIn.resize(__Nevents);
 timeIn.resize(__Nevents);
 string line;
 while (getline(fin, line)) {  // reading the initial partons, line by line
  istringstream sline(line);
  int eventId, jetId, __i, __type, __col, __acol;
  double __px, __py, __pz, __E, __time;
  sline >> eventId >> jetId >> __i >> __type >> __col >> __acol
  >> __E >> __px >> __py >> __pz >> __time;
  if(eventId!=evIdCurr) nEvents++;
  evIdCurr = eventId;
  jetIdIn[nEvents].push_back(jetId);
  typeIn[nEvents].push_back(__type);
  colIn[nEvents].push_back(__col);
  acolIn[nEvents].push_back(__acol);
  EIn[nEvents].push_back(__E);
  pxIn[nEvents].push_back(__px);
  pyIn[nEvents].push_back(__py);
  pzIn[nEvents].push_back(__pz);
  timeIn[nEvents].push_back(__time);
 } // end file read
 //----------------------------------------------------------
 for(int iEvent=0; iEvent<nEvents+1; iEvent++) { // find jets in each event
  vector<fastjet::PseudoJet> input_particles;
  input_particles.clear();

    // input_particles.push_back(fastjet::PseudoJet(px,py,pz,E));
  for(int i=0; i<EIn[iEvent].size(); i++) {
   //double pt = sqrt(pow(pxIn[iEvent].at(i),2) + pow(pyIn[iEvent].at(i),2));
   //if(pt > 0.5)  // no difference for vacuum jets!
   //double pmod = sqrt(pow(pxIn[iEvent].at(i),2) + pow(pyIn[iEvent].at(i),2) + 
   //  pow(pzIn[iEvent].at(i),2));
   //double rap = 0.5*log((pmod+pzIn[iEvent].at(i))/(pmod-pzIn[iEvent].at(i)));
   //if(fabs(rap)<1.5)  // no unexpected difference for vacuum jets!
   fastjet::PseudoJet particle (pxIn[iEvent].at(i),
     pyIn[iEvent].at(i), pzIn[iEvent].at(i), EIn[iEvent].at(i));
   particle.set_user_index(jetIdIn[iEvent].at(i));
   input_particles.push_back(particle);
  }
  // create a jet definition: 
  // a jet algorithm with a given radius parameter
  double R = 0.4;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);


  // get the resulting jets ordered in pt
  double ptmin = 3.0;
  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
  // tell the user what was done
  //  - the description of the algorithm used
  //  - extract the inclusive jets with pt > 5 GeV
  //    show the output as 
  //      {index, rap, phi, pt}
  //cout << "Ran " << jet_def.description() << endl;

  // label the columns
  //printf("%5s %15s %15s %15s\n","jet #", "rapidity", "phi", "pt");
 
  // print out the details for each jet
  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
   //originally: printf("%5u %15.8f %15.8f %15.8f\n",
   //i, inclusive_jets[i].rap(), inclusive_jets[i].phi(),
   //inclusive_jets[i].perp());
   vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
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
   fout << setw(8) << iEvent << setw(14) << inclusive_jets[i].E()
      << setw(14) << inclusive_jets[i].px() << setw(14) << inclusive_jets[i].py()
      << setw(14) << inclusive_jets[i].pz()
      << setw(14) << listFractions[0] << setw(14) << listFractions[1]
      << setw(14) << listFractions[2] << setw(14) << listFractions[3] << endl;
   // print selected jets on screen
   //if(leadingTrigg && inclusive_jets[i].perp()>20.0 && inclusive_jets[i].perp()<24.0) {
    //cout << setw(8) << constituents.size() << setw(14) << inclusive_jets[i].E()
      //<< setw(14) << inclusive_jets[i].perp() << setw(14) << inclusive_jets[i].phi()
      //<< setw(14) << inclusive_jets[i].rap() << endl;
    //vector<fastjet::PseudoJet> constituents = inclusive_jets[i].constituents();
    //for(fastjet::PseudoJet ptl : constituents)
     //cout << "  * " << setw(14) << ptl.E() << setw(14) << ptl.perp()
      //<< setw(14) << ptl.phi() << setw(14) << ptl.rap() << endl;
   //} // end print
  }
 } // end event loop
  return 0;
}
