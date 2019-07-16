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

void outputJetTotalsWfracs(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets,
       std::map<int, std::map<int, double> >& genJetPt);
void outputFullJets(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets);
void outputFullJets_WTA(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets, double R);

int main(int argc, char* argv[]) {

 if(argc != 5) {
  cout << "usage: ./findJets  R  in_file  jet_totals_file  out_file\n";
  return 1;
 }
 ifstream fin(argv[2]);
 if(!fin) {
  cout << "cannot open input file: " << argv[2] << endl;
  return 1;
 }
 ifstream fin_full(argv[3]);
 if(!fin_full) {
  cout << "cannot open input file: " << argv[3] << endl;
  return 1;
 }
 ofstream fout(argv[4]);
 if(!fout) {
  cout << "cannot open output file: " << argv[4] << endl;
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
 evIdCurr = 2147483646;
 int nEvents2 = -1; // similar event count for the 2nd file read loop
 std::map<int, std::map<int, double> > genJetPt; // pt of generator level jets
 while (getline(fin_full, line)) {  // reading the full jets
  istringstream sline(line);
  int eventId, jetId;
  double pxtot, pytot, pztot, Etot;
  sline >> eventId >> jetId >> Etot >> pxtot >> pytot >> pztot;
  if(eventId!=evIdCurr) nEvents2++;
  evIdCurr = eventId;
  genJetPt[nEvents2][jetId] = sqrt(pxtot*pxtot+pytot*pytot);
 } // end 2nd file (full jets) read
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
  if(R = std::atof(argv[1])) {
   ; //cout << "R = " << R << endl;
  } else {
   cout << "R value undefined: " << argv[1] << endl;
   return 1;
  }
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
  // run the jet clustering with the above jet definition
  fastjet::ClusterSequence clust_seq(input_particles, jet_def);

  // get the resulting jets ordered in pt
  double ptmin = 3.0;
  vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));
 
  outputJetTotalsWfracs(fout, iEvent, inclusive_jets, genJetPt);
  //outputFullJets(fout, iEvent, inclusive_jets);
  //outputFullJets_WTA(fout, iEvent, inclusive_jets, R);

 } // end event loop
  return 0;
}


void outputJetTotalsWfracs(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets,
       std::map<int, std::map<int, double> >& genJetPt)
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
  int leadOrigin=99999;
  double maxContrib = 0.0;
  for(auto it = origins.begin(); it != origins.end(); ++it) {
   if(it->second > maxContrib) {
    maxContrib = it->second;
    leadOrigin = it->first;
   }
  }
  double genLevelPt = 0.;
  if(genJetPt.find(iEvent)!=genJetPt.end() &&
    genJetPt[iEvent].find(leadOrigin)!=genJetPt[iEvent].end())
   genLevelPt = genJetPt[iEvent][leadOrigin];
  else
   cout << "origin_not_found[" << iEvent << " " << leadOrigin << "]";
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
     << setw(14) << listFractions[2] << setw(14) << listFractions[3]
     << setw(14) << genLevelPt << endl;
 }
}


void outputFullJets(ofstream& fout, int iEvent, vector<fastjet::PseudoJet>& jets)
// output for jet structure plots
{
 for (unsigned int i = 0; i < jets.size(); i++) {
  vector<fastjet::PseudoJet> constituents = jets[i].constituents();
  bool leadingTrigg = false;
  map<int,double> origins; // map containing the indexes of jets the partons are coming from
  for(fastjet::PseudoJet ptl : constituents) { // loop1: find the leading contribution
   if(ptl.perp() > 5.0) leadingTrigg = true;
   int origJet = ptl.user_index(); // get back the index of original jet
   if(origins.count(origJet)==0)
    origins[origJet] = sqrt(ptl.px()*ptl.px()+ptl.py()*ptl.py());
   else
    origins[origJet] += sqrt(ptl.px()*ptl.px()+ptl.py()*ptl.py());//ptl.perp();
  }
  int leadOrigin=99999;
  double maxContrib = 0.0;
  for(auto it = origins.begin(); it != origins.end(); ++it) {
   if(it->second > maxContrib) {
    maxContrib = it->second;
    leadOrigin = it->first;
   }
  }
  //if(leadingTrigg)
  for(fastjet::PseudoJet ptl : constituents) {
   double phi = ptl.delta_R(jets[i]);
   fout << setw(8) << iEvent << setw(8) << i
     << setw(14) << jets[i].perp() << setw(14) << jets[i].rap()
     << setw(14) << ptl.perp() << setw(14) << phi
     << setw(3) << (ptl.user_index()==leadOrigin ? 1 : 0) << endl;
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
