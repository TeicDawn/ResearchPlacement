#include <iostream>

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"

#ifndef __CINT__
#include "../PremModel.h"
#include "../PMNS_Fast.h"
//bool isCINT = false;

bool isCINT = true;
#endif



// Macro to load OscProb library
#include "LoadOscProb.C"

// Some functions to make nice plots
#include "SetNiceStyle.C"


//========== Make oscillogram for P(flvi->flvf, MH, Chirality) =============
// Input:
//      flvf, flvi: final and initial flavour of the transition
//                 0: nue
//                 1: numu
//                 2: nutau
//
//      mh: mass hierarchy, 1=normal, -1 = inverted.
//      ch: chirality, 1=neutrino, -1= antineutrino
//  Output:
//      oscillogram, cos(theta) vs E
//  Inside Oscprob/tutorial, run as root -l LoadOscProb.C 'MakeOscillogram.C(flvf,flvi,mh,ch)'.
//  Example: root -l LoadOscProb.C 'MakeOscillogram.C(2,1,1,1)' gives oscillogram for anit-numu->anit-nutau, NH.

TH2D* GetOscHist(int flvf, int mh){

  // Use 200 x bins and 100 y bins
  int nbinsx = 200;
  int nbinsy = 100;

  // Set parameters to nuFit 5.0
  double dm21 = 7.42e-5;
  double dm31 = mh>0 ? 2.457e-3 : -2.449e-3 + dm21;
  double th12 = asin(sqrt(0.304));
  double th13 = asin(sqrt(0.02219));
  double th23 = asin(sqrt(mh>0 ? 0.553 : 0.589));
  double dcp  = (mh>0 ? 306 : 254)*TMath::Pi()/180;

  // Create PMNS object
  OscProb::PMNS_Fast myPMNS;

  // Set PMNS parameters
  myPMNS.SetDm(2, dm21);
  myPMNS.SetDm(3, dm31);
  myPMNS.SetAngle(1,2, th12);
  myPMNS.SetAngle(1,3, th13);
  myPMNS.SetAngle(2,3, th23);
  myPMNS.SetDelta(1,3, dcp);



  //  Make equally spaced log bins in E
  // double e_min = 1.0, e_max = 100.0, le_min = TMath::Log10(e_min), le_max = TMath::Log10(e_max);
  // double vx[nbinsx+1];
  // for(int i = 0; i <= nbinsx; i++) vx[i] = TMath::Power(10., le_min+i*(le_max-le_min)/nbinsx);

  TH2D* h2 = new TH2D("","",nbinsx,0,50*nbinsx,nbinsy,-1,0);

  // Create default PREM Model
  OscProb::PremModel prem;

  // Loop over cos(theta_z) and L/E
  for(int ct=1; ct<=nbinsy; ct++){//loops over bins.

    // Get cos(theta_z) from bin center
    double cosT = h2->GetYaxis()->GetBinCenter(ct);

    // Set total path length L
    double L = prem.GetTotalL(cosT);

    // Skip if cosT is unphysical
    if(cosT < -1 || cosT > 1) continue;

    // Fill paths from PREM model
    prem.FillPath(cosT);

    // Set paths in OscProb
    myPMNS.SetPath(prem.GetNuPath());

    // Loop over E for our case, up to 200 GeV
    for(int le = 1; le <= nbinsx; le++){

      // Set E from bin center
      double loe  = h2->GetXaxis()->GetBinCenter(le);

      // Initialize probability
      double prob = 0;


      // Get Probabilities
      for (int j = 0; j <= 1; ++j){
        if(j == 0) myPMNS.SetIsNuBar(false);
        else myPMNS.SetIsNuBar(true);
        prob += (1-j*0.5)*(myPMNS.Prob(1, 1, L/loe) + 0.5*myPMNS.Prob(0, 1, L/loe));
      }
    
    // Fill probabilities in histogram
    h2->SetBinContent(le,ct,prob);

    }// loe loop
  }// cosT loop


  // Set nice histogram
  SetHist(h2);

  // Set titles
  h2->SetTitle(";L/E (km/GeV);cos#theta_{z};P_{#mu#mu} + 0.5#timesP_{#bar{#mu#mu}} + 0.5#timesP_{e#mu} + 0.25#timesP_{#bar{e#mu}}");
    
  return h2;
}


// Make an oscillogram for NH
// nue (0), numu (1) or nutau (2)
void MakeOscillogramLoE(int flvf = 1){// main function

  // Load OscProb library
  if(isCINT) LoadOscProb();

  // Set a nice overall style
  SetNiceStyle();

  TCanvas *c = new TCanvas();
  // Make the oscillogram
  TH2D* hNH = GetOscHist(1, 1);

  // Draw the oscillogram
  hNH->Draw("colz");

  // Add space for the colz bar
  gPad->SetRightMargin(0.18);
  // gPad->SetLogx();
  

}




// Make oscillogram for given final flavour and MH

// Get pad NDC from x value
double GetNDCx(double x) {

  return (x - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());

}

// Get pad NDC from y value
double GetNDCy(double y) {

  return (y - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());

}

// Draw energy lines
void DrawEnergyLines(TH2* hNH){

  // Get max value of x-axis
  double xmax = hNH->GetXaxis()->GetXmax();

  // Define earth diameter
  double dEarth = 2*6371;

  // Define constant energy line
  TF1* f= new TF1("f","-x*[0]/[1]",0,xmax);

  // Second parameter is earth's diameter
  f->SetParameter(1,dEarth);

  // Dashed black lines
  f->SetLineStyle(7);
  f->SetLineColor(kBlack);

  // Define min and max energies
  // for drawing lines
  double minE = 0.1 * dEarth/xmax;
  double maxE = 10. * dEarth/xmax;

  // Start from a rounded power of 10 energy
  double sampleE = pow(10, floor(log10(minE)) );

  // Define how to increase energies
  int idx = 0;
  double scale[3] = {2,2.5,2};

  // Make a set of energies to draw
  vector<double> energies;
  while(sampleE<maxE){

    if(sampleE>minE){
      energies.push_back(sampleE);
    }

    sampleE *= scale[idx%3];
    idx++;
  }

  // Update pad to get correct axis ranges
  gPad->Update();

  // Loop over energies and draw lines
  for(int i=0; i<int(energies.size()); i++){

    double energy = energies[i];

    // Set first parameter to energy
    f->SetParameter(0, energy);

    // Draw line
    f->DrawCopy("same");

    // Get position in pad to draw label
    double xval = 2*xmax*6371/(2*6371+xmax*energy);
    double ndcx = GetNDCx(xval);
    double ndcy = GetNDCy(f->Eval(xval));

    // Draw labels
    if(energy>=1) MiscText(ndcx,ndcy,0.03,TString::Format("%d GeV", int(energy)));
    if(energy<1)  MiscText(ndcx,ndcy,0.03,TString::Format("%0.1f GeV", energy));

  }// energies loop

}
