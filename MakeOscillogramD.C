#include <iostream>

#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TSystem.h"

#ifndef __CINT__
#include "../PremModel.h"
#include "../PMNS_Fast.h"
#include "../PMNS_Decay.h"
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

TH2D* GetOscHist(int flvf, int flvi, int ch){

  // Use 200 x bins and 100 y bins
  int nbinsx = 100;
  int nbinsy = 100;

  // Set parameters to GlobalFit 2018
  double dm21 = 7.55e-5;
  double dm31 = 2.5e-3;
  double th12 = asin(sqrt(0.32));
  double th13 = asin(sqrt(0.0216));
  double th23 = asin(sqrt(0.547));
  double dcp  = 1.32*TMath::Pi();
    
  double a3   = 1e-5;

  // Create PMNS object
  OscProb::PMNS_Fast myPMNS;
  OscProb::PMNS_Decay yourPMNS;
  // Set PMNS parameters
  myPMNS.SetDm(2, dm21);
  myPMNS.SetDm(3, dm31);
  myPMNS.SetAngle(1,2, th12);
  myPMNS.SetAngle(1,3, th13);
  myPMNS.SetAngle(2,3, th23);
  myPMNS.SetDelta(1,3, dcp);
    
  yourPMNS.SetAlpha3(a3);



  //  Make equally spaced log bins in E
  double e_min = 1.0, e_max = 20.0, le_min = TMath::Log10(e_min), le_max = TMath::Log10(e_max); // 20.0 was 100.0
  double vx[nbinsx+1];
  for(int i = 0; i <= nbinsx; i++) vx[i] = TMath::Power(10., le_min+i*(le_max-le_min)/nbinsx);

	TH2D* h2 = new TH2D("","",nbinsx,&vx[0],nbinsy,-1,0);
    TH2D* hD = new TH2D("","",nbinsx,&vx[0],nbinsy,-1,0);

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
    for(int le=1; le<=nbinsx; le++){

      // Set E from bin center
      double E  = h2->GetXaxis()->GetBinCenter(le);

      // Initialize probability
      double prob = 0;


	// Get Probabilities
	if(ch != 1 && ch != -1){cout << " Wrong chirality! " <<endl;}
	if(ch == 1) myPMNS.SetIsNuBar(false);
        else myPMNS.SetIsNuBar(true);
        prob = (-1) * myPMNS.Prob(flvi, flvf, E);

    
      // Fill probabilities in histogram
      h2->SetBinContent(le,ct,prob);

    }// loe loop
  }// cosT loop
    
  // Loop over cos(theta_z) and L/E
  for(int ct=1; ct<=nbinsy; ct++){//loops over bins.

    // Get cos(theta_z) from bin center
    double cosT = hD->GetYaxis()->GetBinCenter(ct);

    // Set total path length L
    double L = prem.GetTotalL(cosT);

    // Skip if cosT is unphysical
    if(cosT < -1 || cosT > 1) continue;

    // Fill paths from PREM model
    prem.FillPath(cosT);

    // Set paths in OscProb
    yourPMNS.SetPath(prem.GetNuPath());

    // Loop over E for our case, up to 200 GeV
    for(int le=1; le<=nbinsx; le++){

      // Set E from bin center
      double E  = hD->GetXaxis()->GetBinCenter(le);

      // Initialize probability
      double prob = 0;


      // Get Probabilities
      if(ch != 1 && ch != -1){cout << " Wrong chirality! " <<endl;}
      if(ch == 1) yourPMNS.SetIsNuBar(false);
          else yourPMNS.SetIsNuBar(true);
          prob = yourPMNS.Prob(flvi, flvf, E);

      
        // Fill probabilities in histogram
        hD->SetBinContent(le,ct,prob);

      }// loe loop
    }// cosT loop

  h2->Add(hD);

  // Set nice histogram
  SetHist(h2);
  //h2->SetMinimum(-0.06);
  //h2->SetMaximum(0.02);
  // Set titles
  string flv1; //flvi
  string flv2; //flvf
  switch(flvi){
    case 0: flv1 = "e"; break;
    case 1: flv1 = "#mu"; break;
    case 2: flv1 = "#tau"; break;
  }
  switch(flvf){
    case 0: flv2 = "e"; break;
    case 1: flv2 = "#mu"; break;
    case 2: flv2 = "#tau"; break;
  } //TString::Format("%0.1f GeV", energy))
  //h2->SetTitle(";E (GeV);cos#theta_{z}; P[#bar{#nu}_{e} #rightarrow #bar{#nu}_{e}]");
  if (ch==1){h2->SetTitle(TString::Format(";E (GeV);cos#theta_{z}; P[#nu_{%s} #rightarrow #nu_{%s}]", flv1.c_str(), flv2.c_str()));}
  else {h2->SetTitle(TString::Format(";E (GeV);cos#theta_{z}; P[#bar{#nu}_{%s} #rightarrow #bar{#nu}_{%s}]", flv1.c_str(), flv2.c_str()));}

  return h2;

}


// Make an oscillogram for NH
// nue (0), numu (1) or nutau (2)
void MakeOscillogramD(int flvf = 0, int flvi = 1, int ch = 1){// main function

  // Load OscProb library
  if(isCINT) LoadOscProb();

  // Set a nice overall style
  SetNiceStyle();

    const Int_t NRGBs = 8; const Int_t NCont = 256;
    Double_t stops[NRGBs] = { 0.00, 0.33, 0.50, 0.66, 0.83, 0.84, 0.92, 1.00 };
    Double_t red[NRGBs] = { 0.00, 0.00, 0.21, 0.63, 1.00, 1.00, 0.87, 0.30 };
    Double_t green[NRGBs] = { 0.00, 0.05, 0.26, 0.63, 1.00, 1.00, 0.18, 0.00 };
    Double_t blue[NRGBs] = { 0.30, 0.60, 0.96, 0.98, 1.00, 1.00, 0.13, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    
  TCanvas *c = new TCanvas();
  // Make the oscillogram
  TH2D* hNH = GetOscHist(flvf, flvi, ch);

  // Draw the oscillogram
  hNH->Draw("colz");
  //gStyle->SetPalette();

  // Add space for the colz bar
  gPad->SetRightMargin(0.18);
  //gPad->SetLogx();

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
