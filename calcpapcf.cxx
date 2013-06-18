#include "fortranc.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TGraph.h"

using namespace std;

// --- Prototype of the function used in the weight calculator
//     (in FsiWeightLedinicky.F)
#define fsiini F77_NAME(fsiini,FSIINI)
extern "C" {void type_of_call F77_NAME(fsiini,FSIINI)(const int &itest,const int &ll,const int &ns,const int &ich, const int &iqs, const int &isi,const int &i3c);}
#define fsinucl F77_NAME(fsinucl,FSINUCL)
extern "C" {void type_of_call  F77_NAME(fsinucl,FSINUCL)(const double &mn,const double &cn);}
#define fsimomentum F77_NAME(fsimomentum,FSIMOMENTUM)
extern "C" {void type_of_call F77_NAME(fsimomentum,FSIMOMENTUM)(double &p1,double &p2);}
#define fsiposition F77_NAME(fsiposition,FSIPOSITION)
extern "C" {void type_of_call F77_NAME(fsiposition,FSIPOSITION)(double &x1,double &x2);}
#define fsiw F77_NAME(fsiw,FSIW)
extern "C" {void type_of_call F77_NAME(fsiw,FSIW)(const int &i,double &weif,
                                                  double &wei,double &wein);}
#define ltran12 F77_NAME(ltran12,LTRAN12)
extern "C" {void type_of_call ltran12_();}


// --- global variables...
int mLL=2;
int mNs=4;//4
int mItest=1;
int mIch=0;
int mIqs=0;
int mIsi=1;
int mI3c=0;

double mNuclMass=1.;
double mNuclCharge=1.;
double mNuclChargeSign=1.0;

const int nPart = 1e6;

TFile* fin;
TH1D* hptprotons;
TH1D* hptantiprotons;
TFile* finxyzt;
TH1D* hx;
TH1D* hy;
TH1D* hz;
TH1D* ht;

// --- calculating weights...
double getWeight(double* momentum1, double* momentum2, double* position1, double* position2) {

  // cout << momentum1[0] << " asd2 " << endl;
  // cout << momentum2[0] << " asd2 " << endl;

  // cout << position1[0] << " asd2 " << endl;
  // cout << position2[0] << " asd2 " << endl;

  // // Fsi weight output
  double mWei;
  double mWein;
  double mWeif;
  double mWeightDen;
  fsimomentum(*momentum1,*momentum2);
  fsiposition(*position1,*position2);
  ltran12();
  fsiw(1,mWeif,mWei,mWein);
  return mWei;
}

// --- initialize fsi...
void InitFsi() {
  fsiini(mItest,mLL,mNs,mIch,mIqs,mIsi,mI3c);
  fsinucl(mNuclMass,mNuclCharge);
}

double* generatePosition() {

  double* position1 = new double[4];

  position1[0] = hx->GetRandom();
  position1[1] = hy->GetRandom();
  position1[2] = hz->GetRandom();
  position1[3] = ht->GetRandom();

  // cout << position1[0] << " qwe " << endl;

  return position1;
}

double* generateMomentum() {

  double* p1 = new double[4];

//  double p1[4];
//  double p2[4];

  double pt1 = hptprotons->GetRandom();
//  double pt2 = hptprotons->GetRandom();

  double phi1 = gRandom->Uniform(0,2*TMath::Pi());
//  double phi2 = gRandom->Uniform(0,2*TMath::Pi());

  double eta1 = gRandom->Uniform(-0.8,0.8);
//  double eta2 = gRandom->Uniform(-0.8,0.8);

  if ( gRandom->Uniform(0,1) < 0.5 ) {
    p1[0] = TMath::Sqrt( TMath::Power(pt1,2) / ( 1+(TMath::Power(TMath::Tan(phi1),2)) ) );
//    p2[0] = TMath::Sqrt( TMath::Power(pt2,2) / ( 1+(TMath::Power(TMath::Tan(phi2),2)) ) );
  }
  else {
    p1[0] = -TMath::Sqrt( TMath::Power(pt1,2) / ( 1+(TMath::Power(TMath::Tan(phi1),2)) ) );
//    p2[0] = -TMath::Sqrt( TMath::Power(pt2,2) / ( 1+(TMath::Power(TMath::Tan(phi2),2)) ) );
  }

  p1[1] = p1[0]*TMath::Tan(phi1);
//  p2[1] = p2[0]*TMath::Tan(phi2);

  double theta1 = 2*TMath::ATan(TMath::Exp(-eta1));
//  double theta2 = 2*TMath::ATan(TMath::Exp(-eta2));

  p1[2] = pt1/TMath::Tan(theta1);
//  p2[2] = pt2/TMath::Tan(theta2);

  double mom1 = TMath::Sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
//  double mom2 = TMath::Sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);

  p1[3] = TMath::Sqrt(0.938*0.938+mom1*mom1);
//  p2[3] = TMath::Sqrt(0.938*0.938+mom2*mom2);


//  cout << p1[0] << " qwe " << endl;

  return p1;

}

void getDist() {
  fin = new TFile("ptdist.root","read");
  hptprotons = (TH1D*)fin->Get("PtPhicutPass1PPtpcM0Psi6_px");
  hptantiprotons = (TH1D*)fin->Get("PtPhicutPass1APAPtpcM0Psi6_px");
  finxyzt = new TFile("xyztdist.root","read");
  hx = (TH1D*)finxyzt->Get("hx");
  hy = (TH1D*)finxyzt->Get("hy");
  hz = (TH1D*)finxyzt->Get("hz");
  ht = (TH1D*)finxyzt->Get("ht");

}


double calckstar(double* p1, double* p2) {

  double tPx = p1[0]+p2[0];
  double tPy = p1[1]+p2[1];
  double tPz = p1[2]+p2[2];
  double tE  = p1[3]+p2[3];
  double tPt = tPx*tPx + tPy*tPy;

  double tMt = tE*tE - tPz*tPz;
  double tM = TMath::Sqrt(tMt - tPt);
  tMt = TMath::Sqrt(tMt);
  tPt = TMath::Sqrt(tPt);

  // Boost to LCMS
  double tBeta = tPz/tE;
  double tGamma = tE/tMt;
  double mKStarLong = tGamma * (p1[2] - tBeta * p1[3]);
  double tE1L = tGamma * (p1[3]  - tBeta * p1[2]);

  // Rotate in transverse plane
  double mKStarOut  = ( p1[0]*tPx + p1[1]*tPy)/tPt;
  double mKStarSide = (-p1[0]*tPy + p1[1]*tPx)/tPt;

  // Boost to pair cms
  mKStarOut = tMt/tM * (mKStarOut - tPt/tMt * tE1L);
  double mKStarSigned = sqrt(mKStarSide*mKStarSide + mKStarOut*mKStarOut + mKStarLong*mKStarLong);

  return mKStarSigned;
}

// --- main code...
int main() {

  srand(time(NULL));
  gRandom->SetSeed(0);

  getDist();

  // double momentum1[4];
  // double momentum2[4];
  // double position1[4];
  // double position2[4];
  // double p1[4];
  // double p2[4];

  // double pt1 = hptprotons->GetRandom();
  // double pt2 = hptprotons->GetRandom();

  // double phi1 = gRandom->Uniform(0,2*TMath::Pi());
  // double phi2 = gRandom->Uniform(0,2*TMath::Pi());

  // double eta1 = gRandom->Uniform(-0.8,0.8);
  // double eta2 = gRandom->Uniform(-0.8,0.8);

  // if ( gRandom->Uniform(0,1) < 0.5 ) {
  //   p1[0] = TMath::Sqrt( TMath::Power(pt1,2) / ( 1+(TMath::Power(TMath::Tan(phi1),2)) ) );
  //   p2[0] = TMath::Sqrt( TMath::Power(pt2,2) / ( 1+(TMath::Power(TMath::Tan(phi2),2)) ) );
  // }
  // else {
  //   p1[0] = -TMath::Sqrt( TMath::Power(pt1,2) / ( 1+(TMath::Power(TMath::Tan(phi1),2)) ) );
  //   p2[0] = -TMath::Sqrt( TMath::Power(pt2,2) / ( 1+(TMath::Power(TMath::Tan(phi2),2)) ) );
  // }

  // p1[1] = p1[0]*TMath::Tan(phi1);
  // p2[1] = p2[0]*TMath::Tan(phi2);

  // double theta1 = 2*TMath::ATan(TMath::Exp(-eta1));
  // double theta2 = 2*TMath::ATan(TMath::Exp(-eta2));

  // p1[2] = pt1/TMath::Tan(theta1);
  // p2[2] = pt2/TMath::Tan(theta2);

  // double mom1 = TMath::Sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]);
  // double mom2 = TMath::Sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);

  // p1[3] = TMath::Sqrt(0.938*0.938+mom1*mom1);
  // p2[3] = TMath::Sqrt(0.938*0.938+mom2*mom2);


  // momentum1[0] = p1[0];
  // momentum1[1] = p1[1];
  // momentum1[2] = p1[2];
  // momentum1[3] = p1[3];

  // momentum2[0] = p2[0];
  // momentum2[1] = p2[1];
  // momentum2[2] = p2[2];
  // momentum2[3] = p2[3];

  // position1[0] = hx->GetRandom();
  // position1[1] = hy->GetRandom();
  // position1[2] = hz->GetRandom();
  // position1[3] = ht->GetRandom();

  // position2[0] = hx->GetRandom();
  // position2[1] = hy->GetRandom();
  // position2[2] = hz->GetRandom();
  // position2[3] = ht->GetRandom();

  // // Fsi weight output
  // double mWei;
  // double mWein;
  // double mWeif;
  // double mWeightDen;


  InitFsi();

  double* momentum1 = new double[4];
  double* momentum2 = new double[4];
  double* position1 = new double[4];
  double* position2 = new double[4];

  TGraph* grcfpap = new TGraph();
  TH1D* numhcfpap = new TH1D("numhcfpap","numhcfpap",100,0,0.5);
  TH1D* denhcfpap = new TH1D("denhcfpap","denhcfpap",100,0,0.5);

  numhcfpap->Sumw2();
  denhcfpap->Sumw2();

  for (int iPart = 0; iPart < nPart; ++iPart) {

    momentum1 = generateMomentum();
    momentum2 = generateMomentum();
    position1 = generatePosition();
    position2 = generatePosition();

    // cout << momentum1[0] << " asd " << endl;
    // cout << momentum2[0] << " asd " << endl;

    // cout << position1[0] << " asd " << endl;
    // cout << position2[0] << " asd " << endl;

    double kstar = calckstar(momentum1,momentum2);
    if (kstar > 0.5) {
      --iPart;
      continue;
    }
    double mWei = getWeight(momentum1,momentum2,position1,position2);

    //grcfpap->SetPoint(iPart,kstar,mWei);
    numhcfpap->Fill(kstar,mWei);
    denhcfpap->Fill(kstar,1);

    // fsimomentum(*momentum1,*momentum2);
    // fsiposition(*position1,*position2);
    // ltran12();
    // fsiw(1,mWeif,mWei,mWein);

    // cout << "kstar = " << kstar << endl;
    // cout << "mWei = " << mWei << endl;
    // cout << "asd = " << hx->GetRandom() << endl;
    // cout << "asd = " << hx->GetRandom() << endl;
  }

  // grcfpap->Draw("acp");

  TH1D* cfpap = new TH1D("cfpap","cfpap",100,0,0.5);
  cfpap->Sumw2();
  cfpap->Divide(numhcfpap,denhcfpap,1.0,1.0);

  TFile* fout = new TFile("CFpapkpap.root","recreate");
  //grcfpap->Write();
  numhcfpap->Write();
  denhcfpap->Write();
  cfpap->Write();

  fout->Close();

  delete momentum1;
  delete momentum2;
  delete position1;
  delete position2;

  return 0;
}
