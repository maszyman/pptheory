#include <iostream>
#include <cfloat>
#include <cmath>
#include <complex>
using std::complex;
#include "TMath.h"
#include "TComplex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"

using namespace std;



// //pp
double ac = 57.5;
double fc = .0;
double f0s = 7.77 ;
double d0s = 2.77 ;
double d0t = 1.7;
double f0t = 5.4;

// //ppbar
// double ac = -57.5;
// double fc = .0;

//  double f0re = -0.8/0.197327;
//  double f0im = 0.8/0.197327;

// // //double f0s = 7.77 ;
// double d0s = 0. ;
// double d0t = 0.;
// // //double f0t = 5.4;

// double rstar = 3.;

// double ac = 57.5 / 0.197327;
// double fc = .0;
// double f0 = 7.77  / 0.197327;
// double d0 = 2.77  / 0.197327;
// double rstar = 4. / 0.197327;
double pi = TMath::Pi();

double Ac(double kstar) {
  return 2*pi/(kstar*ac)*(1./(TMath::Exp(2*pi/(kstar*ac))-1));
}

double calchx(double ksa) {
  double h = 0;
  double C = 0.5772;
  for (int n = 1; n < 10000; ++n) {
    h += 1./(n*(n*n+1./(ksa*ksa)));
  }
  h = h/(ksa*ksa)-C+TMath::Log(TMath::Abs(ksa));
  return h;
}

Double_t F1( Double_t z)
{
	Double_t retval = 0.;

	// This has been Taylor expanded via Mathematica 8 about several regions to
	//  get the (absolute) relative error (|exact-expansion|/|exact|) to <10^-12
	//  in each region

	if(-0.5 < z && z <= 0.5) // expand about z=0
	{
		retval = 1.
			- 0.6666666666666666 * pow(z,2)
			+ 0.26666666666666666 * pow(z,4)
			- 0.0761904761904762 * pow(z,6)
			+ 0.016931216931216932 * pow(z,8)
			- 0.0030784030784030783 * pow(z,10)
			+ 0.0004736004736004736 * pow(z,12)
			- 0.00006314672981339648 * pow(z,14)
			+ 7.4290270368701745e-6 * pow(z,16)
			- 7.820028459863341e-7 * pow(z,18)
			+ 7.447646152250801e-8 * pow(z,20);
	}
	else if(0.5 < z && z <= 1.5)  // expand about z=1
	{
		Double_t Z = z - 1.;

		retval = 0.5380795069127684
			- 0.6142385207383052 * Z
			+ 0.15231802765107363 * pow(Z,2)
			+ 0.20640164362410562 * pow(Z,3)
			- 0.15480123271807905 * pow(Z,4)
			- 0.009326800154403014 * pow(Z,5)
			+ 0.04683600747655492 * pow(Z,6)
			- 0.010659200176460604 * pow(Z,7)
			- 0.0077623034791009515 * pow(Z,8)
			+ 0.003816680446982551 * pow(Z,9)
			+ 0.0006567448905534512 * pow(Z,10)
			- 0.0007527089460839242 * pow(Z,11)
			+ 0.000023132065749664754 * pow(Z,12)
			+ 0.0001038742320756958 * pow(Z,13)
			- 0.000017792720288706265 * pow(Z,14)
			- 0.000010618987659610779 * pow(Z,15)
			+ 3.4102621797593713e-6 * pow(Z,16)
			+ 7.803770470382609e-7 * pow(Z,17)
			- 4.4503413004770564e-7 * pow(Z,18)
			- 3.138504300892464e-8 * pow(Z,19)
			+ 4.5492668632185485e-8 * pow(Z,20);
	}
	else if(1.5 < z && z <= 2.5)  // expand about z=2
	{
		Double_t Z = z - 2.;

		retval = 0.15067019446189603
			- 0.1780158750785322 * Z
			+ 0.1436992987725385 * pow(Z,2)
			- 0.07631761246557614 * pow(Z,3)
			+ 0.015281088695458992 * pow(Z,4)
			+ 0.012448814913856642 * pow(Z,5)
			- 0.011991407785542712 * pow(Z,6)
			+ 0.003551315720097834 * pow(Z,7)
			+ 0.0008882863084414465 * pow(Z,8)
			- 0.0010849209685112123 * pow(Z,9)
			+ 0.00026598277627372933 * pow(Z,10)
			+ 0.0000840510173639067 * pow(Z,11)
			- 0.00006829335918523089 * pow(Z,12)
			+ 8.837955824274462e-6 * pow(Z,13)
			+ 6.564636093571255e-6 * pow(Z,14)
			- 2.836785279195758e-6 * pow(Z,15)
			- 6.59423030131484e-8 * pow(Z,16)
			+ 3.298108124044789e-7 * pow(Z,17)
			- 6.594367055610375e-8 * pow(Z,18)
			- 1.9108494426177143e-8 * pow(Z,19)
			+ 1.0074139589320833e-8 * pow(Z,20);
	}
	else if(2.5 < z && z <= 3.5)  // expand about z=3
	{
		Double_t Z = z - 3.;

		retval = 0.05942367687018614
			- 0.04301662017784554 * Z
			+ 0.02454138018577943 * pow(Z,2)
			- 0.013112988389733158 * pow(Z,3)
			+ 0.006668535225039065 * pow(Z,4)
			- 0.003006880658043776 * pow(Z,5)
			+ 0.0010204827706693652 * pow(Z,6)
			- 0.00013174151525177814 * pow(Z,7)
			- 0.0001169480223082913 * pow(Z,8)
			+ 0.00009990849051875193 * pow(Z,9)
			- 0.000037685947877400034 * pow(Z,10)
			+ 3.875352934429673e-6 * pow(Z,11)
			+ 3.7820501516223997e-6 * pow(Z,12)
			- 2.26604856817076e-6 * pow(Z,13)
			+ 4.61386821264982e-7 * pow(Z,14)
			+ 9.783816473460321e-8 * pow(Z,15)
			- 9.023004832948148e-8 * pow(Z,16)
			+ 2.080823670953649e-8 * pow(Z,17)
			+ 2.555328313725995e-9 * pow(Z,18)
			- 2.873436692340084e-9 * pow(Z,19)
			+ 6.151696854189207e-10 * pow(Z,20);
	}
	else if(3.5 < z && z <= 4.5) // expand about z=4
	{
		Double_t Z = z - 4.;

		retval = 0.032337000309001274
			- 0.016780252549260462 * Z
			+ 0.006642072716354607 * pow(Z,2)
			- 0.0023885420751870896 * pow(Z,3)
			+ 0.000829678521473416 * pow(Z,4)
			- 0.0002882788762123043 * pow(Z,5)
			+ 0.0001023677126197116 * pow(Z,6)
			- 0.00003711556484687208 * pow(Z,7)
			+ 0.000013228029512003321 * pow(Z,8)
			- 4.256544380350717e-6 * pow(Z,9)
			+ 1.0339380366350342e-6 * pow(Z,10)
			- 6.387919351110839e-8 * pow(Z,11)
			- 1.0873406894575189e-7 * pow(Z,12)
			+ 7.398507934346277e-8 * pow(Z,13)
			- 2.7425181421793622e-8 * pow(Z,14)
			+ 5.378173253263366e-9 * pow(Z,15)
			+ 5.106316854411601e-10 * pow(Z,16)
			- 8.267847325147262e-10 * pow(Z,17)
			+ 3.1128865486429603e-10 * pow(Z,18)
			- 4.826880327270791e-11 * pow(Z,19)
			- 1.0213390532130392e-11 * pow(Z,20);
	}
	else // use the Milne numerical integration technique (because I already had it coded)
	{
		Double_t z2 = z*z;

		// Perform 128+1 point (must be 4n+1) Milne integration technique
		Double_t sa = 0., sb = 0., sc = 0., sd = 0.;

		Double_t width = z / 128.;
		Double_t x = width;

		for(int i=0; i<32; i++)
		{
			sa += exp( -z2 + x*x );
			sb += exp( -z2 + (x+   width)*(x+   width) );
			sc += exp( -z2 + (x+2.*width)*(x+2.*width) );
			sd += exp( -z2 + (x+3.*width)*(x+3.*width) );

			x += 4. * width;
		}

		// TODO check to make sure the final terms are right
		retval = 64.*sa + 24.*sb + 64.*sc + 28.*sd + 14.*(exp(-z2) - 1.);
		retval *= width / (45. * z);
	}

	return retval;
}

Double_t F2( Double_t z)
{
	return (1. - exp(-z*z)) / z;
}



// MAIN
// ___________________________________________________________________
void calcppcf () {


  int i = 0;

  double Bo =.0;
  double Bi =.0;

  TFile* fout = new TFile("CFppkpp.root","recreate");

  for (double Rinv = 2.0; Rinv <= 2.0; Rinv += 0.1) {

    TGraphErrors* ppcf = new TGraphErrors();

    TGraphErrors* ppcfStrong = new TGraphErrors();
    TGraphErrors* ppcfCoulomb = new TGraphErrors();
    TGraphErrors* ppcfQS = new TGraphErrors();

    for (double x = 0.01; x < 0.5; x += 0.002){

      // cout << x << endl;

      double kstar = x / 0.197327;
      double rstar = Rinv * sqrt(2);

      // Ac = 2*pi/(kstar*ac)*(1./(TMath::Exp(2*pi/(kstar*ac))-1));
      // cfval = Ac;
      // cout << "h = " << calchx(kstar*ac) << endl;

      Bo  = -0.5 * exp( -4. * kstar*kstar * rstar*rstar );

      TComplex fss;
      TComplex fst;
      TComplex impartfs;

      impartfs.fRe = 0.;
      impartfs.fIm = -kstar*Ac(kstar);

// // //ppbar

//       double AA = f0re/(f0re*f0re+f0im*f0im) + .5*d0s*kstar*kstar - 2./ac  * calchx(kstar*ac) ;
//       double BB = kstar*Ac(kstar) + (f0im/(f0re*f0re+f0im*f0im));

//       fss.fRe = AA / (AA*AA+BB*BB);
//       fss.fIm = BB / (AA*AA+BB*BB);

//       Bi = .5 * fss.Rho2()/(rstar*rstar) * (1-d0s/(2.*sqrt(pi)*rstar)) + 2*fss.Re()*F1(2.*kstar*rstar)/(sqrt(pi)*rstar) - fss.Im()*F2(2.*kstar*rstar)/(rstar);

//pp
      // fss = 1. / ( 1./f0s + .5*d0s*kstar*kstar - 2./ac * calchx(kstar*ac) + impartfs );
      // fst = 1. / ( 1./f0t + .5*d0t*kstar*kstar - 2./ac * calchx(kstar*ac) + impartfs );

      double AA = 1./f0s + .5*d0s*kstar*kstar - 2./ac * calchx(kstar*ac);
      double BB = -kstar*Ac(kstar);

      fss.fRe =  AA / (AA*AA+BB*BB) ;
      fss.fIm =  -BB / (AA*AA+BB*BB) ;

      AA = 1./f0t + .5*d0t*kstar*kstar - 2./ac * calchx(kstar*ac);
      BB = -kstar*Ac(kstar);

      fst.fRe =  AA / (AA*AA+BB*BB) ;
      fst.fIm =  -BB / (AA*AA+BB*BB) ;

      double Bis = .5 * fss.Rho2()/(rstar*rstar) * (1-d0s/(2.*sqrt(pi)*rstar)) + 2*fss.Re()*F1(2.*kstar*rstar)/(sqrt(pi)*rstar) - fss.Im()*F2(2.*kstar*rstar)/(rstar);

      double Bit = .5 * fst.Rho2()/(rstar*rstar) * (1-d0t/(2.*sqrt(pi)*rstar)) + 2*fst.Re()*F1(2.*kstar*rstar)/(sqrt(pi)*rstar) - fst.Im()*F2(2.*kstar*rstar)/(rstar);

      Bi = 0.25 * Bis + 0.75 * Bit;

      double Coulomb  = Ac(kstar);

       // cout << "Bi = " << Bi << endl;
       // cout << "re = " << fss.fRe << endl;
       // cout << "im = " << fss.fIm << endl;
       // cout << "AA = " << AA << endl;
       // cout << "BB = " << BB << endl;
      // cout << "Bo = " << Bo << endl;
      // cout << "Coulomb = " << Coulomb << endl << endl;

      // Coulomb = 1;
      //Bo = 0;
       //  Bi = 0;

      ppcfStrong->SetPoint(i,x,1+Bi);
      ppcfCoulomb->SetPoint(i,x,Coulomb);
      ppcfQS->SetPoint(i,x,1+Bo);

      double scale = 0.96 / 0.875;
      //double scale = 1;
      ppcf->SetPoint(i++,x,scale*Coulomb*(1+Bo+Bi));

    }

    // ppcf->Draw("ac*");
    // ppcfStrong->Draw("csame");
    // ppcfCoulomb->Draw("csame");
    // ppcfQS->Draw("csame");
    // ppcf->RemovePoint(0);
    ppcf->SetName(Form("cppR%.1f",Rinv));
    // ppcf->SetName(Form("cppR%d",(int)(Rinv*10)));
    ppcf->Write();

    // delete ppcf;
    // delete ppcfStrong;
    // delete ppcfCoulomb;
    // delete ppcfQS;

    }
 fout->Close();
}
