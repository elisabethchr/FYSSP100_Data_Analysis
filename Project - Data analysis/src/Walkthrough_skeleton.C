/*
* Project:        Exercises 11.1-11.3
* File:           Walkthrough_skeleton.C
* Author:         Ivo van Vulpen, Aart Heijboer
* Version (date): 1.0 (23.06.2013)
*
* Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
* All rights reserved.
*
* Description:
* A code skeleton for the searches part.
*
* This code is distributed with the solution manual to the book
*
* Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
* Wiley-VCH (2013),
* O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
*/


#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TGraph.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h" // for ROOT::Math::gaussian_cdf

#include <iostream>
using namespace std;

//-------------------------------------------------------------------------------------------------------------------------
/* DECLARATION OF FUNCTIONS (HEADER) */
//-- full functions
TH1D * GetMassDistribution(int Itype = 1, double scalefactor = 1.00);		// 1D-histogram with a double per channel
void MassPlot(int Irebin = 20);

//-- skeleton functions
void SideBandFit(string type, int Irebin = 10);
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig);
TH1D * TestStatistic_DistributionMC(TH1D *h_mass_template, TH1D *h_bgr, TH1D *h_sig, int Ncycle);
TH1D * GenerateToyDataSet(TH1D *h_mass_template);

//-- some useful functions
double IntegratePoissonFromRight(double mu, int N_obs);
double IntegrateFromRight(TH1D * h_X_bgr, double X_value);
double ExpectedSignificance_ToyMC(double Nbgr, double Nsig, int Ndata, double dNbgr, int N, bool obs);
void Significance_Likelihood_ToyMC(string mass_distribution_type, double signal_sf=1.0, double lum_sf=1.0, int Ncycle=10000);
void Get_ConfidenceLevels(double signal_sf=1.00, double lum_sf=1.0);
void Plot_ScaleFactors();
void SignalScaleFactors();
void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045,
Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color =1 );
vector<double> Get_Quantiles( TH1D* hist );

//-- write to file functions
// void WriteToFile();
//-------------------------------------------------------------------------------------------------------------------------



//========================================
// S O M E   F I N A L   F U N C T I O N S
//========================================


//========================================================
TH1D * GetMassDistribution(int Itype, double scalefactor){
	//========================================================
	//----------------------------------------------------------
	// Goal: return the histogram of the 4-lepton invariant mass
	//       for given type with an optional scale factor
	//
	//       Itype 1 = ZZ SM background
	//             2 = data
	//           125 = Higgs 125
	//           200 = Higgs 200
	//
	//      scalefactor: histograms will be scaled with this number
	//
	//  Note: Histograms have ~200 MeV bins, so need to rebin
	//---------------------------------------------------------

	//-- [1] Get histogram from the file
	TH1D *h_mass = 0;
	TDirectory* dir = gDirectory;
	TFile *file = new TFile("Histograms_fake.root", "READ");
	dir->cd();

	//-- Higgs 125
	if(Itype == 125){
		h_mass  = (TH1D*) file->Get("h_m4l_Higgs125_fake")->Clone("h_mass");
	}
	//-- Higgs 200
	if(Itype == 200){
		h_mass  = (TH1D*) file->Get("h_m4l_Higgs200_fake")->Clone("h_mass");
	}
	//-- ZZ SM background
	if(Itype == 1){
		h_mass  = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");
	}
	//-- data
	if(Itype == 2){
		h_mass  = (TH1D*) file->Get("h_m4l_data_fake")->Clone("h_mass");
	}

	//-- [2] scale histograms
	int Nbins = h_mass->GetNbinsX();
	for (int i_bin = 1; i_bin < Nbins; i_bin++){
		double mu_bin = h_mass->GetBinContent(i_bin);
		h_mass -> SetBinContent( i_bin, scalefactor * mu_bin);
	}


	file->Close();
	//-- [3] return histogram
	return h_mass;

	//===========================
} // end GetMassDistribution()
//===========================




//========================
void MassPlot(int Irebin){
	//========================
	// ------------------------------------------
	// Goal: produce SM+Higgs+data plot
	//       Note: rebinning is only for plotting
	// ------------------------------------------

	//------------------------------------
	//-- Standard stuff and prepare canvas
	//------------------------------------
	gROOT->Clear();
	gROOT->Delete();

	//-- Prepare canvas and plot histograms
	TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
	canvas1->SetLeftMargin(0.125);
	canvas1->SetBottomMargin(0.125);
	canvas1->cd();

	//------------------------------------------------------------------
	//-- [1] Prepare histograms
	//--     o Get histograms from the files (signal, background and data)
	//--     o Make cumulative histograms (for signal and background)
	//------------------------------------------------------------------

	//-- Get histograms from the files (higgs, zz and data)
	TH1D *h_sig, *h_bgr, *h_data;
	h_sig  = GetMassDistribution(125);
	h_bgr  = GetMassDistribution(1);
	h_data = GetMassDistribution(2);

	//-----------------------------------
	//-- [2] Plot histograms and make gif
	//--     o rebin histograms
	//--     o prepare cumulative histogram
	//--     o make plot + opsmuk + gif
	//-----------------------------------

	//-- Rebin histograms (only for plotting)
	h_sig->Rebin(Irebin);
	h_bgr->Rebin(Irebin);
	h_data->Rebin(Irebin);

	//-- Prepare cumulative histogram for signal + background
	TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
	h_sig_plus_bgr->Reset();
	for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
		h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
		printf("  REBINNED HISTOGRAM:  bin %d, Ndata = %d\n",i_bin,(int)h_data->GetBinContent(i_bin));
	}

	//-- prepare histograms and plot them on canvas
	double Data_max = h_data->GetBinContent(h_data->GetMaximumBin());
	double Ymax_plot = 1.10* (Data_max + TMath::Sqrt(Data_max));
	h_sig_plus_bgr->SetFillColor(7);
	h_sig_plus_bgr->SetAxisRange(0.,Ymax_plot,"Y");
	h_sig_plus_bgr->SetAxisRange(0.,400.,"X");
	h_bgr->SetFillColor(2);
	h_sig_plus_bgr->Draw("hist");
	h_bgr->Draw("same");
	h_bgr->Draw("axis same");
	h_data->Draw("e same");

	//-- some nice axes and add legend
	AddText( 0.900, 0.035, "4-lepton invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
	AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
	TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.85);
	leg1->SetBorderSize(0); leg1->SetFillColor(0);
	TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.04);
	TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.04);
	leg1->Draw();


	//-- prepare gif
	canvas1->Print(Form("../../data/11_1/MassPlot_rebin%d.pdf",Irebin));

	return;
}



//===============================================
// S O M E   S K E L E T O N    F U N C T I O N S
//===============================================


//=============================================================
void Significance_Optimization(double Lumi_scalefactor){
	//=============================================================

	printf("\n Significance_Optimization()\n\n");

	//------------------------------------------------------------------
	//-- [1] Prepare histograms
	//--     o Get histograms from the files (signal, background and data)
	//--     o scale to correct luminosity
	//------------------------------------------------------------------


	//-------------------------------------------------------
	//-- [2] Compute significance for various mass windows
	//--     o try various options for window (use histogram)
	//--     o compute expected and observed significance
	//-------------------------------------------------------

	int lum_max = Lumi_scalefactor;

	vector<double> lumi_scalefactors;
	lumi_scalefactors.push_back(1.0);
	int n = 30;
	double dLum = (lum_max - 1.0)/((double) n);
	for(int i = 1; i<n+1; i++){
		lumi_scalefactors.push_back(lumi_scalefactors[i-1]+dLum);
		printf("lumi_scalefactors[%d] = %5.2f\n", i, lumi_scalefactors[i]);
	}

	TH1D *h_lumi_significances = new TH1D("h_lumi_significances","Luminosity scale factors",n,1.,lum_max);  // histogram to hold results for observed
	vector<double> lumi_significances;

	for (int i_lum = 0; i_lum < n+1; i_lum++){

		printf("Lumi_scalefactor = %5.2f\n", lumi_scalefactors[i_lum]);
		//-- Get histograms from the files (higgs, zz and data)
		TH1D *h_sig, *h_bgr, *h_data;
		printf ("\n  INFO: Mass distribution in the 4 lepton channel\n");
		h_sig  = GetMassDistribution(125, lumi_scalefactors[i_lum]);
		h_bgr  = GetMassDistribution(1, lumi_scalefactors[i_lum]);
		h_data = GetMassDistribution(2, lumi_scalefactors[i_lum]);

		//-- Define histogram that defines mass windows to try
		TH1D *h_masswindow          = new TH1D("h_masswindow","",250,0.,25.);           // make a mass window - full width between 0 and 25 GeV
		TH1D *h_masswindow_expected = new TH1D("h_masswindow_expected","",250,0.,25.);  // histogram to hold results for expected
		TH1D *h_masswindow_observed = new TH1D("h_masswindow_observed","",250,0.,25.);  // histogram to hold results for observed

		//---------------------------------
		//-- Loop over various mass windows
		//---------------------------------
		for (int i_bin = 1; i_bin<h_masswindow->GetNbinsX(); i_bin++ ){	//loop over 250 bins in range 0-25GeV, such that e.g. i_bin = 1 corresponds to starting at 0.1GeV, and that the mass window is then 0.05
			// i_bin = 2 corresponds to 0.3Gev, with the mass window 0.15+-0.15

			//-- get full width of mass window (bin center of our histogram) and the number of events in mass window for each event type
			double masswindow_fullwidth = h_masswindow->GetBinCenter(i_bin);
			// printf("   Trying as mass window: %5.2f GeV\n", masswindow_fullwidth);

			// [a] determine the number of events in the mass window for each event type
			// Ndata_win, Nbgr_win and Nsig_win (Nsig_win as in mass window for H1 = sig+bkg)
			int Ndata_win, bin_min, bin_max;
			double x_min, x_max, Nbgr_win, Nsig_win;
			TAxis *xAxis;

			xAxis = h_bgr->GetXaxis();
			x_min = 125.-0.5*masswindow_fullwidth;
			x_max = 125+0.5*masswindow_fullwidth;
			bin_min = xAxis->FindBin(x_min);
			bin_max = xAxis->FindBin(x_max);

			Ndata_win = h_data->Integral(bin_min, bin_max);
			Nbgr_win = h_bgr->Integral(bin_min, bin_max);
			Nsig_win = h_sig->Integral(bin_min, bin_max);


			// [b] compute EXPECTED significance and save in histogram
			//    double pvalue_expected       = ROOT::Math::poisson_cdf_c(Nbgr_win+Nsig_win, Nbgr_win); // _cdf calculates cumulative distribution function, i.e. int(Poisson)_(-inf)^(Nbgr+Nsig) (lower tail),
			double pvalue_expected = IntegratePoissonFromRight(Nbgr_win, Nbgr_win+Nsig_win);
			/* If number observed under b+s < 1 then skip, (otherwise p = 1.0000 and Z goes to -inf) */
			if (pvalue_expected == 1.00){ printf("pvalue_expected = %5.9f\n", pvalue_expected); pvalue_expected -= 1e-9; }//continue; }

			double significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1); // can use a gaussian approximation -> see the Central Limit Theorem
			h_masswindow_expected->SetBinContent(i_bin, significance_expected);

			// [c] compute OBSERVED significance and save in histogram
			double pvalue_observed       = IntegratePoissonFromRight(Nbgr_win, Ndata_win);//IntegratePoissonFromRight(Nbgr_win, Ndata_win); // you need to do this yourself
			double significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
			h_masswindow_observed->SetBinContent(i_bin, significance_observed);
			//	    printf("           expected significance = %5.14f || observed significance = %5.14f \n", significance_expected, significance_observed);

		}

		// calculate the optimum significance
		double masswindow_exp, masswindow_obs, Z_exp, Z_obs;
		int bin_max;

		// expected
		bin_max = h_masswindow_expected->GetMaximumBin();
		masswindow_exp = h_masswindow_expected->GetBinCenter(bin_max);
		Z_exp = h_masswindow_expected->GetBinContent(bin_max);

		// observed
		bin_max = h_masswindow_observed->GetMaximumBin();
		masswindow_obs = h_masswindow_observed->GetBinCenter(bin_max);
		Z_obs = h_masswindow_observed->GetBinContent(bin_max);

		// print optimum to the screen
		printf("\n-------------------------------------------\n");
		printf("Optimum values for lumi_scalefactor = %5.2f\n", lumi_scalefactors[i_lum]);
		printf("masswindow_exp = %5.2f, Z_exp = %5.2f\n", masswindow_exp, Z_exp);
		printf("masswindow_obs = %5.2f, Z_obs = %5.2f\n", masswindow_obs, Z_obs);
		printf("-------------------------------------------\n");

		lumi_significances.push_back(Z_exp);
		h_lumi_significances->SetBinContent(i_lum, Z_exp);

		//----------------------------------
		//-- [3] Plot histogram and make gif
		//----------------------------------
		TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
		canvas1->SetLeftMargin(0.125);
		canvas1->SetBottomMargin(0.125);
		canvas1->cd();

		h_masswindow_expected->SetLineColor(1);
		h_masswindow_expected->SetLineWidth(2);
		h_masswindow_observed->SetLineColor(4);
		h_masswindow_observed->SetLineWidth(2);

		h_masswindow_expected->SetAxisRange(0,6.,"Y");
		h_masswindow_expected->Draw("l");
		if(fabs(lumi_scalefactors[i_lum]-1.00)<0.01){
			h_masswindow_observed->Draw("l same");
		}
		//-- axes
		AddText( 0.900, 0.035, "Mass window GeV",0.060, 0.,"right"); // X-axis
		AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis
		AddText( 0.225, 0.825, Form("Luminosity scalefactor = %5.2f",lumi_scalefactors[i_lum]),0.050, 0.,"left");

		AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);
		if(fabs(lumi_scalefactors[i_lum]-1.00)<0.01){
			AddText( 0.700, 0.300, "Observed significance",0.050, 0.,"right",1,4);
		}
		//-- prepare pdf
		//	  canvas1->Print(Form("../../data/11_1/Significance_Optimization_lumiscalefactor%5.2f.pdf", lumi_scalefactors[i_lum]));

		/*
		// write values to file
		TString filepath = "../../data/11_1/";
		TString filename = Form("../../data/11_1/optimum_masswindow_lumiscalefactor%5.2f.txt", lumi_scalefactors[i_lum]);

		ofstream myfile;
		myfile.open(filename); //, ios::trunc | ios::out);
		myfile << setw(7) << "	" << setw(15) << "masswin" << setw(15) << "Z" << endl;
		myfile << setw(10) << setprecision(8) << "exp";
		myfile << setw(10) << setprecision(8) << masswindow_exp;
		myfile << setw(10) << setprecision(8) << Z_exp << endl;
		if(fabs(lumi_scalefactors[i_lum] - 1.00)<0.01){
		myfile << setw(10) << setprecision(8) << "obs";
		myfile << setw(10) << setprecision(8) << masswindow_obs;
		myfile << setw(10) << setprecision(8) << Z_obs << endl;
	}
	myfile.close();
	*/

}

//-- print optimum luminosity scalefactor to  screen
double lumi_sf, Z_sf;
int bin;
vector<double> lum_sf_greater_5sigma;
for (int i=0; i<n;i++){
	Z_sf = lumi_significances[i];
	printf("lumi_scalefactor = %5.2f, Z = %5.2f\n", lumi_scalefactors[i], Z_sf);
	if(Z_sf > 5){ lum_sf_greater_5sigma.push_back(lumi_scalefactors[i]); }
}

printf("L_sf = %5.2f gives Z_exp = %5.2f\n", lum_sf_greater_5sigma[0], h_lumi_significances->GetBinContent(h_lumi_significances->FindBin(lum_sf_greater_5sigma[0])-1));

// Print luminosity scale factors vs. optimum significance
TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
canvas2->SetLeftMargin(0.125);
canvas2->SetBottomMargin(0.125);
gStyle->SetOptTitle(1);

h_lumi_significances->SetAxisRange(1, lum_max-0.1);
h_lumi_significances->GetXaxis()->SetTitle("Luminosity scale factors");
h_lumi_significances->GetXaxis()->SetTitleSize(0.045);
h_lumi_significances->GetXaxis()->SetTitleOffset(1.0);
h_lumi_significances->GetYaxis()->SetTitle("Z_{exp} [#sigma]");
h_lumi_significances->GetYaxis()->SetTitleSize(0.045);
h_lumi_significances->GetYaxis()->SetTitleOffset(1.0);
canvas2->cd();
h_lumi_significances->Draw("l");

canvas2->Update();
TLine *MinimumBound = new TLine(1, h_lumi_significances->GetBinContent(h_lumi_significances->FindBin(lum_sf_greater_5sigma[0])), lum_max-0.1, h_lumi_significances->GetBinContent(h_lumi_significances->FindBin(lum_sf_greater_5sigma[0])));
MinimumBound->SetLineColor(kBlue); MinimumBound->SetLineWidth(2); MinimumBound->SetLineStyle(2);
MinimumBound->Draw("same");

//canvas2->Print("../../data/11_1/Luminosity_scalefactors_vs_optimum_significance.pdf");

return;

//================================
} // end Significance_Optimization()
//================================


void Get_ConfidenceLevels(double signal_sf, double lum_sf){
	TH1D *h_data, *h_bgr, *h_sig;
	TH1D *h_bgr_testStat, *h_sb_testStat, *h_sig_scaled;
	int Nbins;
	double s_scaled;

	h_sig = GetMassDistribution(125, lum_sf);
	h_bgr = GetMassDistribution(1, lum_sf);
	h_data = GetMassDistribution(2, lum_sf);

	// scale signal
	h_sig_scaled = (TH1D*) h_sig->Clone();
	for (int bin=1; bin<h_sig->GetNbinsX()+1;bin++){
		s_scaled = signal_sf*h_sig->GetBinContent(bin);
		h_sig_scaled->SetBinContent(bin, s_scaled);
	}

	TFile *file_bgr = new TFile("../../data/11_3/hist_test_statistic_bgr.root", "READ");
	TFile *file_sb = new TFile("../../data/11_3/hist_test_statistic_sb.root", "READ");

	h_sb_testStat = (TH1D*) file_sb->Get("h_test_statistic")->Clone();
	h_bgr_testStat = (TH1D*) file_bgr->Get("h_test_statistic")->Clone();

	Nbins = h_bgr_testStat->GetNbinsX();


	/*******************************************/
	/* Checking p-values for b-only hypothesis */
	/*******************************************/

	// confidence level (b-only)
	double median_bgr = Get_Quantiles(h_bgr_testStat)[2];
	double CLb_b = h_bgr_testStat->Integral(h_bgr_testStat->FindBin(median_bgr),Nbins) / h_bgr_testStat->Integral();
	printf("\nmedian_bgr = %5.2f\n", median_bgr);

	// confidence level (s+b)
	double median_sb = Get_Quantiles(h_sb_testStat)[2];
	double CLb_sb = h_bgr_testStat->Integral(h_bgr_testStat->FindBin(median_sb),Nbins) / h_bgr_testStat->Integral();
	printf("median_sb = %5.2f\n", median_sb);

	// confidence level (data)
	double t_obs = Get_TestStatistic(h_data, h_bgr, h_sig_scaled);
	double CLb_data = h_bgr_testStat->Integral(h_bgr_testStat->FindBin(t_obs),Nbins) / h_bgr_testStat->Integral();
	printf("t_obs = %5.2f\n", t_obs);

	// output to terminal
	printf("\n==================\n");
	printf("b-only hypothesis:\n");
	printf("==================\n");
	printf("------------------------------------------------\n");
	printf("type     CL_b     p-value     Z\n");
	printf("bgr:     %5.4f     %5.4f     %5.4f\n", CLb_b, 1-CLb_b, ROOT::Math::gaussian_quantile_c(1-CLb_b, 1));
	printf("s+b:     %5.4f     %5.4f     %5.4f\n", CLb_sb, 1-CLb_sb, ROOT::Math::gaussian_quantile_c(1-CLb_sb, 1));
	printf("data:     %5.4f     %5.4f     %5.4f\n", CLb_data, 1-CLb_data, ROOT::Math::gaussian_quantile_c(1-CLb_data, 1));
	printf("------------------------------------------------\n");

	/****************************************/
	/* Checking p-values for s+b hypothesis */
	/****************************************/

	// confidence level (b-only)
	double CLsb_b = h_sb_testStat->Integral(h_bgr_testStat->FindBin(median_bgr),Nbins) / h_sb_testStat->Integral();

	// confidence level (s+b)
	double CLsb_sb = h_sb_testStat->Integral(h_sb_testStat->FindBin(median_sb),Nbins) / h_sb_testStat->Integral();

	// confidence level (data)
	double CLsb_data = h_sb_testStat->Integral(h_sb_testStat->FindBin(t_obs),Nbins) / h_sb_testStat->Integral();

	// output to terminal
	// printf("\n==================\n");
	// printf("s+b hypothesis:\n");
	// printf("==================\n");
	// printf("------------------------------------------------\n");
	// printf("type     CL_sb     1-CLsb_b     Z		CL_s\n");
	// printf("bgr:     %5.4f     %5.4f     %5.4f		%5.4f\n", CLsb_b, 1-CLsb_b, ROOT::Math::gaussian_quantile_c(1-CLsb_b, 1), CLsb_b/((double)CLb_b));
	// printf("s+b:     %5.4f     %5.4f     %5.4f		%5.4f\n", CLsb_sb, 1-CLsb_sb, ROOT::Math::gaussian_quantile_c(1-CLsb_sb, 1), CLsb_sb/((double)CLb_sb));
	// printf("data:    %5.4f     %5.4f     %5.4f		%5.4f\n", CLsb_data, 1-CLsb_data, ROOT::Math::gaussian_quantile_c(1-CLsb_data, 1), CLsb_data/((double)CLb_data));
	// printf("------------------------------------------------\n");

	printf("\n==================\n");
	printf("s+b hypothesis v2:\n");
	printf("==================\n");
	printf("------------------------------------------------\n");
	printf("type		1-CL_b	CL_sb	Z	CL_s\n");
	printf("median bgr:     %5.4f	%5.4f	%5.4f	%5.4f\n", 1-CLb_b, CLsb_b, ROOT::Math::gaussian_quantile_c(1-CLb_b, 1), CLsb_b/((double)CLb_b));
	printf("median s+b:     %5.4f	%5.4f	%5.4f	%5.4f\n", 1-CLb_sb, CLsb_sb, ROOT::Math::gaussian_quantile_c(1-CLb_sb, 1), CLsb_sb/((double)CLb_sb));
	printf("median data:    %5.4f	%5.4f	%5.4f	%5.4f\n", 1-CLb_data, CLsb_data, ROOT::Math::gaussian_quantile_c(1-CLb_data, 1), CLsb_data/((double)CLb_data));
	printf("------------------------------------------------\n");

	return;
}


void SignalScaleFactors(){
	int N = 2;
	double lum, mu;
	double sf_min = 1.0, lum_min = 1.0;
	double sf_max = 3.0, lum_max = 7.0;
	double dSf = (sf_max-sf_min)/((double) N);
	double dLum = (lum_max-lum_min)/((double) N);

	vector<double> signal_sf;
	vector<double> elements;
	// set signal scalefactor vector
	// for (int i=0; i<N; i++){
	// 	signal_sf.push_back(sf_min + i*dSf);
	// 	cout << signal_sf[i] << endl;
	// }

	// testing
	// signal_sf.push_back(2.75);
	// signal_sf.push_back(3.50);
	signal_sf.push_back(1.00);

	// set luminosity scalefactor vector
	vector<double> lum_sf;
	lum_sf.push_back(1.50); lum_sf.push_back(2.00);	lum_sf.push_back(2.50);	lum_sf.push_back(3.00);
	lum_sf.push_back(1.00);
	// set signal scalefactor vector
	// for (int i=0; i<N; i++){
	// 	lum_sf.push_back(lum_min + i*dLum);
	// 	cout << lum_sf[i] << endl;
	// }

	if (signal_sf.size() > 1){ elements = signal_sf; }
	else if (lum_sf.size() > 1){ elements = lum_sf; }

	// find the confidence levels given signal scalefactors or luminosity scalefactors
	for (int i=0; i<elements.size(); i++){
		if (lum_sf.size() > 1){ lum = lum_sf[i]; }
		else{ lum = 1.0; }

		if (signal_sf.size() > 1){ mu = signal_sf[i]; }
		else{ mu = 1.0; }

		printf("\n----------------------------------\n");
		printf("Running signal sf = %5.2f\n", mu);
		printf("Running luminosity sf = %5.2f\n", lum);
		printf("----------------------------------\n");

		// run significance levels with new signal scalefactor
		Significance_Likelihood_ToyMC("bgr", mu, lum, 10000);
		Significance_Likelihood_ToyMC("sb", mu, lum, 10000);

		// get confidence levels based on output of significance likelihood
		Get_ConfidenceLevels();
	}
return;
}


//===========================
void SideBandFit(string type, int Irebin){
//===========================
	/* We will here find the expected and observed (although this should in practice not be done -> only used as a comparison to the previous task) significance after rescaling.
	Do this by:
	i) Computing multiple possible scalefactors (between Nbgr_rescaled = 10% of original Nbgr to Nbgr_rescaled = 300% greater than Nbgr) in the sideband regions, so
	as not to desturb the possible signal region.
	ii) The sideband region is in this case taken to be 150 <= m_ll <= 400.
	iii) Calculate which scalefactor gives the maximum log(L). This is the scalefactor that will be used in any further calculations.
	2) Reset the histogram for the background events, and obtain the rescaled histogram for the background events in its stead.
	3) Calculate the uncertainty of the background events (b ± sigma_b) using the rescaled background events vs. -2ln(L)
	4) Use this to calculate the expected and observed significances.
	*/


	printf("\n SideBandFit()\n\n");

	//-------------------------
	//-- [1] Prepare histograms
	//-------------------------
	TH1D *h_bgr, *h_data, *h_sig;
	h_sig = GetMassDistribution(125);//, 1.00);
	h_bgr  = GetMassDistribution(1, 1.00);
	h_data = GetMassDistribution(2, 1.00);

	//-- rebin histograms if necessary
	// h_bgr->Rebin(Irebin);
	// h_data->Rebin(Irebin);
	// printf(" INFO: Rebinning the histograms with a factor %d. Binwidth is now %5.3f GeV\n", Irebin, h_data->GetBinWidth(1));

	//-------------------------------
	//-- [2] Loop over scale factor
	//-------------------------------
	TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","Scalefactor SM",1000,0.1,3.1); // you'll want to put more here I guess
	TH1D *h_scalefactor_sig = new TH1D("h_scalefactor_sig","Scalefactor Higgs",1000,0.1,3.1); // you'll want to put more here I guess
	TH1D *h_loglikelihood_bgr = new TH1D("h_loglikelihood_bgr", "Log-likelihood SM", 1000, 0.1, 3.1);
	TH1D *h_loglikelihood_sig = new TH1D("h_loglikelihood_sig", "Log-likelihood Higgs", 1000, 0.1, 3.1);
	TH1D *h_loglikelihood_bgr_scaled = new TH1D("h_loglikelihood_bgr_scale", "Log-likelihood rescaled SM", 1000, 0.1, 3.1);
	TH1D *h_loglikelihood_sig_scaled = new TH1D("h_loglikelihood_sig_scale", "Log-likelihood rescaled Higgs", 1000, 0.1, 3.1);
	double scale_factor;//  = 1.00;		// scalefactors (respectively SM and Higgs);

	//---------------------------------------------
	// [2a] Loop 1: loop over scalefactors in alpha
	//---------------------------------------------
	for (int i_bin_sf = 0; i_bin_sf <= h_scalefactor_bgr->GetNbinsX()+1; i_bin_sf++){

		// get scale factor
		scale_factor = h_scalefactor_bgr->GetBinCenter(i_bin_sf);

		//-----------------------------------------------------------------------------------
		// [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
		//-----------------------------------------------------------------------------------
		double loglik_bgr = 1e-10, loglik_bgr_scale = 1e-10;
		double mass, Nbgr, Nbgr_rescale, Nsig;
		int Ndata;
		for (int i_bin = 0; i_bin <=  h_data->GetNbinsX()+1; i_bin++){
			Ndata = h_data->GetBinContent(i_bin);
			Nbgr = h_bgr->GetBinContent(i_bin);
			Nsig = h_sig->GetBinContent(i_bin);
			mass = h_data->GetBinCenter(i_bin);

			if (mass >= 150 && mass <= 400){
				// // background
				// Nbgr_rescale = scale_factor*Nbgr;
				// loglik_bgr += log(TMath::Poisson(Ndata, Nbgr));
				// loglik_bgr_scale += log(TMath::Poisson(Ndata, Nbgr_rescale));

				if (type == "bgr"){
					// printf("Ndata = %d, Nbgr = %5.2f\n", Ndata, Nbgr);
					Nbgr_rescale = scale_factor*Nbgr;
					loglik_bgr += log(TMath::Poisson(Ndata, Nbgr));
					// loglik_bgr_scale += log(TMath::Poisson(Ndata, Nbgr_rescale));
					// printf("loglik_bgr = %5.2f\n", loglik_bgr);
				}
				else if (type == "sig" && Nsig > 0){
					// printf("Ndata = %d, Nsig = %5.2f\n", Ndata, Nsig);
					Nbgr_rescale = scale_factor*Nsig;
					loglik_bgr += log(TMath::Poisson(Ndata, Nsig));
					// printf("loglik_bgr = %5.2f\n", loglik_bgr);
				}
				loglik_bgr_scale += log(TMath::Poisson(Ndata, Nbgr_rescale));
			}


		}
		h_loglikelihood_bgr->SetBinContent(i_bin_sf, loglik_bgr);
		h_loglikelihood_bgr_scaled->SetBinContent(i_bin_sf, loglik_bgr_scale);
		h_scalefactor_bgr->SetBinContent(i_bin_sf,-2.*loglik_bgr_scale);
	}


	//----------------------------------------------------
	// [3] Interpret the -2Log (Likelihood distribution)
	//----------------------------------------------------
	/* The optimal scalefactor (alpha, mu) is the scalefactor which gives the largest log likelihood (maximum of log(L)), i.e.
	the minimum of -2log(L). */

	/*************************************************/
	/* Find optimal scale factor for SM (background) */
	/*************************************************/
	TH1D *h_scalefactor_bgr_rescaled  = (TH1D*) h_scalefactor_bgr->Clone("h_scalefactor_bgr_rescaled");
	h_scalefactor_bgr_rescaled->Reset();

	// Find maximum of log(L) of background
	int bin_max_bgr = h_loglikelihood_bgr->GetMaximumBin();
	double loglik_max_bgr = h_loglikelihood_bgr->GetBinContent(bin_max_bgr);

	// Find maximum of log(L) of rescaled background
	int bin_max_scale_bgr = h_loglikelihood_bgr_scaled->GetMaximumBin();
	double alpha = h_loglikelihood_bgr_scaled->GetBinCenter(bin_max_scale_bgr);
	double loglik_max_scale_bgr = h_loglikelihood_bgr_scaled->GetBinContent(bin_max_scale_bgr);

	// Rescale and find +/- 1 sigma errors by plotting -2log(L) vs. alpha with the rescaled plot for the background events
	double d_logL;
	vector<double> dAlpha;
	for (int i_bin = 0; i_bin <= h_loglikelihood_bgr->GetNbinsX(); i_bin++){
		d_logL = h_loglikelihood_bgr_scaled->GetBinContent(i_bin) - loglik_max_scale_bgr;
		if(d_logL >= -0.5*1.05 && d_logL <= -0.5*0.95){
			dAlpha.push_back(-(alpha - h_loglikelihood_bgr->GetBinCenter(i_bin)));
		}
	}

	//----------------------------------------------------
	// [4] Calculate rescaled events with uncertainty
	//----------------------------------------------------

	/***************************************/
	/* Rescaled events for SM (background) */
	/***************************************/

	int Ndata_win;
	int bin1, bin2;
	double Nbgr_win, Nbgr_win_rescale, dNbgr_rescale_central, Nsig_win;
	double x_min, x_max, masswindow_fullwidth;
	TAxis *xAxis;
	vector<double> dNbgr_win_rescale;

	masswindow_fullwidth = 7.15;	//GeV
	xAxis = h_bgr->GetXaxis();
	x_min = 125.-0.5*masswindow_fullwidth;
	x_max = 125+0.5*masswindow_fullwidth;
	bin1 = xAxis->FindBin(x_min);
	bin2 = xAxis->FindBin(x_max);

	Ndata_win = h_data->Integral(bin1, bin2);
	Nbgr_win = h_bgr->Integral(bin1, bin2);
	Nsig_win = h_sig->Integral(bin1, bin2);

	// calculate rescaled b with uncertainty (b +/- db)
	if (type=="bgr"){
		Nbgr_win_rescale = Nbgr_win*alpha;
		dNbgr_win_rescale.push_back(Nbgr_win*dAlpha[0]); dNbgr_win_rescale.push_back(Nbgr_win*dAlpha[1]);
	}
	else if(type=="sig"){
		Nbgr_win_rescale = Nsig_win*alpha;
		dNbgr_win_rescale.push_back(Nsig_win*dAlpha[0]); dNbgr_win_rescale.push_back(Nsig_win*dAlpha[1]);
	}

	// Finding the gaussian equivalence of the uncertainty (not poisson distributed around alpha or mu, but rather gaussian)
	dNbgr_rescale_central = 0.5*(dNbgr_win_rescale[0]+dNbgr_win_rescale[1]);
	dNbgr_rescale_central = 0.5*(dNbgr_win_rescale[0]+dNbgr_win_rescale[1]);

	/* Computing the expected and observed significances, based on the new background estimate */
	double Z_exp = ExpectedSignificance_ToyMC(Nbgr_win_rescale, Nsig_win, Ndata_win, 0, 1e4, false);
	double Z_obs = ExpectedSignificance_ToyMC(Nbgr_win_rescale, Nsig_win, Ndata_win, 0, 1e4,  true);

	/* Computing the expected and observed significances, based on the new background estimate with uncertainty */
	double Z_exp_unc = ExpectedSignificance_ToyMC(Nbgr_win_rescale, Nsig_win, Ndata_win, dNbgr_rescale_central, 1e4, false);
	double Z_obs_unc = ExpectedSignificance_ToyMC(Nbgr_win_rescale, Nsig_win, Ndata_win, dNbgr_rescale_central, 1e4,  true);


	/* Finding the test statistic (log ratio) (Part 3) */
	double test_statistic = Get_TestStatistic(h_data, h_bgr, h_sig);

	// print summary to screen
	printf("\n----------------------------------------------------------\n");
	printf("alpha: %5.2f + (%5.2f), %5.2f + (%5.2f)\n", alpha, dAlpha[0], alpha, dAlpha[1]);
	printf("Masswindow = %5.2f:\n", masswindow_fullwidth);
	if (type=="bgr"){
		printf("Nsig = %5.2f, Ndata = %d\n", Nsig_win, Ndata_win);
		printf("Nbgr = %5.2f, Nbgr scaled = %5.2f +(%5.2f), %5.2f + (%5.2f)\n", Nbgr_win, Nbgr_win_rescale, dNbgr_win_rescale[0], Nbgr_win_rescale, dNbgr_win_rescale[1]);
		printf("\nExpected significance (toy MC) = %5.2f\n", Z_exp);
		printf("Observed significance (toy MC) = %5.2f\n", Z_obs);
		printf("\nExpected significance with background uncertainty (toy MC) = %5.2f\n", Z_exp);
		printf("Observed significance with background uncertainty (toy MC) = %5.2f\n", Z_obs);
	}
	else if (type=="sig"){
		printf("Nbgr = %5.2f, Ndata = %d\n", Nbgr_win, Ndata_win);
		printf("Nsig = %5.2f, Nsig scaled = %5.2f +(%5.2f), %5.2f + (%5.2f)\n", Nbgr_win, Nbgr_win_rescale, dNbgr_win_rescale[0], Nbgr_win_rescale, dNbgr_win_rescale[1]);
	}
	printf("\n Test statistic = %5.4f\n", test_statistic);
	printf("----------------------------------------------------------\n");

	//	h_scalefactor_bgr->Print("all");

	/*
	//---------------------------------------------------
	// PLOTTING:
	// print scalefactor vs. -2log(L) to pdf
	TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
	canvas1->SetLeftMargin(0.125);
	canvas1->SetBottomMargin(0.125);
	gStyle->SetOptTitle(1);
	//	gStyle->SetTitleFont(0.01, "t");

	h_scalefactor_bgr->SetAxisRange(1.0,1.2, "X");
	h_scalefactor_bgr->SetAxisRange(h_scalefactor_bgr->GetMinimum(h_scalefactor_bgr->GetBin(alpha_optimal))-2,1367.5, "Y");
	h_scalefactor_bgr->GetXaxis()->SetTitle("#alpha");
	h_scalefactor_bgr->GetXaxis()->SetTitleSize(0.04);
	h_scalefactor_bgr->GetXaxis()->SetTitleOffset(1.0);
	h_scalefactor_bgr->GetYaxis()->SetTitle("-2log(L)");
	h_scalefactor_bgr->GetYaxis()->SetTitleOffset(1.3);
	h_scalefactor_bgr->GetYaxis()->SetTitleSize(0.04);
	h_scalefactor_bgr->SetNameTitle("", "Background scalefactor");
	// canvas1->SetGrid();
	canvas1->cd();

	h_scalefactor_bgr->Draw("l");

	canvas1->Update();
	TLine *sigma1_min = new TLine(alpha_optimal+dAlpha[0], h_scalefactor_bgr->GetMinimum(h_scalefactor_bgr->GetBin(alpha_optimal+dAlpha[0])), alpha_optimal+dAlpha[0],h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal+dAlpha[0])));
	TLine *sigma1_max = new TLine(alpha_optimal+dAlpha[1], h_scalefactor_bgr->GetMinimum(h_scalefactor_bgr->GetBin(alpha_optimal+dAlpha[1])), alpha_optimal+dAlpha[1],h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal+dAlpha[1])));
	TLine *alpha_line = new TLine(alpha_optimal, h_scalefactor_bgr->GetMinimum(h_scalefactor_bgr->GetBin(alpha_optimal)), alpha_optimal,h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal)));
	TLine *max 		  = new TLine(1.0, h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal+dAlpha[0])), 1.2, h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal+dAlpha[0])));
	TLine *min 		  = new TLine(1.0, h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal)), 1.2, h_scalefactor_bgr->GetBinContent(h_scalefactor_bgr->FindBin(alpha_optimal)));
	sigma1_min->SetLineColor(616-8); sigma1_min->SetLineWidth(3);	//kMagenta
	sigma1_max->SetLineColor(616-8); sigma1_max->SetLineWidth(3);	//kMagenta
	alpha_line->SetLineColor(kRed); alpha_line->SetLineWidth(3);
	min->SetLineColor(kGray+2); min->SetLineWidth(1); min->SetLineStyle(2);
	max->SetLineColor(kGray+2); max->SetLineWidth(1); max->SetLineStyle(2);
	sigma1_min->Draw();
	sigma1_max->Draw();
	alpha_line->Draw();
	min->Draw();
	max->Draw();

	//	canvas1->Print("../../data/11_2/SideBand_fit.pdf");
	//---------------------------------------------------
	*/
	return;

}


//=========================================================================================
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig){
	//=========================================================================================

	double test_statistic = 0.;

	double Nbgr, Nsig;
	int Ndata;
	double loglik_bgr = 0.0, loglik_sig_plus_bgr = 0.0;

	for (int i_bin = 1; i_bin <= h_mass_dataset->GetNbinsX()+1; i_bin++){
		Ndata = h_mass_dataset->GetBinContent(i_bin);
		Nbgr = h_template_bgr->GetBinContent(i_bin);
		Nsig = h_template_sig->GetBinContent(i_bin);

		// likelihood for the mu=0 (b-only) scenario
		loglik_bgr += log(TMath::Poisson(Ndata, Nbgr));

		// likelihood for the mu=1 (s+b) scenario
		loglik_sig_plus_bgr += log(TMath::Poisson(Ndata, Nbgr+Nsig));
	}

	// printf("\nL(x|H_0) = %5.4f\n", loglik_bgr);
	// printf("L(x|H_1) = %5.4f", loglik_sig_plus_bgr);

	// compute test statistic
	test_statistic = 2*(loglik_bgr - loglik_sig_plus_bgr);

	return test_statistic;
}



//=========================================================================================
TH1D * TestStatistic_DistributionMC(TH1D *h_mass_template, TH1D *h_bgr, TH1D *h_sig, int Ncycle){
	//=========================================================================================
	TH1D *h_toy_dataset;
	double test_statistic;
	// generate histogram with test statistics
	TH1D *h_test_statistic = new TH1D("h_test_statistic", "", h_mass_template->GetNbinsX(), -40, 40);

	for (int cycle=0; cycle<Ncycle; cycle++){
		// generate toy dataset for given mass distribution
		h_toy_dataset = GenerateToyDataSet(h_mass_template);
		// generate test statistic on toy data set
		test_statistic = Get_TestStatistic(h_toy_dataset, h_bgr, h_sig);
		h_test_statistic->Fill(test_statistic);

		//		printf("t-test = %5.2f\n", test_statistic);
	}
	// h_test_statistic->Print("all");

	return h_test_statistic;
}


//===============================================
TH1D * GenerateToyDataSet(TH1D *h_mass_template){
	//===============================================
	//-------------------------------------------------
	// Goal: Generate Toy data set from input histogram
	// How:  dumb way -> draw random Poisson in each bin
	//--------------------------------------------------
	TRandom *RNG = new TRandom(0);
	TH1D *h_mass_toydataset = (TH1D*) h_mass_template->Clone("h_mass_toydataset");

	//-- Create new histogram for the data-set
	h_mass_toydataset->Reset();

	int NbinsX = h_mass_template->GetNbinsX();
	double mu, N_bin;


	//-- Loop over bins and draw Poisson number of event in each bin
	for (int i_bin =0; i_bin <= NbinsX+1; i_bin++){
		// generate mean of Poisson equal to bin content at central value
		mu = h_mass_template->GetBinContent(i_bin);
		// generate random Poisson number of events around mu
		N_bin = RNG->Poisson(mu);
		//		printf("N_bin = %5.9f\n", N_bin);
		// fill toy histogram
		h_mass_toydataset->SetBinContent(i_bin, N_bin);
	}

	// return histogram of toy data-set
	return h_mass_toydataset;
}

void Plot_ScaleFactors(int NBins, int irebin=10){
	// get histograms
	TH1D *h_data, *h_bgr, *h_sig;

	h_sig = GetMassDistribution(125);
	h_bgr = GetMassDistribution(1);
	h_data = GetMassDistribution(2);

	h_bgr->Rebin(irebin);
	h_sig->Rebin(irebin);
	h_data->Rebin(irebin);

	TH2D *h_scalefactor_L = new TH2D("h_scalefactor", "Not to be used other than for sigma", NBins, 0.5, 2., NBins, 0., 5.);
	TH2D *h_scalefactor = new TH2D("h_scalefactor", "Scale factors", NBins, 0.5, 2., NBins, 0., 5.);

	double loglik;
	double sf_bgr, sf_sig;
	double Nsb_bin, NObs_bin, Nsb_rescaled;

	for(int i=0; i<h_scalefactor->GetNbinsX()+1; i++){
		for(int j=0; j<h_scalefactor->GetNbinsY()+1; j++){

			sf_bgr = h_scalefactor->GetXaxis()->GetBinCenter(i);
			sf_sig = h_scalefactor->GetYaxis()->GetBinCenter(j);

			loglik = 0.;

			for(int bin_data=0; bin_data < h_data->GetNbinsX()+1; bin_data++){
				Nsb_bin = h_bgr->GetBinContent(bin_data) + h_sig->GetBinContent(bin_data);
				NObs_bin = h_data->GetBinContent(bin_data);

				Nsb_rescaled = sf_bgr*h_bgr->GetBinContent(bin_data) + sf_sig*h_sig->GetBinContent(bin_data);

				if(Nsb_rescaled>0){ loglik += TMath::Log(TMath::Poisson(NObs_bin, Nsb_rescaled)); }
			}
			h_scalefactor_L->SetBinContent(i, j, loglik);
			h_scalefactor->SetBinContent(i, j, -2*loglik);
		}
	}

	int x, y, z, minimumBin;
	double alpha, mu, min_logL;
	minimumBin = h_scalefactor->GetMinimumBin();		//get the minimum bin in plot
	h_scalefactor->GetBinXYZ(minimumBin, x, y, z);		// from the minimum bin, find the equivalent bins along the x-, y-, and z-axis
	mu = h_scalefactor->GetYaxis()->GetBinCenter(y);
	alpha = h_scalefactor->GetXaxis()->GetBinCenter(x);
	min_logL = h_scalefactor->GetBinContent(x, y);

	// rescale histogram
	for(int i=0; i<h_scalefactor->GetNbinsX()+1; i++){
		for(int j=0; j<h_scalefactor->GetNbinsY()+1; j++){
			h_scalefactor->SetBinContent(i, j, h_scalefactor->GetBinContent(i, j)-min_logL);
		}
	}

	// obtain +/-1sigma errors
	TH2D *h_sigma = (TH2D*) h_scalefactor->Clone();
	h_sigma->Reset();
	for(int i=0; i<h_sigma->GetNbinsX()+1; i++){
		for(int j=0; j<h_sigma->GetNbinsY()+1; j++){
			if(h_scalefactor->GetBinContent(i, j) <= 1.00){
				h_sigma->SetBinContent(i, j, 1.);
				if (h_scalefactor->GetXaxis()->GetBinCenter(i) == alpha){
					printf("dMu = %5.2f\n", h_sigma->GetYaxis()->GetBinCenter(j)-mu);
				}
				if (h_scalefactor->GetYaxis()->GetBinCenter(j) == mu){
					// printf("dAlpha = %5.2f\n", h_sigma->GetXaxis()->GetBinCenter(i)-alpha);
				}
			}
			// if(0.99<=h_scalefactor->GetBinContent(i, j)  && h_scalefactor->GetBinContent(i, j) <= 1.00){
			// 	printf("dAlpha: %5.2f\n", h_scalefactor->GetXaxis()->GetBinCenter(i)-alpha);
			// 	printf("dMu: %5.2f\n", h_scalefactor->GetYaxis()->GetBinCenter(j)-mu);
			// }
		}
	}

	double maximumBin = h_sigma->GetMaximumBin();
	minimumBin = h_sigma->GetMinimumBin();
	h_sigma->GetBinXYZ(maximumBin, x, y, z);
	double mu1 = h_sigma->GetYaxis()->GetBinCenter(y);
//	double alpha1 = h_
	h_sigma->GetBinXYZ(minimumBin, x, y, z);
	double mu2 = h_sigma->GetYaxis()->GetBinCenter(y);

	printf("mu = %5.2f+(%5.2f), %5.2f+(%5.2f)", mu, mu1 - mu, mu, mu2 - mu);


//	printf("\n Test: mu = %5.2f+(%5.2f), %5.2f+(%5.2f)\n", mu, h_sigma->GetYaxis()->GetMinimumBin() - mu, mu, h_sigma->GetYaxis()->GetMaximum());

	//print to terminal
	printf("\n--------------------------\n");
	printf("x = %d, y = %d, z = %d\n", x, y, z);
	printf("Minimum -2logL: %5.3f\n", min_logL);
	printf("alpha = %5.3f, mu = %5.3f\n", alpha, mu);
	printf("----------------------------\n");


	TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
	canvas1->SetLeftMargin(0.125);
	canvas1->SetBottomMargin(0.125);
	canvas1->cd();
	gStyle->SetOptTitle(1);
	gStyle->SetTitleFontSize(0.001);
	gStyle->SetTitleFont(42, "t");

	TMarker *sf = new TMarker(alpha, mu, 29);
	sf->SetMarkerSize(2);
	sf->SetMarkerColor(kWhite);

	h_scalefactor->GetXaxis()->SetTitle("#alpha_{bgr}");
	h_scalefactor->GetYaxis()->SetTitle("#mu_{s}");
	h_scalefactor->GetZaxis()->SetTitle("-2log(L)");
	h_scalefactor->GetXaxis()->SetTitleSize(0.045);
	h_scalefactor->GetYaxis()->SetTitleSize(0.045);
	h_scalefactor->GetXaxis()->SetTitleOffset(1.0);
	h_scalefactor->GetYaxis()->SetTitleOffset(1.0);
	h_sigma->SetLineColor(632+1);	//kRed

	h_scalefactor->Draw("colz");
	sf->Draw("SAME");
	h_sigma->Draw("SAME cont3");

	canvas1->Print("../../data/11_3/optimal_scale_factors_-2logL.pdf");

	return;
}















//================================================
// S O M E   U S E F U L    F U N C T I O N S
//================================================



//====================================================
double IntegratePoissonFromRight(double mu, int N_obs){
	//====================================================

	// --------------------------------------------------------------------
	// Compute p-value for case zero background uncertainty, i.e.
	//         just integrate Poisson from the right from N_obs to infinity
	// --------------------------------------------------------------------

	double integral = 1.;
	for(int i_obs = 0; i_obs < N_obs; i_obs++){
		integral -= TMath::Poisson(i_obs,mu);
	}

	return integral;

} // end IntegratePoissonFromRight()

//========================================================
double IntegrateFromRight(TH1D * h_X_bgr, double X_value){
	//========================================================
	// --------------------------------------------------------------------
	// Compute p-value: integrate number of events from X_value to infinity
	// --------------------------------------------------------------------

	//-- Integrate distributions
	int Nbins = h_X_bgr->GetNbinsX();
	int X_bin = h_X_bgr->FindBin(X_value);

	//-- Compute integral from X-value to infinity
	double pvalue = h_X_bgr->Integral(X_bin,Nbins) / h_X_bgr->Integral();

	return pvalue;

} // end IntegrateFrom Right()

//===============================================================================
double ExpectedSignificance_ToyMC(double Nbgr, double Nsig, int Ndata, double dNbgr, int N, bool obs){
	//===============================================================================
	// --------------------------------------------------------------------
	// Compute significance: Generate toy MC datasets to compute the significance.
	// Steps:
	// 1) Draw a random number background events and s+b events based on input
	// (maximum value)
	// 2) Calculate the cdf for the Poisson function (integrate from right) in
	// order to obtain the p-value
	// 3) Calculate significance from gaussian quantiles (Central limit theorem)
	// --------------------------------------------------------------------

	double Nb_rand, Nbs_rand;
	double p_exp=0.0, p_obs=0.0;
	double Z_exp=0.0, Z_obs=0.0;

	// initialize seed
	TRandom *RNG = new TRandom();

	for (int i=0; i<N; i++){
		// obtain random (Poisson) generated number for s+b, and gaussian for b-only (central limit theorem)
		Nb_rand = RNG->Gaus(Nbgr, dNbgr);
		Nbs_rand = Nb_rand + RNG->Poisson(Nsig);

		// calculate p-values
		p_exp = IntegratePoissonFromRight(Nb_rand, Nbs_rand);
		p_obs = IntegratePoissonFromRight(Nb_rand, Ndata);

		// disregard any non-sensible values
		if (p_exp < 0 || p_exp >= 1){ continue; }
		if (p_obs < 0 || p_obs >= 1){ continue; }

		// calculate significances
		Z_exp += ROOT::Math::gaussian_quantile_c(p_exp, 1);
		Z_obs += ROOT::Math::gaussian_quantile_c(p_obs, 1);
	}

	Z_exp = Z_exp/((double) N);
	Z_obs = Z_obs/((double) N);

	if( obs == true){ return Z_obs; }
	else { return Z_exp; }
}

//===========================================================================
void Significance_Likelihood_ToyMC(string mass_distribution_type, double signal_sf=1.0, double lum_sf=1.0, int Ncycle=10000){
	//===========================================================================
	// Generate the test statistic based on toy MC datasets
	TH1D *h_mass_template, *h_bgr, *h_sig, *h_toy_dataset;
	TH1D *h_test_statistic;
	double test_statistic, s_scaled;

	// get signal and background histograms
	h_bgr = GetMassDistribution(1, lum_sf);
	h_sig = GetMassDistribution(125, lum_sf);

	printf("Number MC cycles: %d\n", Ncycle);

	// account for potential signal scalefactors
	TH1D *h_sig_scaled = (TH1D*) h_sig->Clone();
	for (int bin=0; bin<h_sig->GetNbinsX(); bin++){
		s_scaled = signal_sf*h_sig->GetBinContent(bin);
		h_sig_scaled->SetBinContent(bin, s_scaled);
	}

	// determine which mass distribution to base toy data on
	if (mass_distribution_type=="bgr"){ h_mass_template = (TH1D*) h_bgr->Clone();
		h_test_statistic = TestStatistic_DistributionMC(h_mass_template, h_bgr, h_sig_scaled, Ncycle);

		// print to file
		TFile *f = new TFile(("../../data/11_3/hist_test_statistic_bgr.root"), "RECREATE");
		h_test_statistic->Write();
		f->Close();

		printf("hist_test_statistic_bgr.root file created\n");
	}

	else if (mass_distribution_type=="sb"){
		h_mass_template = (TH1D*) h_bgr->Clone("h_mass_template");
		for (int i_bin=1; i_bin<(h_mass_template->GetNbinsX()+1); i_bin++){
			h_mass_template->SetBinContent(i_bin, h_bgr->GetBinContent(i_bin) + h_sig_scaled->GetBinContent(i_bin));
		}
		// printf("so far so good");
		h_test_statistic = TestStatistic_DistributionMC(h_mass_template, h_bgr, h_sig_scaled, Ncycle);
		// printf("so far so good");

		TFile *f = new TFile("../../data/11_3/hist_test_statistic_sb.root", "RECREATE");
		h_test_statistic->Write();
		f->Close();

		printf("hist_test_statistic_sb.root file created\n");
	}

	else if (mass_distribution_type=="data"){
		h_mass_template = GetMassDistribution(2);
		double test_statistic = Get_TestStatistic(h_mass_template, h_bgr, h_sig_scaled);
		printf("Test-statistic data = %5.2f\n", test_statistic);
	}

	return;
}



//=========================================
vector<double> Get_Quantiles( TH1D* hist ){
	//=========================================
	// Quantiles returns a vector<double> with 5 entries.
	// Entries 0 and 4 are the values on the histogram x-axis
	// so that 95% of the content lies between these values.
	// Entries 1 and 3 bracket 68% in the same way.
	// Entry 2 is the median of the histogram.

	// define quantiles
	double fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // 15.8655 %
	double fraction_2sigma = ROOT::Math::gaussian_cdf(-2.,1.,0.); //  2.2750 %
	double probs[5] = {fraction_2sigma, fraction_1sigma, 0.50, 1.00-fraction_1sigma, 1.00-fraction_2sigma };

	// output of the quantiles
	double Xvalues[5];

	// extract quantiles
	hist->GetQuantiles( 5, Xvalues, probs );

	vector<double> Xvalues_output(5);
	for (int i=0; i<5; i++)
	{
		Xvalues_output[i] = Xvalues[i];
	}

	return Xvalues_output;
}





//=======================================================================================================================
void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
	Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color)
	//=======================================================================================================================
	{
		Int_t txt_align = 12;
		if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left
		if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
		if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center

		TLatex* t1 = new TLatex( txt_x, txt_y, txt);
		if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
		t1->SetTextSize(txt_size);
		t1->SetTextAlign(txt_align);
		t1->SetTextAngle(txt_angle);
		t1->SetTextColor(txt_color);
		t1->Draw();

	} // end AddText()


	/*
	//=======================================================================================================================
	void WriteToFile(filepath, filename):
	//=======================================================================================================================
	// Write values to file
	return;
	*/
