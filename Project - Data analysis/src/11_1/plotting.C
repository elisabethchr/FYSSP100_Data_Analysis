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

// declaration of functions
void Plot_TestStatistic(bool plot_tData, double tData=-11.53);
vector<double> Get_Quantiles( TH1D* hist );


void Plot_TestStatistic(bool plot_tData, double tData){
    // [1] Get histograms from file
    TH1D *h_data = 0, *h_bgr = 0, *h_sb = 0;
    //    TDirectory* dir = gDirectory;
    TFile *file_bgr = new TFile("../../data/11_3/hist_test_statistic_bgr.root", "READ");
    TFile *file_sb = new TFile("../../data/11_3/hist_test_statistic_sb.root", "READ");

    h_sb = (TH1D*) file_sb->Get("h_test_statistic")->Clone();
    h_bgr = (TH1D*) file_bgr->Get("h_test_statistic")->Clone();

    h_sb->Rebin(5);     // merge five bins into one
    h_bgr->Rebin(5);    // merge five bins into one

    // [2] Get the +/- 1sigma and 2sigma errors on the b-only distribution:
    TH1D *h_bgr_sigma1 = (TH1D*) h_bgr->Clone();
    TH1D *h_bgr_sigma2_min = (TH1D*) h_bgr->Clone();
    TH1D *h_bgr_sigma2_max = (TH1D*) h_bgr->Clone();

    // (range of which 68% of the data fall in)
    double minSigma_68 = Get_Quantiles(h_bgr)[1];
    double maxSigma_68 = Get_Quantiles(h_bgr)[3];
    // (range of which 95% of the data fall in)
    double minSigma_95 = Get_Quantiles(h_bgr)[0];
    double maxSigma_95 = Get_Quantiles(h_bgr)[4];

    h_bgr_sigma1->GetXaxis()->SetRangeUser(minSigma_68, maxSigma_68);
    h_bgr_sigma2_min->GetXaxis()->SetRangeUser(maxSigma_68, maxSigma_95);
    h_bgr_sigma2_max->GetXaxis()->SetRangeUser(minSigma_95, minSigma_68);

    // Get the +/- 1sigma and 2sigma errors on the s+b distribution:
    TH1D *h_sb_sigma1 = (TH1D*) h_sb->Clone();
    TH1D *h_sb_sigma2_min = (TH1D*) h_sb->Clone();
    TH1D *h_sb_sigma2_max = (TH1D*) h_sb->Clone();

    // (range of which 68% of the data fall in)
    minSigma_68 = Get_Quantiles(h_sb)[1];
    maxSigma_68 = Get_Quantiles(h_sb)[3];
    // (range of which 95% of the data fall in)
    minSigma_95 = Get_Quantiles(h_sb)[0];
    maxSigma_95 = Get_Quantiles(h_sb)[4];

    h_sb_sigma1->GetXaxis()->SetRangeUser(minSigma_68, maxSigma_68);
    h_sb_sigma2_min->GetXaxis()->SetRangeUser(maxSigma_68, maxSigma_95);
    h_sb_sigma2_max->GetXaxis()->SetRangeUser(minSigma_95, minSigma_68);


    // [3] Prepare canvas for histograms
    TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
    canvas1->SetLeftMargin(0.125) ;
    canvas1->SetBottomMargin(0.125);
    canvas1->SetGrid();
    canvas1->cd();
    gStyle->SetOptTitle(1);

    h_bgr->SetLineColor(kBlack);
    h_bgr->SetLineWidth(2);
    h_sb->SetLineColor(kBlue);
    h_sb->SetLineWidth(2);
    h_bgr->SetNameTitle("", "Toy MC - test statistic");
    h_bgr->GetXaxis()->SetTitleSize(0.040);
    h_bgr->GetXaxis()->SetTitle("Test statistic");
    h_bgr->GetXaxis()->SetTitleOffset(1.0);
    h_bgr->GetYaxis()->SetTitle("MC cycles");
    h_bgr->GetYaxis()->SetTitleSize(0.040);
    h_bgr->GetYaxis()->SetTitleOffset(1.0);
    h_bgr_sigma1->SetFillColor(kGreen);
    h_sb_sigma1->SetFillColor(kGreen);
    h_bgr_sigma2_min->SetFillColor(kYellow);
    h_bgr_sigma2_max->SetFillColor(kYellow);
    h_sb_sigma2_min->SetFillColor(kYellow);
    h_sb_sigma2_max->SetFillColor(kYellow);

    h_bgr->Draw("l");
    h_bgr_sigma1->Draw("same");
    h_bgr_sigma2_min->Draw("same");
    h_bgr_sigma2_max->Draw("same");
    h_sb->Draw("l same");

    // Add test-statistic line for data
    TArrow *line_tData;
    if(plot_tData == 1){
        canvas1->Update();
        line_tData = new TArrow(tData, 35, tData, h_bgr->GetMaximum(), 0.02, "<");
//        line_tData->SetDefaultArrowSize(0.001);
        line_tData->SetLineColor(kRed); //kYellow
        line_tData->SetLineWidth(2); //line_tData->SetLineStyle(2);
        line_tData->Draw();
    }

    // Add legend
    TLegend *leg = new TLegend(0.15, 0.75,0.35,0.85);
    leg->SetBorderSize(0); leg->SetFillColor(0);// leg->SetLineColor(kBlack);
    TLegendEntry *leg1a = leg->AddEntry(h_bgr, " b-only", "f");  leg1a->SetTextSize(0.04);
    TLegendEntry *leg1b = leg->AddEntry(h_sb, " s+b" , "f");  leg1b->SetTextSize(0.04);
    TLegendEntry *leg2a = leg->AddEntry(line_tData, "data", "l"); leg2a->SetTextSize(0.04);
    leg->Draw();

    // Save to pdf
//    canvas1->Print("../../data/11_3/test_statistic_distribution.pdf");

    // Calculate median of histograms
    double x, q, median;
    q = 0.5;
    h_bgr->ComputeIntegral();
    h_bgr->GetQuantiles(1, &x, &q);
    median = x;

    printf("Median test-statistic for b-only = %5.2f\n", median);

    h_sb->ComputeIntegral();
    h_sb->GetQuantiles(1, &x, &q);
    median = x;

    printf("Median test-statistic for s+b = %5.2f\n", median);
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
