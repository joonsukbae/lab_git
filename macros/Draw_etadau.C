#include "BSHelper.cxx"
#include "Filipad2.h"

enum {
  kECbegin = 0,
  kDATA = 1,
  kINEL,
  kECend
};
enum {
  kTrigbegin = 0,
  kMBAND = 1,
  kTrigend
};
enum {
  kSpeciesbegin = 0,
  kK0short = 1,
  kLambda,
  kAntilambda,
  kSpeciesend
};
enum {
  kSignbegin = 0,
  kPositive = 1,
  kNegative,
  kSignend
};
enum {
  kStepbegin = 0,
  kAll = 1,
  kBasiccut,
  kMasscut,
  kStepend
};

const char *MCFile = "~/cernbox/O2Mult/mc/AnalysisResults_mc_pbpb536_test_0201.root";
const char *DataFile = "~/cernbox/O2Mult/data/AnalysisResults_data_pbpb536_test_0201.root";
const char *MCdir = "multiplicity-counter";
const char *Datadir = "multiplicity-counter";
const char *MCDataset = "LHC22k3b2";
const char *DataDataset = "LHC22s_pass4";

Int_t n=0;
Int_t nn=0;
const char *mstring = "d#it{N}_{ch}/d#it{#eta}";
const char *avmstring = "LT d#it{N}_{ch}/d#it{#eta} #GT";
const char *avmstringincl = "#LT d#it{N}_{ch}/d#it{#eta} #GT / #LT d#it{N}_{ch}/d#it{#eta} #GT_{Inclusive}";
const char *mstringincl = "(d#it{N}_{ch}/d#it{#eta})/(d#it{N}_{ch}/d#it{#eta})_{Inclusive}";

void setpad(TVirtualPad *pad)
{
  pad->SetTopMargin(0.02);
  pad->SetLeftMargin(0.13);
  pad->SetRightMargin(0.2);
  pad->SetBottomMargin(0.15);
  pad->SetName(Form("c%d", ++n));
}

void hset(TH1 &hid, TString xtit = "", TString ytit = "",
          double titoffx = 0.9, double titoffy = 1.2,
          double titsizex = 0.06, double titsizey = 0.06,
          double labeloffx = 0.01, double labeloffy = 0.001,
          double labelsizex = 0.05, double labelsizey = 0.05,
          int divx = 510, int divy = 510)
{
	hid.SetStats(0);

	hid.GetXaxis()->CenterTitle(1);
	hid.GetYaxis()->CenterTitle(1);

	hid.GetXaxis()->SetTitleOffset(titoffx);
	hid.GetYaxis()->SetTitleOffset(titoffy);

	hid.GetXaxis()->SetTitleSize(titsizex);
	hid.GetYaxis()->SetTitleSize(titsizey);

	hid.GetXaxis()->SetLabelOffset(labeloffx);
	hid.GetYaxis()->SetLabelOffset(labeloffy);

	hid.GetXaxis()->SetLabelSize(labelsizex);
	hid.GetYaxis()->SetLabelSize(labelsizey);

	hid.GetXaxis()->SetNdivisions(divx);
	hid.GetYaxis()->SetNdivisions(divy);

	hid.GetXaxis()->SetTitle(xtit);
	hid.GetYaxis()->SetTitle(ytit);
}

void Draw_etadau() {

// MC File open and find objects
    auto fmc = TFile::Open(MCFile, "open");
    auto dir = (TDirectory *) fmc->Get(MCdir);
    dir->cd();

    // reconstructed zvtx (MC)
    auto HRZ = (THnSparse *) gROOT->FindObject("hreczvtx");
    auto Hreczvtx = BSTHnSparseHelper(HRZ);
    auto hreczvtx = Hreczvtx.GetTH1("hrecz", 3, {kINEL, kMBAND, -1, -1});
    Double_t nrecevent = hreczvtx->Integral(hreczvtx->GetXaxis()->FindBin(-10), hreczvtx->GetXaxis()->FindBin(10));

    // V0DauEta (MC)
    auto HETADAU = (THnSparse *) gROOT->FindObject("hV0DauEta");
    auto Hetadau = BSTHnSparseHelper(HETADAU);
    auto hetadau_k0pos = Hetadau.GetTH1("hetadau_K0short_positive", 3, {kINEL, kPositive, kK0short, -1});
    auto hetadau_k0neg = Hetadau.GetTH1("hetadau_K0short_negative", 3, {kINEL, kNegative, kK0short, -1});
    auto hetadau_lambdapos = Hetadau.GetTH1("hetadau_Lambda_positive", 3, {kINEL, kPositive, kLambda, -1});
    auto hetadau_lambdaneg = Hetadau.GetTH1("hetadau_Lambda_negative", 3, {kINEL, kNegative, kLambda, -1});
    auto hetadau_antilambdapos = Hetadau.GetTH1("hetadau_AntiLambda_positive", 3, {kINEL, kPositive, kAntilambda, -1});
    auto hetadau_antilambdaneg = Hetadau.GetTH1("hetadau_AntiLambda_negative", 3, {kINEL, kPositive, kAntilambda, -1});


// Data File open and find objects
    auto fdata = TFile::Open(DataFile, "open");
    auto Ddir = (TDirectory *) fdata->Get(Datadir);
    Ddir->cd();

    // reconstructed zvtx (Data)
    auto DHRZ = (THnSparse *) gROOT->FindObject("hreczvtx");
    auto DHreczvtx = BSTHnSparseHelper(DHRZ);
    auto Dhreczvtx = DHreczvtx.GetTH1("data_hrecz", 3, {kDATA, kMBAND, -1, -1});
    Double_t Dnrecevent = Dhreczvtx->Integral(Dhreczvtx->GetXaxis()->FindBin(-10), Dhreczvtx->GetXaxis()->FindBin(10));

    // v0DauEta (Data)
    auto DHETADAU = (THnSparse *) gROOT->FindObject("hV0DauEta");
    auto DHetadau = BSTHnSparseHelper(DHETADAU);
    auto Dhetadau_k0pos = DHetadau.GetTH1("Dhetadau_K0short_positive", 3, {kDATA, kPositive, kK0short, -1});
    auto Dhetadau_k0neg = DHetadau.GetTH1("Dhetadau_K0short_negative", 3, {kDATA, kNegative, kK0short, -1});
    auto Dhetadau_lambdapos = DHetadau.GetTH1("Dhetadau_Lambda_positive", 3, {kDATA, kPositive, kLambda, -1});
    auto Dhetadau_lambdaneg = DHetadau.GetTH1("Dhetadau_Lambda_negative", 3, {kDATA, kNegative, kLambda, -1});
    auto Dhetadau_antilambdapos = DHetadau.GetTH1("Dhetadau_AntiLambda_positive", 3, {kDATA, kPositive, kAntilambda, -1});
    auto Dhetadau_antilambdaneg = DHetadau.GetTH1("Dhetadau_AntiLambda_negative", 3, {kDATA, kPositive, kAntilambda, -1});
//--------------------------------------------------------------------------------------------------------------------//

// draw etadau_k0pos

    auto canvas1 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas1->Draw();
    TPad *etadau_k0posPad = canvas1->GetPad(1);
    etadau_k0posPad->SetTickx();
    etadau_k0posPad->SetGridx();
    etadau_k0posPad->SetGridy();
    etadau_k0posPad->SetLogy();
    etadau_k0posPad->cd();

    // etadua_K0shortMasscut_positive (MC)
    hetadau_k0pos->Scale(1. / nrecevent, "width");
    hetadau_k0pos->SetLineColor(kBlack);
    hetadau_k0pos->SetMarkerColor(kBlack);
    hetadau_k0pos->SetMarkerStyle(20);
    hetadau_k0pos->GetXaxis()->SetRangeUser(-2, 2);
    hetadau_k0pos->GetYaxis()->SetRangeUser(1e-2,1e2);
    hset(*hetadau_k0pos, "#it{#eta}", mstring, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hetadau_k0pos->Draw("e");

    // etadau_K0shortMasscut_positive (Data)
    Dhetadau_k0pos->Scale(1. / Dnrecevent, "width");
    Dhetadau_k0pos->SetLineColor(kRed);
    Dhetadau_k0pos->SetMarkerColor(kRed);
    Dhetadau_k0pos->SetMarkerStyle(24);
    Dhetadau_k0pos->Draw("samee");

    // legend of etadau_K0shortMasscut_positive histogram
    TLegend *legend_etadauK0pos = new TLegend(0.51,0.81,0.8,0.93, NULL, "brNDC");
    legend_etadauK0pos->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_etadauK0pos->AddEntry(hetadau_k0pos, "MC", "lep");
    legend_etadauK0pos->AddEntry(Dhetadau_k0pos, "Data", "lep");
    legend_etadauK0pos->SetTextSize(0.04);
    legend_etadauK0pos->SetBorderSize(0);
    legend_etadauK0pos->SetTextAlign(22);
    legend_etadauK0pos->Draw();

    // draw ratio of etadau_K0shortMasscut (MC/Data)
    TPad *etadau_k0pos_ratioPad = canvas1->GetPad(2);
    etadau_k0pos_ratioPad->SetTickx();
    etadau_k0pos_ratioPad->SetGridx();
    etadau_k0pos_ratioPad->SetGridy();
    etadau_k0pos_ratioPad->cd();

    auto etadau_k0posRatio = (TH1D *) hetadau_k0pos->Clone();
    etadau_k0posRatio->Divide(Dhetadau_k0pos);
    etadau_k0posRatio->SetLineColor(kRed);
    etadau_k0posRatio->SetMarkerColor(kRed);
    etadau_k0posRatio->SetMarkerStyle(24);
    etadau_k0posRatio->GetYaxis()->SetRangeUser(-10, 50);
    hset(*etadau_k0posRatio, "#it{#eta}", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    etadau_k0posRatio->Draw();

// draw etadau_k0neg

    auto canvas2 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas2->Draw();
    TPad *etadau_k0negPad = canvas2->GetPad(1);
    etadau_k0negPad->SetTickx();
    etadau_k0negPad->SetGridx();
    etadau_k0negPad->SetGridy();
    etadau_k0negPad->SetLogy();
    etadau_k0negPad->cd();

    // etadua_K0shortMasscut_positive (MC)
    hetadau_k0neg->Scale(1. / nrecevent, "width");
    hetadau_k0neg->SetLineColor(kBlack);
    hetadau_k0neg->SetMarkerColor(kBlack);
    hetadau_k0neg->SetMarkerStyle(20);
    hetadau_k0neg->GetXaxis()->SetRangeUser(-2, 2);
    hetadau_k0neg->GetYaxis()->SetRangeUser(1e-2,1e2);
    hset(*hetadau_k0neg, "#it{#eta}", mstring, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hetadau_k0neg->Draw("e");

    // etadau_K0shortMasscut_positive (Data)
    Dhetadau_k0neg->Scale(1. / Dnrecevent, "width");
    Dhetadau_k0neg->SetLineColor(kRed);
    Dhetadau_k0neg->SetMarkerColor(kRed);
    Dhetadau_k0neg->SetMarkerStyle(24);
    Dhetadau_k0neg->Draw("samee");

    // legend of etadau_K0shortMasscut_positive histogram
    TLegend *legend_etadauK0neg = new TLegend(0.51,0.81,0.8,0.93, NULL, "brNDC");
    legend_etadauK0neg->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_etadauK0neg->AddEntry(hetadau_k0neg, "MC", "lep");
    legend_etadauK0neg->AddEntry(Dhetadau_k0neg, "Data", "lep");
    legend_etadauK0neg->SetTextSize(0.04);
    legend_etadauK0neg->SetBorderSize(0);
    legend_etadauK0neg->SetTextAlign(22);
    legend_etadauK0neg->Draw();

    // draw ratio of etadau_K0shortMasscut (MC/Data)
    TPad *etadau_k0neg_ratioPad = canvas2->GetPad(2);
    etadau_k0neg_ratioPad->SetTickx();
    etadau_k0neg_ratioPad->SetGridx();
    etadau_k0neg_ratioPad->SetGridy();
    etadau_k0neg_ratioPad->cd();

    auto etadau_k0negRatio = (TH1D *) hetadau_k0neg->Clone();
    etadau_k0negRatio->Divide(Dhetadau_k0neg);
    etadau_k0negRatio->SetLineColor(kRed);
    etadau_k0negRatio->SetMarkerColor(kRed);
    etadau_k0negRatio->SetMarkerStyle(24);
    etadau_k0negRatio->GetYaxis()->SetRangeUser(-10, 50);
    hset(*etadau_k0negRatio, "#it{#eta}", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    etadau_k0negRatio->Draw();

// draw etadau_Lambdapos

    auto canvas3 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas3->Draw();
    TPad *etadau_lambdaposPad = canvas3->GetPad(1);
    etadau_lambdaposPad->SetTickx();
    etadau_lambdaposPad->SetGridx();
    etadau_lambdaposPad->SetGridy();
    etadau_lambdaposPad->SetLogy();
    etadau_lambdaposPad->cd();

    // etadua_LambdaMasscut_positive (MC)
    hetadau_lambdapos->Scale(1. / nrecevent, "width");
    hetadau_lambdapos->SetLineColor(kBlack);
    hetadau_lambdapos->SetMarkerColor(kBlack);
    hetadau_lambdapos->SetMarkerStyle(20);
    hetadau_lambdapos->GetXaxis()->SetRangeUser(-2, 2);
    hetadau_lambdapos->GetYaxis()->SetRangeUser(1e-2,1e2);
    hset(*hetadau_lambdapos, "#it{#eta}", mstring, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hetadau_lambdapos->Draw("e");

    // etadau_LambdaMasscut_positive (Data)
    Dhetadau_lambdapos->Scale(1. / Dnrecevent, "width");
    Dhetadau_lambdapos->SetLineColor(kRed);
    Dhetadau_lambdapos->SetMarkerColor(kRed);
    Dhetadau_lambdapos->SetMarkerStyle(24);
    Dhetadau_lambdapos->Draw("samee");

    // legend of etadau_LambdaMasscut_positive histogram
    TLegend *legend_etadauLambdapos = new TLegend(0.51,0.81,0.8,0.93, NULL, "brNDC");
    legend_etadauLambdapos->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_etadauLambdapos->AddEntry(hetadau_lambdapos, "MC", "lep");
    legend_etadauLambdapos->AddEntry(Dhetadau_lambdapos, "Data", "lep");
    legend_etadauLambdapos->SetTextSize(0.04);
    legend_etadauLambdapos->SetBorderSize(0);
    legend_etadauLambdapos->SetTextAlign(22);
    legend_etadauLambdapos->Draw();

    // draw ratio of etadau_lambdashortMasscut (MC/Data)
    TPad *etadau_lambdapos_ratioPad = canvas3->GetPad(2);
    etadau_lambdapos_ratioPad->SetTickx();
    etadau_lambdapos_ratioPad->SetGridx();
    etadau_lambdapos_ratioPad->SetGridy();
    etadau_lambdapos_ratioPad->cd();

    auto etadau_lambdaposRatio = (TH1D *) hetadau_lambdapos->Clone();
    etadau_lambdaposRatio->Divide(Dhetadau_lambdapos);
    etadau_lambdaposRatio->SetLineColor(kRed);
    etadau_lambdaposRatio->SetMarkerColor(kRed);
    etadau_lambdaposRatio->SetMarkerStyle(24);
    etadau_lambdaposRatio->GetYaxis()->SetRangeUser(-10, 50);
    hset(*etadau_lambdaposRatio, "#it{#eta}", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    etadau_lambdaposRatio->Draw();

// draw etadau_Lambdaneg

    auto canvas4 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas4->Draw();
    TPad *etadau_lambdanegPad = canvas4->GetPad(1);
    etadau_lambdanegPad->SetTickx();
    etadau_lambdanegPad->SetGridx();
    etadau_lambdanegPad->SetGridy();
    etadau_lambdanegPad->SetLogy();
    etadau_lambdanegPad->cd();

    // etadua_LambdaMasscut_negative (MC)
    hetadau_lambdaneg->Scale(1. / nrecevent, "width");
    hetadau_lambdaneg->SetLineColor(kBlack);
    hetadau_lambdaneg->SetMarkerColor(kBlack);
    hetadau_lambdaneg->SetMarkerStyle(20);
    hetadau_lambdaneg->GetXaxis()->SetRangeUser(-2, 2);
    hetadau_lambdaneg->GetYaxis()->SetRangeUser(1e-2,1e2);
    hset(*hetadau_lambdaneg, "#it{#eta}", mstring, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hetadau_lambdaneg->Draw("e");

    // etadau_LambdaMasscut_positive (Data)
    Dhetadau_lambdaneg->Scale(1. / Dnrecevent, "width");
    Dhetadau_lambdaneg->SetLineColor(kRed);
    Dhetadau_lambdaneg->SetMarkerColor(kRed);
    Dhetadau_lambdaneg->SetMarkerStyle(24);
    Dhetadau_lambdaneg->Draw("samee");

    // legend of etadau_LambdaMasscut_positive histogram
    TLegend *legend_etadauLambdaneg = new TLegend(0.51,0.81,0.8,0.93, NULL, "brNDC");
    legend_etadauLambdaneg->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_etadauLambdaneg->AddEntry(hetadau_lambdaneg, "MC", "lep");
    legend_etadauLambdaneg->AddEntry(Dhetadau_lambdaneg, "Data", "lep");
    legend_etadauLambdaneg->SetTextSize(0.04);
    legend_etadauLambdaneg->SetBorderSize(0);
    legend_etadauLambdaneg->SetTextAlign(22);
    legend_etadauLambdaneg->Draw();

    // draw ratio of etadau_K0shortMasscut (MC/Data)
    TPad *etadau_lambdaneg_ratioPad = canvas4->GetPad(2);
    etadau_lambdaneg_ratioPad->SetTickx();
    etadau_lambdaneg_ratioPad->SetGridx();
    etadau_lambdaneg_ratioPad->SetGridy();
    etadau_lambdaneg_ratioPad->cd();

    auto etadau_lambdanegRatio = (TH1D *) hetadau_lambdaneg->Clone();
    etadau_lambdanegRatio->Divide(Dhetadau_lambdaneg);
    etadau_lambdanegRatio->SetLineColor(kRed);
    etadau_lambdanegRatio->SetMarkerColor(kRed);
    etadau_lambdanegRatio->SetMarkerStyle(24);
    etadau_lambdanegRatio->GetYaxis()->SetRangeUser(-10, 50);
    hset(*etadau_lambdanegRatio, "#it{#eta}", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    etadau_lambdanegRatio->Draw();

// draw etadau_AntiLambdapos

    auto canvas5 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas5->Draw();
    TPad *etadau_antilambdaposPad = canvas5->GetPad(1);
    etadau_antilambdaposPad->SetTickx();
    etadau_antilambdaposPad->SetGridx();
    etadau_antilambdaposPad->SetGridy();
    etadau_antilambdaposPad->SetLogy();
    etadau_antilambdaposPad->cd();

    // etadua_AntiLambdaMasscut_positive (MC)
    hetadau_antilambdapos->Scale(1. / nrecevent, "width");
    hetadau_antilambdapos->SetLineColor(kBlack);
    hetadau_antilambdapos->SetMarkerColor(kBlack);
    hetadau_antilambdapos->SetMarkerStyle(20);
    hetadau_antilambdapos->GetXaxis()->SetRangeUser(-2, 2);
    hetadau_antilambdapos->GetYaxis()->SetRangeUser(1e-2,1e2);
    hset(*hetadau_antilambdapos, "#it{#eta}", mstring, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hetadau_antilambdapos->Draw("e");

    // etadau_LambdaMasscut_positive (Data)
    Dhetadau_antilambdapos->Scale(1. / Dnrecevent, "width");
    Dhetadau_antilambdapos->SetLineColor(kRed);
    Dhetadau_antilambdapos->SetMarkerColor(kRed);
    Dhetadau_antilambdapos->SetMarkerStyle(24);
    Dhetadau_antilambdapos->Draw("samee");

    // legend of etadau_LambdaMasscut_positive histogram
    TLegend *legend_etadauAntiLambdapos = new TLegend(0.51,0.81,0.8,0.93, NULL, "brNDC");
    legend_etadauAntiLambdapos->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_etadauAntiLambdapos->AddEntry(hetadau_antilambdapos, "MC", "lep");
    legend_etadauAntiLambdapos->AddEntry(Dhetadau_antilambdapos, "Data", "lep");
    legend_etadauAntiLambdapos->SetTextSize(0.04);
    legend_etadauAntiLambdapos->SetBorderSize(0);
    legend_etadauAntiLambdapos->SetTextAlign(22);
    legend_etadauAntiLambdapos->Draw();

    // draw ratio of etadau_lambdashortMasscut (MC/Data)
    TPad *etadau_antilambdapos_ratioPad = canvas5->GetPad(2);
    etadau_antilambdapos_ratioPad->SetTickx();
    etadau_antilambdapos_ratioPad->SetGridx();
    etadau_antilambdapos_ratioPad->SetGridy();
    etadau_antilambdapos_ratioPad->cd();

    auto etadau_antilambdaposRatio = (TH1D *) hetadau_antilambdapos->Clone();
    etadau_antilambdaposRatio->Divide(Dhetadau_antilambdapos);
    etadau_antilambdaposRatio->SetLineColor(kRed);
    etadau_antilambdaposRatio->SetMarkerColor(kRed);
    etadau_antilambdaposRatio->SetMarkerStyle(24);
    etadau_antilambdaposRatio->GetYaxis()->SetRangeUser(-10, 50);
    hset(*etadau_antilambdaposRatio, "#it{#eta}", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    etadau_antilambdaposRatio->Draw();

// draw etadau_Lambdaneg

    auto canvas6 = new Filipad2(++nn, 2, 0.4, 100, 40, 0.7, 1, 1);
    canvas6->Draw();
    TPad *etadau_antilambdanegPad = canvas6->GetPad(1);
    etadau_antilambdanegPad->SetTickx();
    etadau_antilambdanegPad->SetGridx();
    etadau_antilambdanegPad->SetGridy();
    etadau_antilambdanegPad->SetLogy();
    etadau_antilambdanegPad->cd();

    // etadua_LambdaMasscut_negative (MC)
    hetadau_antilambdaneg->Scale(1. / nrecevent, "width");
    hetadau_antilambdaneg->SetLineColor(kBlack);
    hetadau_antilambdaneg->SetMarkerColor(kBlack);
    hetadau_antilambdaneg->SetMarkerStyle(20);
    hetadau_antilambdaneg->GetXaxis()->SetRangeUser(-2, 2);
    hetadau_antilambdaneg->GetYaxis()->SetRangeUser(1e-2,1e2);
    hset(*hetadau_antilambdaneg, "#it{#eta}", mstring, 0.9, 1.3, 0.05, 0.05, 0.01, 0.01, 0.05, 0.05, 510, 510);
    hetadau_antilambdaneg->Draw("e");

    // etadau_LambdaMasscut_positive (Data)
    Dhetadau_antilambdaneg->Scale(1. / Dnrecevent, "width");
    Dhetadau_antilambdaneg->SetLineColor(kRed);
    Dhetadau_antilambdaneg->SetMarkerColor(kRed);
    Dhetadau_antilambdaneg->SetMarkerStyle(24);
    Dhetadau_antilambdaneg->Draw("samee");

    // legend of etadau_LambdaMasscut_positive histogram
    TLegend *legend_etadauAntiLambdaneg = new TLegend(0.51,0.81,0.8,0.93, NULL, "brNDC");
    legend_etadauAntiLambdaneg->SetHeader(Form("MC: %s \n Data: %s", MCDataset, DataDataset), "C");
    legend_etadauAntiLambdaneg->AddEntry(hetadau_antilambdaneg, "MC", "lep");
    legend_etadauAntiLambdaneg->AddEntry(Dhetadau_antilambdaneg, "Data", "lep");
    legend_etadauAntiLambdaneg->SetTextSize(0.04);
    legend_etadauAntiLambdaneg->SetBorderSize(0);
    legend_etadauAntiLambdaneg->SetTextAlign(22);
    legend_etadauAntiLambdaneg->Draw();

    // draw ratio of etadau_K0shortMasscut (MC/Data)
    TPad *etadau_antilambdaneg_ratioPad = canvas6->GetPad(2);
    etadau_antilambdaneg_ratioPad->SetTickx();
    etadau_antilambdaneg_ratioPad->SetGridx();
    etadau_antilambdaneg_ratioPad->SetGridy();
    etadau_antilambdaneg_ratioPad->cd();

    auto etadau_antilambdanegRatio = (TH1D *) hetadau_antilambdaneg->Clone();
    etadau_antilambdanegRatio->Divide(Dhetadau_antilambdaneg);
    etadau_antilambdanegRatio->SetLineColor(kRed);
    etadau_antilambdanegRatio->SetMarkerColor(kRed);
    etadau_antilambdanegRatio->SetMarkerStyle(24);
    etadau_antilambdanegRatio->GetYaxis()->SetRangeUser(-10, 50);
    hset(*etadau_antilambdanegRatio, "#it{#eta}", "Ratio (MC/Data)", 1.0, 0.9, 0.08, 0.07, 0.01, 0.01, 0.07, 0.07, 510, 510);
    etadau_antilambdanegRatio->Draw();


    
}