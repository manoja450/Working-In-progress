#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TLatex.h"
#include <map>
#include <TSystem.h>
#include <sys/stat.h>
#include <TGaxis.h>
#include <TStyle.h>

using namespace std;

void processLowLightEvents(const char *fileName) {
    // Create output directory
    const char* outDir = "area_plots";
    
    // Check if directory exists, if not create it
    struct stat info;
    if (stat(outDir, &info) != 0) {
        // Directory doesn't exist, try to create it
        if (gSystem->mkdir(outDir, kTRUE) != 0) {
            cerr << "Error: Could not create directory " << outDir << endl;
            return;
        }
        cout << "Created directory: " << outDir << endl;
    } 
    else if (!(info.st_mode & S_IFDIR)) {
        cerr << "Error: " << outDir << " exists but is not a directory" << endl;
        return;
    }

    TFile *file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        cerr << "Error opening file: " << fileName << endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("tree");
    if (!tree) {
        cerr << "Error accessing TTree 'tree'!" << endl;
        file->Close();
        return;
    }

    Short_t adcVal[23][45];
    Double_t area[23];
    Int_t triggerBits;

    tree->SetBranchAddress("adcVal", adcVal);
    tree->SetBranchAddress("area", area);
    tree->SetBranchAddress("triggerBits", &triggerBits);

    Long64_t nEntries = tree->GetEntries();

    // PMT and SiPM channel mappings
    int pmtChannelMap[12] = {0, 10, 7, 2, 6, 3, 8, 9, 11, 4, 5, 1};
    int sipmChannelMap[10] = {12, 13, 14, 15, 16, 17, 18, 19, 20, 21};

    // Create histograms for PMTs and SiPMs with area
    TH1F *histPMT[12];
    TH1F *histSiPM[10];
    for (int i = 0; i < 12; i++) {
        histPMT[i] = new TH1F(Form("PMT%d", i + 1), 
                              Form("PMT%d;Area;Events/550 ADC", i + 1),
                              100, -5000, 50000);
        histPMT[i]->SetLineColor(kRed);
    }
    for (int i = 0; i < 10; i++) {
        histSiPM[i] = new TH1F(Form("SiPM%d", i + 1), 
                               Form("SiPM%d;Area;Events/55 ADC", i + 1),
                               100, -500, 5000);
        histSiPM[i]->SetLineColor(kBlue);
    }

    // Process each event
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        tree->GetEntry(entry);
        if (triggerBits == 34) {
            // Fill PMT histograms
            for (int pmt = 0; pmt < 12; pmt++) {
                int adcIdx = pmtChannelMap[pmt];
                histPMT[pmt]->Fill(area[adcIdx]);
            }
            // Fill SiPM histograms
            for (int sipm = 0; sipm < 10; sipm++) {
                int adcIdx = sipmChannelMap[sipm];
                histSiPM[sipm]->Fill(area[adcIdx]);
            }
        }
    }

    // Save individual PMT plots
    TCanvas *individualCanvas = new TCanvas("IndividualCanvas", "Individual Plots", 800, 600);
    for (int i = 0; i < 12; i++) {
        individualCanvas->Clear();
        
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.12);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.12);

        histPMT[i]->Draw();
        
        histPMT[i]->GetXaxis()->SetTitleSize(0.06);
        histPMT[i]->GetYaxis()->SetTitleSize(0.06);
        histPMT[i]->GetYaxis()->SetTitleOffset(1.2);

        individualCanvas->SaveAs(Form("%s/PMT%d_area.png", outDir, i + 1));
    }

    // Save individual SiPM plots
    for (int i = 0; i < 10; i++) {
        individualCanvas->Clear();
        
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.12);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.12);

        histSiPM[i]->Draw();
        
        histSiPM[i]->GetXaxis()->SetTitleSize(0.06);
        histSiPM[i]->GetYaxis()->SetTitleSize(0.06);
        histSiPM[i]->GetYaxis()->SetTitleOffset(1.2);

        individualCanvas->SaveAs(Form("%s/SiPM%d_area.png", outDir, i + 1));
    }

    // Create master canvas with 6x5 layout
    TCanvas *masterCanvas = new TCanvas("MasterCanvas", "Combined PMT and SiPM Area Distributions", 3600, 3000);
    masterCanvas->Divide(5, 6, 0.005, 0.005);

    // Define the custom layout
    int layout[6][5] = {
        {-1,  -1,  20,  21, -1},
        {16,  9,   3,   7,  12},
        {15,  5,   4,   8,   -1},
        {19,  0,   6,   1,  17},
        {-1,  10,  11,  2,  13},
        {-1,  -1,  14,  18, -1}
    };

    // Reverse map for PMT channels (channel -> pmtIndex)
    map<int, int> pmtReverseMap;
    for (int i = 0; i < 12; i++) {
        pmtReverseMap[pmtChannelMap[i]] = i;
    }

    // Configure scientific notation
    TGaxis::SetMaxDigits(3);
    gStyle->SetPaintTextFormat("4.1e");

    // Populate each pad in the master canvas
    for (int row = 0; row < 6; row++) {
        for (int col = 0; col < 5; col++) {
            int padNum = row * 5 + col + 1;
            masterCanvas->cd(padNum);
            
            int channel = layout[row][col];
            if (channel == -1) {
                gPad->SetFillColor(0);
                gPad->Modified();
                gPad->Update();
                continue;
            }

            TH1F *hist = nullptr;
            TString title;

            if (channel >= 0 && channel < 12) { // PMT
                auto it = pmtReverseMap.find(channel);
                if (it == pmtReverseMap.end()) {
                    cerr << "PMT channel " << channel << " not mapped!" << endl;
                    continue;
                }
                int pmtIdx = it->second;
                hist = histPMT[pmtIdx];
                title = Form("PMT%d", pmtIdx + 1);
            } else if (channel >= 12 && channel <= 21) { // SiPM
                int sipmIdx = channel - 12;
                if (sipmIdx < 0 || sipmIdx >= 10) {
                    cerr << "SiPM channel " << channel << " invalid!" << endl;
                    continue;
                }
                hist = histSiPM[sipmIdx];
                title = Form("SiPM%d", sipmIdx + 1);
            } else {
                cerr << "Invalid channel: " << channel << endl;
                continue;
            }

            // Configure histogram
            hist->SetTitle("");
            hist->GetXaxis()->SetTitleSize(0.05);
            hist->GetYaxis()->SetTitleSize(0.07);
            hist->GetXaxis()->SetLabelSize(0.05);
            hist->GetYaxis()->SetLabelSize(0.05);
            hist->GetYaxis()->SetTitleOffset(1.2);
            hist->GetXaxis()->SetTitleOffset(1.1);
            hist->GetXaxis()->SetNdivisions(505);
            hist->GetYaxis()->SetNdivisions(505);

            // Adjust pad margins
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.10);
            gPad->SetBottomMargin(0.12);
            gPad->SetTopMargin(0.12);

            hist->Draw();

            // Force scientific notation
            gPad->Update();
            hist->GetXaxis()->SetMoreLogLabels();
            hist->GetYaxis()->SetMoreLogLabels();

            // Add title
            TLatex *latex = new TLatex();
            latex->SetTextSize(0.12);
            latex->SetTextAlign(22);
            latex->SetNDC(true);
            latex->DrawLatex(0.5, 0.93, title);
        }
    }

    // Add axis labels textbox
    masterCanvas->cd();
    TLatex *textbox = new TLatex();
    textbox->SetTextSize(0.05);
    textbox->SetTextAlign(15);
    textbox->SetNDC(true);
    textbox->DrawLatex(0.01, 0.11, "X axis: Area");
    textbox->DrawLatex(0.01, 0.07, "Y axis: Events");

    // Save and clean up
    masterCanvas->SaveAs(Form("%s/Combined_PMT_SiPM_Area_Distributions.png", outDir));
    for (int i = 0; i < 12; i++) delete histPMT[i];
    for (int i = 0; i < 10; i++) delete histSiPM[i];
    delete individualCanvas;
    delete masterCanvas;
    file->Close();

    cout << "All plots saved in directory: " << outDir << endl;
    cout << "Individual PMT and SiPM area plots saved as PMT1_area.png, PMT2_area.png, etc." << endl;
    cout << "Combined area histogram saved as Combined_PMT_SiPM_Area_Distributions.png" << endl;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <root_file>" << endl;
        return 1;
    }
    processLowLightEvents(argv[1]);
    return 0;
}