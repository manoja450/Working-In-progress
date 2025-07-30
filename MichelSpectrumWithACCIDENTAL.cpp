#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TGaxis.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <string>
#include <map>
#include <set>
#include <sys/stat.h>
#include <unistd.h>
#include <ctime>
#include <TGraph.h>
#include <TPad.h>

using std::cout;
using std::endl;
using namespace std;

// Constants
const int N_PMTS = 12;
const int PMT_CHANNEL_MAP[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const int PULSE_THRESHOLD = 30;     // ADC threshold for pulse detection
const int BS_UNCERTAINTY = 5;       // Baseline uncertainty (ADC)
const int EV61_THRESHOLD = 1200;    // Beam on if channel 22 > this (ADC)
const double MUON_ENERGY_THRESHOLD = 50; // Min PMT energy for muon (p.e.)
const double MICHEL_ENERGY_MIN = 40;    // Min PMT energy for Michel (p.e.)
const double MICHEL_ENERGY_MAX = 1000;  // Max PMT energy for Michel (p.e.)
const double MICHEL_ENERGY_MAX_DT = 400; // Max PMT energy for dt plots (p.e.)
const double MICHEL_DT_MIN = 0.76;       // Min time after muon for Michel (µs)
const double MICHEL_DT_MAX = 16.0;      // Max time after muon for Michel (µs)
const double MIN_MICHEL_SEPARATION = 0.5; // Minimum time between Michel candidates (µs)
const int ADCSIZE = 45;                 // Number of ADC samples per waveform

// Generate unique output directory with timestamp
string getTimestamp() {
    time_t now = time(nullptr);
    struct tm *t = localtime(&now);
    char buffer[20];
    strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", t);
    return string(buffer);
}
const string OUTPUT_DIR = "./AnalysisOutput_" + getTimestamp();

const std::vector<double> SIDE_VP_THRESHOLDS = {750, 950, 1200, 1375, 525, 700, 700, 500}; // Channels 12-19 (ADC)
const double TOP_VP_THRESHOLD = 450; // Channels 20-21 (ADC)
const double FIT_MIN = 1.0; // Fit range min (µs)
const double FIT_MAX = 16.0; // Fit range max (µs)

// Pulse structure
struct pulse {
    double start;          // Start time (µs)
    double end;            // End time (µs)
    double peak;           // Max amplitude (p.e. for PMTs, ADC for SiPMs)
    double energy;         // Energy (p.e. for PMTs, ADC for SiPMs)
    double number;         // Number of channels with pulse
    bool single;           // Timing consistency
    bool beam;             // Beam status
    double trigger;        // Trigger type
    double side_vp_energy; // Side veto energy (ADC)
    double top_vp_energy;  // Top veto energy (ADC)
    double all_vp_energy;  // All veto energy (ADC)
    double last_muon_time; // Time of last muon (µs)
    bool is_muon;          // Muon candidate flag
    bool is_michel;        // Michel electron candidate flag
    bool veto_hit[10];     // Which veto panels were hit (channels 12-21)
};

// Temporary pulse structure
struct pulse_temp {
    double start;  // Start time (µs)
    double end;    // End time (µs)
    double peak;   // Max amplitude
    double energy; // Energy
};

// SPE fitting function
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0] - sqrt(2) * par[4]) / sqrt(2 * pow(par[5], 2) - pow(par[2], 2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0] - sqrt(3) * par[4]) / sqrt(3 * pow(par[5], 2) - 2 * pow(par[2], 2)), 2));
    return term1 + term2 + term3 + term4;
}

// Exponential fit function: N0 * exp(-t/tau) + C (t, tau in µs)
Double_t ExpFit(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0] / par[1]) + par[2];
}

// Exponential fit function without constant term (for accidental)
Double_t ExpFitNoC(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0] / par[1]);
}

// Double exponential fit function: N1*exp(-t/tau1) + N2*exp(-t/tau2)
Double_t DoubleExpFit(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0] / par[1]) + par[2] * exp(-x[0] / par[3]);
}

// Utility functions
template<typename T>
double getAverage(const std::vector<T>& v) {
    if (v.empty()) return 0;
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

template<typename T>
double mostFrequent(const std::vector<T>& v) {
    if (v.empty()) return 0;
    std::map<T, int> count;
    for (const auto& val : v) count[val]++;
    T most_common = v[0];
    int max_count = 0;
    for (const auto& pair : count) {
        if (pair.second > max_count) {
            max_count = pair.second;
            most_common = pair.first;
        }
    }
    return max_count > 1 ? most_common : getAverage(v);
}

template<typename T>
double variance(const std::vector<T>& v) {
    if (v.size() <= 1) return 0;
    double mean = getAverage(v);
    double sum = 0;
    for (const auto& val : v) {
        sum += (val - mean) * (val - mean);
    }
    return sum / (v.size() - 1);
}

// Create output directory
void createOutputDirectory(const string& dirName) {
    struct stat st;
    if (stat(dirName.c_str(), &st) != 0) {
        if (mkdir(dirName.c_str(), 0755) != 0) {
            cerr << "Error: Could not create directory " << dirName << endl;
            exit(1);
        }
        cout << "Created output directory: " << dirName << endl;
    } else {
        cout << "Output directory already exists: " << dirName << endl;
    }
}

// SPE calibration function
void performCalibration(const string &calibFileName, Double_t *mu1, Double_t *mu1_err) {
    TFile *calibFile = TFile::Open(calibFileName.c_str());
    if (!calibFile || calibFile->IsZombie()) {
        cerr << "Error opening calibration file: " << calibFileName << endl;
        exit(1);
    }

    TTree *calibTree = (TTree*)calibFile->Get("tree");
    if (!calibTree) {
        cerr << "Error accessing tree in calibration file" << endl;
        calibFile->Close();
        exit(1);
    }

    TCanvas *c = new TCanvas("c", "SPE Fits", 1200, 800);
    TH1F *histArea[N_PMTS];
    Long64_t nLEDFlashes[N_PMTS] = {0};
    for (int i = 0; i < N_PMTS; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area", i + 1),
                               Form("PMT %d;ADC Counts;Events", i + 1), 150, -50, 400);
    }

    Int_t triggerBits;
    Double_t area[23];
    calibTree->SetBranchAddress("triggerBits", &triggerBits);
    calibTree->SetBranchAddress("area", area);

    Long64_t nEntries = calibTree->GetEntries();
    cout << "Processing " << nEntries << " calibration events from " << calibFileName << "..." << endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        calibTree->GetEntry(entry);
        if (triggerBits != 16) continue;
        for (int pmt = 0; pmt < N_PMTS; pmt++) {
            histArea[pmt]->Fill(area[PMT_CHANNEL_MAP[pmt]]);
            nLEDFlashes[pmt]++;
        }
    }

    for (int i = 0; i < N_PMTS; i++) {
        TF1 *fitFunc = new TF1("fitFunc", SPEfit, -50, 400, 8);
        Double_t histMean = histArea[i]->GetMean();
        Double_t histRMS = histArea[i]->GetRMS();

        fitFunc->SetParameters(1000, histMean - histRMS, histRMS / 2,
                              1000, histMean, histRMS,
                              500, 200);

        histArea[i]->Fit(fitFunc, "Q", "", -50, 400);

        mu1[i] = fitFunc->GetParameter(4);
        Double_t sigma_mu1 = fitFunc->GetParError(4);
        Double_t sigma1 = fitFunc->GetParameter(5);
        mu1_err[i] = sqrt(pow(sigma_mu1, 2) + pow(sigma1 / sqrt(nLEDFlashes[i]), 2));

        // Plot SPE fit
        c->Clear();
        histArea[i]->Draw();
        fitFunc->Draw("same");
        TLegend *leg = new TLegend(0.6, 0.7, 0.9, 0.9);
        leg->AddEntry(histArea[i], Form("PMT %d Data", i + 1), "l");
        leg->AddEntry(fitFunc, "SPE Fit", "l");
        leg->AddEntry((TObject*)0, Form("mu1 = %.2f #pm %.2f", mu1[i], mu1_err[i]), "");
        leg->Draw();
        string plotName = OUTPUT_DIR + Form("/SPE_Fit_PMT%d.png", i + 1);
        c->Update();
        c->SaveAs(plotName.c_str());
        cout << "Saved SPE plot: " << plotName << endl;
        delete leg;
        delete fitFunc;
        delete histArea[i];
    }

    delete c;
    calibFile->Close();
}

void createVetoPanelPlots(TH1D* h_veto_panel[10], const string& outputDir) {
    for (int i = 0; i < 10; i++) {
        TCanvas *c = new TCanvas(Form("c_veto_%d", i+12), Form("Veto Panel %d", i+12), 1200, 800);
        gStyle->SetOptStat(1111);
        gStyle->SetOptTitle(1);
        gStyle->SetStatX(0.9);
        gStyle->SetStatY(0.9);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.15);
        
        h_veto_panel[i]->SetLineColor(kBlack);
        h_veto_panel[i]->SetLineWidth(2);
        h_veto_panel[i]->Draw("hist");
        
        string plotName = outputDir + Form("/Veto_Panel_%d.png", i+12);
        c->SaveAs(plotName.c_str());
        cout << "Saved veto panel plot: " << plotName << endl;
        delete c;
    }

    TCanvas *c_combined = new TCanvas("c_veto_combined", "Combined Veto Panel Energies", 1600, 1200);
    c_combined->Divide(4, 3);

    for (int i = 0; i < 10; i++) {
        c_combined->cd(i+1);
        h_veto_panel[i]->SetLineColor(kBlack);
        h_veto_panel[i]->SetLineWidth(2);
        h_veto_panel[i]->SetTitle("");
        h_veto_panel[i]->Draw("hist");
    }

    string combinedPlotName = outputDir + "/Combined_Veto_Panels.png";
    c_combined->SaveAs(combinedPlotName.c_str());
    cout << "Saved combined veto panel plot: " << combinedPlotName << endl;
    delete c_combined;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <calibration_file> <input_file1> [<input_file2> ...]" << endl;
        return -1;
    }

    string calibFileName = argv[1];
    vector<string> inputFiles;
    for (int i = 2; i < argc; i++) {
        inputFiles.push_back(argv[i]);
    }

    createOutputDirectory(OUTPUT_DIR);

    cout << "Calibration file: " << calibFileName << endl;
    cout << "Input files:" << endl;
    for (const auto& file : inputFiles) {
        cout << "  " << file << endl;
    }

    if (gSystem->AccessPathName(calibFileName.c_str())) {
        cerr << "Error: Calibration file " << calibFileName << " not found" << endl;
        return -1;
    }

    bool anyInputFileExists = false;
    for (const auto& file : inputFiles) {
        if (!gSystem->AccessPathName(file.c_str())) {
            anyInputFileExists = true;
            break;
        }
    }
    if (!anyInputFileExists) {
        cerr << "Error: No input files found" << endl;
        return -1;
    }

    Double_t mu1[N_PMTS] = {0};
    Double_t mu1_err[N_PMTS] = {0};
    performCalibration(calibFileName, mu1, mu1_err);

    cout << "SPE Calibration Results (from " << calibFileName << "):\n";
    for (int i = 0; i < N_PMTS; i++) {
        cout << "PMT " << i + 1 << ": mu1 = " << mu1[i] << " ± " << mu1_err[i] << " ADC counts/p.e.\n";
    }

    int num_muons = 0;
    int num_michels = 0;
    int num_events = 0;
    int num_accidental = 0;

    std::map<int, int> trigger_counts;

    TH1D* h_muon_energy = new TH1D("muon_energy", "Muon Energy Distribution (with Michel Electrons);Energy (p.e.);Counts/100 p.e.", 550, -500, 5000);
    TH1D* h_michel_energy = new TH1D("michel_energy", "Michel Electron Energy Distribution;Energy (p.e.);Counts/8 p.e.", 100, 0, 800);
    TH1D* h_dt_michel = new TH1D("DeltaT", "Muon-Michel Time Difference ;Time to Previous event(Muon)(#mus);Counts/0.08 #mus", 200, 0, MICHEL_DT_MAX);
    TH1D* h_dt_michel_double = new TH1D("DeltaT_double", "Muon-Michel Time Difference (Double Exp Fit);Time to Previous event(Muon)(#mus);Counts/0.08 #mus", 200, 0, MICHEL_DT_MAX);
    TH2D* h_energy_vs_dt = new TH2D("energy_vs_dt", "Michel Energy vs Time Difference;dt (#mus);Energy (p.e.)", 160, 0, 16, 200, 0, 1000);
    TH1D* h_side_vp_muon = new TH1D("side_vp_muon", "Side Veto Energy for Muons;Energy (ADC);Counts", 200, 0, 5000);
    TH1D* h_top_vp_muon = new TH1D("top_vp_muon", "Top Veto Energy for Muons;Energy (ADC);Counts", 200, 0, 1000);
    TH1D* h_trigger_bits = new TH1D("trigger_bits", "Trigger Bits Distribution;Trigger Bits;Counts", 36, 0, 36);
    
    TH1D* h_veto_panel[10];
    const char* veto_names[10] = {
        "Veto Panel 12 ", "Veto Panel 13", "Veto Panel 14", "Veto Panel 15",
        "Veto Panel 16 ", "Veto Panel 17", "Veto Panel 18", "Veto Panel 19",
        "Veto Panel 20 ", "Veto Panel 21"
    };
    
    for (int i = 0; i < 10; i++) {
        h_veto_panel[i] = new TH1D(Form("h_veto_panel_%d", i+12), 
            Form("%s;Energy (ADC);Counts", veto_names[i]), 
            200, 0, 6000);
    }

    TH1D* h_dt_accidental = new TH1D("dt_accidental", "Accidental Events Time Difference;Time to Previous event(Muon)(#mus);Counts/0.08 #mus", 
                                     200, 0, MICHEL_DT_MAX);
    TH1D* h_dt_michel_true = new TH1D("dt_michel_true", "True Michel Events Time Difference;Time to Previous event(Muon)(#mus);Counts/0.08 #mus", 
                                       200, 0, MICHEL_DT_MAX);
    TH1D* h_accidental_energy = new TH1D("accidental_energy", "Accidental Events Energy Distribution;Energy (p.e.);Counts/8 p.e.", 
                                         100, 0, 800);

    for (const auto& inputFileName : inputFiles) {
        if (gSystem->AccessPathName(inputFileName.c_str())) {
            cout << "Could not open file: " << inputFileName << ". Skipping..." << endl;
            continue;
        }

        TFile *f = new TFile(inputFileName.c_str());
        cout << "Processing file: " << inputFileName << endl;

        TTree* t = (TTree*)f->Get("tree");
        if (!t) {
            cout << "Could not find tree in file: " << inputFileName << endl;
            f->Close();
            continue;
        }

        Int_t eventID;
        Int_t nSamples[23];
        Short_t adcVal[23][45];
        Double_t baselineMean[23];
        Double_t baselineRMS[23];
        Double_t pulseH[23];
        Int_t peakPosition[23];
        Double_t area[23];
        Long64_t nsTime;
        Int_t triggerBits;

        t->SetBranchAddress("eventID", &eventID);
        t->SetBranchAddress("nSamples", nSamples);
        t->SetBranchAddress("adcVal", adcVal);
        t->SetBranchAddress("baselineMean", baselineMean);
        t->SetBranchAddress("baselineRMS", baselineRMS);
        t->SetBranchAddress("pulseH", pulseH);
        t->SetBranchAddress("peakPosition", peakPosition);
        t->SetBranchAddress("area", area);
        t->SetBranchAddress("nsTime", &nsTime);
        t->SetBranchAddress("triggerBits", &triggerBits);

        int numEntries = t->GetEntries();
        cout << "Processing " << numEntries << " entries in " << inputFileName << endl;
        double last_muon_time = 0.0;
        std::set<double> michel_muon_times;
        std::vector<std::pair<double, double>> sprung_candidates;
        std::vector<double> michel_times;
        std::vector<double> accidental_times;

        for (int iEnt = 0; iEnt < numEntries; iEnt++) {
            t->GetEntry(iEnt);
            num_events++;

            h_trigger_bits->Fill(triggerBits);
            trigger_counts[triggerBits]++;
            if (triggerBits < 0 || triggerBits > 36) {
                cout << "Warning: triggerBits = " << triggerBits << " out of histogram range (0-31) in file " << inputFileName << ", event " << eventID << endl;
            }

            struct pulse p;
            p.start = nsTime / 1000.0;
            p.end = nsTime / 1000.0;
            p.peak = 0;
            p.energy = 0;
            p.number = 0;
            p.single = false;
            p.beam = false;
            p.trigger = triggerBits;
            p.side_vp_energy = 0;
            p.top_vp_energy = 0;
            p.all_vp_energy = 0;
            p.last_muon_time = last_muon_time;
            p.is_muon = false;
            p.is_michel = false;
            for (int i = 0; i < 10; i++) p.veto_hit[i] = false;

            std::vector<double> all_chan_start, all_chan_end, all_chan_peak, all_chan_energy;
            std::vector<double> side_vp_energy, top_vp_energy;
            std::vector<double> chan_starts_no_outliers;
            TH1D h_wf("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

            bool pulse_at_end = false;
            int pulse_at_end_count = 0;
            std::vector<double> veto_energies(10, 0);

            for (int iChan = 0; iChan < 23; iChan++) {
                for (int i = 0; i < ADCSIZE; i++) {
                    h_wf.SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);
                }

                if (iChan == 22) {
                    double ev61_energy = 0;
                    for (int iBin = 1; iBin <= ADCSIZE; iBin++) {
                        ev61_energy += h_wf.GetBinContent(iBin);
                    }
                    if (ev61_energy > EV61_THRESHOLD) {
                        p.beam = true;
                    }
                }

                std::vector<pulse_temp> pulses_temp;
                bool onPulse = false;
                int thresholdBin = 0, peakBin = 0;
                double peak = 0, pulseEnergy = 0;
                double allPulseEnergy = 0;

                for (int iBin = 1; iBin <= ADCSIZE; iBin++) {
                    double iBinContent = h_wf.GetBinContent(iBin);
                    if (iBin > 15) allPulseEnergy += iBinContent;

                    if (!onPulse && iBinContent >= PULSE_THRESHOLD) {
                        onPulse = true;
                        thresholdBin = iBin;
                        peakBin = iBin;
                        peak = iBinContent;
                        pulseEnergy = iBinContent;
                    } else if (onPulse) {
                        pulseEnergy += iBinContent;
                        if (peak < iBinContent) {
                            peak = iBinContent;
                            peakBin = iBin;
                        }
                        if (iBinContent < BS_UNCERTAINTY || iBin == ADCSIZE) {
                            pulse_temp pt;
                            pt.start = thresholdBin * 16.0 / 1000.0;
                            pt.peak = iChan <= 11 && mu1[iChan] > 0 ? peak / mu1[iChan] : peak;
                            pt.end = iBin * 16.0 / 1000.0;
                            for (int j = peakBin - 1; j >= 1 && h_wf.GetBinContent(j) > BS_UNCERTAINTY; j--) {
                                if (h_wf.GetBinContent(j) > peak * 0.1) {
                                    pt.start = j * 16.0 / 1000.0;
                                }
                                pulseEnergy += h_wf.GetBinContent(j);
                            }
                            if (iChan <= 11) {
                                pt.energy = mu1[iChan] > 0 ? pulseEnergy / mu1[iChan] : 0;
                                all_chan_start.push_back(pt.start);
                                all_chan_end.push_back(pt.end);
                                all_chan_peak.push_back(pt.peak);
                                all_chan_energy.push_back(pt.energy);
                                if (pt.energy > 1) p.number += 1;
                            }
                            pulses_temp.push_back(pt);
                            peak = 0;
                            peakBin = 0;
                            pulseEnergy = 0;
                            thresholdBin = 0;
                            onPulse = false;
                        }
                    }
                }

                if (iChan >= 12 && iChan <= 19) {
                    side_vp_energy.push_back(allPulseEnergy);
                    veto_energies[iChan - 12] = allPulseEnergy;
                    if (allPulseEnergy > SIDE_VP_THRESHOLDS[iChan - 12]) {
                        p.veto_hit[iChan - 12] = true;
                    }
                } else if (iChan >= 20 && iChan <= 21) {
                    double factor = (iChan == 20) ? 1.07809 : 1.0;
                    top_vp_energy.push_back(allPulseEnergy * factor);
                    veto_energies[iChan - 12] = allPulseEnergy * factor;
                    if (allPulseEnergy * factor > TOP_VP_THRESHOLD) {
                        p.veto_hit[iChan - 12] = true;
                    }
                }

                if (iChan <= 11 && h_wf.GetBinContent(ADCSIZE) > 100) {
                    pulse_at_end_count++;
                    if (pulse_at_end_count >= 10) pulse_at_end = true;
                }

                h_wf.Reset();
            }

            p.start += mostFrequent(all_chan_start);
            p.end += mostFrequent(all_chan_end);
            p.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0);
            p.peak = std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0);
            p.side_vp_energy = std::accumulate(side_vp_energy.begin(), side_vp_energy.end(), 0.0);
            p.top_vp_energy = std::accumulate(top_vp_energy.begin(), top_vp_energy.end(), 0.0);
            p.all_vp_energy = p.side_vp_energy + p.top_vp_energy;

            for (const auto& start : all_chan_start) {
                if (fabs(start - mostFrequent(all_chan_start)) < 10 * 16.0 / 1000.0) {
                    chan_starts_no_outliers.push_back(start);
                }
            }
            p.single = (variance(chan_starts_no_outliers) < 5 * 16.0 / 1000.0);

            bool veto_hit = false;
            for (size_t i = 0; i < SIDE_VP_THRESHOLDS.size(); i++) {
                if (veto_energies[i] > SIDE_VP_THRESHOLDS[i]) {
                    veto_hit = true;
                    break;
                }
            }
            if (!veto_hit && p.top_vp_energy > TOP_VP_THRESHOLD) veto_hit = true;

            if ((p.energy > MUON_ENERGY_THRESHOLD && veto_hit) ||
                (pulse_at_end && p.energy > MUON_ENERGY_THRESHOLD / 2 && veto_hit)) {
                p.is_muon = true;
                last_muon_time = p.start;
                num_muons++;
                sprung_candidates.emplace_back(p.start, p.energy);
                h_side_vp_muon->Fill(p.side_vp_energy);
                h_top_vp_muon->Fill(p.top_vp_energy);
                
                for (int i = 0; i < 10; i++) {
                    if (p.veto_hit[i]) {
                        h_veto_panel[i]->Fill(veto_energies[i]);
                    }
                }
            }

            double dt = p.start - last_muon_time;
            bool veto_low = true;
            for (size_t i = 0; i < SIDE_VP_THRESHOLDS.size(); i++) {
                if (veto_energies[i] > SIDE_VP_THRESHOLDS[i]) {
                    veto_low = false;
                    break;
                }
            }
            if (veto_energies[8] > TOP_VP_THRESHOLD || veto_energies[9] > TOP_VP_THRESHOLD) {
                veto_low = false;
            }

            bool is_michel_candidate = p.energy >= MICHEL_ENERGY_MIN &&
                                      p.energy <= MICHEL_ENERGY_MAX &&
                                      dt >= MICHEL_DT_MIN &&
                                      dt <= MICHEL_DT_MAX &&
                                      p.number >= 8 &&
                                      veto_low &&
                                      p.trigger != 1 &&
                                      p.trigger != 4 &&
                                      p.trigger != 8 &&
                                      p.trigger != 16;
            h_energy_vs_dt->Fill(dt, p.energy);

            bool is_michel_for_dt = is_michel_candidate && p.energy <= MICHEL_ENERGY_MAX_DT;

            if (is_michel_candidate) {
                p.is_michel = true;
                num_michels++;
                michel_muon_times.insert(last_muon_time);
                michel_times.push_back(dt);
                h_michel_energy->Fill(p.energy);
            }

            if (is_michel_for_dt) {
                h_dt_michel->Fill(dt);
                h_dt_michel_double->Fill(dt);
            }

            if (is_michel_candidate) {
                // Look ahead for other Michel candidates in the same window
                for (int jEnt = iEnt + 1; jEnt < numEntries; jEnt++) {
                    t->GetEntry(jEnt);
                    double next_time = nsTime / 1000.0;
                    double next_dt = next_time - last_muon_time;
                    
                    // Exit if beyond 16µs window
                    if (next_dt > MICHEL_DT_MAX) break;
                    
                    // Skip if too close in time (potential double-counting)
                    if (next_time - p.start < MIN_MICHEL_SEPARATION) continue;
                    
                    bool next_veto_low = true;
                    for (int iChan = 12; iChan <= 21; iChan++) {
                        double veto_energy = 0;
                        for (int iBin = 16; iBin <= ADCSIZE; iBin++) {
                            veto_energy += adcVal[iChan][iBin] - baselineMean[iChan];
                        }
                        if ((iChan <= 19 && veto_energy > SIDE_VP_THRESHOLDS[iChan-12]) || 
                            (iChan >= 20 && veto_energy > TOP_VP_THRESHOLD)) {
                            next_veto_low = false;
                            break;
                        }
                    }
                    
                    double next_energy = 0;
                    for (int iChan = 0; iChan < 12; iChan++) {
                        double chan_energy = 0;
                        for (int iBin = 1; iBin <= ADCSIZE; iBin++) {
                            chan_energy += adcVal[iChan][iBin] - baselineMean[iChan];
                        }
                        if (mu1[iChan] > 0) next_energy += chan_energy / mu1[iChan];
                    }
                    
                    bool next_is_michel = next_energy >= MICHEL_ENERGY_MIN &&
                                         next_energy <= MICHEL_ENERGY_MAX &&
                                         next_dt >= MICHEL_DT_MIN &&
                                         next_dt <= MICHEL_DT_MAX &&
                                         next_veto_low &&
                                         triggerBits != 1 &&
                                         triggerBits != 4 &&
                                         triggerBits != 8 &&
                                         triggerBits != 16;
                    
                    if (next_is_michel) {
                        cout << "Found accidental-true Michel pair: " 
                             << "First at dt=" << dt << "µs (accidental), "
                             << "Second at dt=" << next_dt << "µs (true)" << endl;
                             
                        h_dt_accidental->Fill(dt);
                        h_accidental_energy->Fill(p.energy);
                        accidental_times.push_back(dt);
                        num_accidental++;
                        
                        h_dt_michel_true->Fill(next_dt);
                        michel_times.push_back(next_dt);
                        break;
                    }
                }
            }

            p.last_muon_time = last_muon_time;
        }

        for (const auto& muon : sprung_candidates) {
            if (michel_muon_times.find(muon.first) != michel_muon_times.end()) {
                h_muon_energy->Fill(muon.second);
            }
        }

        cout << "File " << inputFileName << " Statistics:\n";
        cout << "Total Events: " << num_events << "\n";
        cout << "Muons Detected: " << num_muons << "\n";
        cout << "Michel Electrons Detected: " << num_michels << "\n";
        cout << "Accidental Events Detected: " << num_accidental << "\n";
        cout << "------------------------\n";

        f->Close();

        num_events = 0;
        num_muons = 0;
        num_michels = 0;
        num_accidental = 0;
    }

    cout << "Trigger Bits Distribution (all files):\n";
    for (const auto& pair : trigger_counts) {
        cout << "Trigger " << pair.first << ": " << pair.second << " events\n";
    }
    cout << "------------------------\n";

    // Diagnostic: Check the number of events in h_dt_michel_double
    double michel_double_integral = h_dt_michel_double->Integral();
    cout << "Number of events in h_dt_michel_double (MICHEL_ENERGY_MAX_DT = " << MICHEL_ENERGY_MAX_DT << "): " << michel_double_integral << endl;
    if (michel_double_integral < 100) {
        cout << "Warning: Low statistics in h_dt_michel_double (" << michel_double_integral << " events). Fit may be unstable." << endl;
    }

    TCanvas *c = new TCanvas("c", "Analysis Plots", 1200, 800);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    c->Clear();
    h_muon_energy->SetLineColor(kBlue);
    h_muon_energy->Draw();
    c->Update();
    string plotName = OUTPUT_DIR + "/Muon_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    c->Clear();
    h_michel_energy->SetLineColor(kRed);
    h_michel_energy->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Michel_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Original Michel dt plot with single exponential fit
    c->Clear();
    h_dt_michel->SetLineWidth(2);
    h_dt_michel->SetLineColor(kBlack);
    h_dt_michel->GetXaxis()->SetTitle("Time to previous event (Muon) [#mus]");
    h_dt_michel->SetTitle("Michel Electron Time Difference (Single Exp Fit)");
    h_dt_michel->Draw("HIST");

    TF1* expFit = new TF1("expFit", ExpFit, FIT_MIN, FIT_MAX, 3);
    double integral = h_dt_michel->Integral(h_dt_michel->FindBin(FIT_MIN), h_dt_michel->FindBin(FIT_MAX));
    double bin_width = h_dt_michel->GetBinWidth(1);
    double N0_init = integral * bin_width / (FIT_MAX - FIT_MIN);
    double C_init = 0;
    
    int bin_12 = h_dt_michel->FindBin(12.0);
    int bin_16 = h_dt_michel->FindBin(16.0);
    double min_content = 1e9;
    for (int i = bin_12; i <= bin_16; i++) {
        double content = h_dt_michel->GetBinContent(i);
        if (content > 0 && content < min_content) min_content = content;
    }
    if (min_content < 1e9) C_init = min_content;
    else C_init = 0.1;

    expFit->SetParameters(N0_init, 2.2, C_init);
    expFit->SetParLimits(0, 0, N0_init * 100);
    expFit->SetParLimits(1, 0.1, 20.0);
    expFit->SetParLimits(2, -C_init * 10, C_init * 10);
    expFit->SetParNames("N_{0}", "#tau", "C");
    expFit->SetLineColor(kRed);
    expFit->SetLineWidth(3);

    int fitStatus = h_dt_michel->Fit(expFit, "RE+", "SAME", FIT_MIN, FIT_MAX);
    expFit->Draw("SAME");
    
    gPad->Update();
    TPaveStats *stats = (TPaveStats*)h_dt_michel->FindObject("stats");
    if (stats) {
        stats->SetX1NDC(0.6);
        stats->SetX2NDC(0.9);
        stats->SetY1NDC(0.6);
        stats->SetY2NDC(0.9);
        stats->SetTextColor(kRed);
        stats->Clear();
        stats->AddText("Single Exp Fit");
        if (fitStatus == 0) {
            stats->AddText(Form("#tau = %.4f #pm %.4f #mus", expFit->GetParameter(1), expFit->GetParError(1)));
            stats->AddText(Form("#chi^{2}/NDF = %.4f", expFit->GetChisquare() / expFit->GetNDF()));
        } else {
            stats->AddText(Form("Fit failed (status %d)", fitStatus));
        }
        stats->Draw();
    }
    
    if (fitStatus == 0) {
        double N0 = expFit->GetParameter(0);
        double N0_err = expFit->GetParError(0);
        double tau = expFit->GetParameter(1);
        double tau_err = expFit->GetParError(1);
        double C = expFit->GetParameter(2);
        double C_err = expFit->GetParError(2);
        double chi2 = expFit->GetChisquare();
        int ndf = expFit->GetNDF();
        double chi2_ndf = ndf > 0 ? chi2 / ndf : 0;

        cout << "Single Exponential Fit Results (Michel dt, " << FIT_MIN << "-" << FIT_MAX << " µs):\n";
        cout << "Fit Status: " << fitStatus << " (0 = success)\n";
        cout << Form("τ = %.4f ± %.4f µs", tau, tau_err) << endl;
        cout << Form("N₀ = %.1f ± %.1f", N0, N0_err) << endl;
        cout << Form("C = %.1f ± %.1f", C, C_err) << endl;
        cout << Form("χ²/NDF = %.4f", chi2_ndf) << endl;
    } else {
        cout << "Single exponential fit failed with status " << fitStatus << endl;
    }
    cout << "----------------------------------------" << endl;

    c->Update();
    plotName = OUTPUT_DIR + "/Michel_dt_single_exp.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Fit accidental dt with simple exponential (no constant term) to get tau_accidental
    c->Clear();
    h_dt_accidental->SetLineColor(kBlue);
    h_dt_accidental->SetLineWidth(2);
    h_dt_accidental->GetXaxis()->SetTitle("Time to previous event (Muon) [#mus]");
    h_dt_accidental->SetTitle("Accidental Events Time Difference");
    h_dt_accidental->Draw("HIST");

    TF1* expFit_accidental = new TF1("expFit_accidental", ExpFitNoC, FIT_MIN, FIT_MAX, 2);
    double integral_acc = h_dt_accidental->Integral(h_dt_accidental->FindBin(FIT_MIN), h_dt_accidental->FindBin(FIT_MAX));
    double bin_width_acc = h_dt_accidental->GetBinWidth(1);
    double N0_init_acc = integral_acc * bin_width_acc / (FIT_MAX - FIT_MIN);
    
    expFit_accidental->SetParameters(N0_init_acc, 2.2);
    expFit_accidental->SetParNames("N_{0}", "#tau");
    expFit_accidental->SetParLimits(0, 0, N0_init_acc * 100);
    expFit_accidental->SetParLimits(1, 0.1, 20.0);
    expFit_accidental->SetLineColor(kRed);
    expFit_accidental->SetLineWidth(3);

    int fitStatus_acc = h_dt_accidental->Fit(expFit_accidental, "RE+", "SAME", FIT_MIN, FIT_MAX);
    expFit_accidental->Draw("SAME");
    
    gPad->Update();
    TPaveStats *stats_acc = (TPaveStats*)h_dt_accidental->FindObject("stats");
    if (stats_acc) {
        stats_acc->SetX1NDC(0.6);
        stats_acc->SetX2NDC(0.9);
        stats_acc->SetY1NDC(0.6);
        stats_acc->SetY2NDC(0.9);
        stats_acc->SetTextColor(kRed);
        stats_acc->Clear();
        stats_acc->AddText("Accidental DeltaT");
        if (fitStatus_acc == 0) {
            stats_acc->AddText(Form("#tau = %.4f #pm %.4f #mus", expFit_accidental->GetParameter(1), expFit_accidental->GetParError(1)));
            stats_acc->AddText(Form("#chi^{2}/NDF = %.4f", expFit_accidental->GetChisquare() / expFit_accidental->GetNDF()));
        } else {
            stats_acc->AddText(Form("Fit failed (status %d)", fitStatus_acc));
        }
        stats_acc->Draw();
    }
    
    double tau_accidental = 0;
    if (fitStatus_acc == 0) {
        tau_accidental = expFit_accidental->GetParameter(1);
        cout << "Accidental Events Exponential Fit Results (" << FIT_MIN << "-" << FIT_MAX << " µs):\n";
        cout << "Fit Status: " << fitStatus_acc << " (0 = success)\n";
        cout << Form("τ = %.4f ± %.4f µs", tau_accidental, expFit_accidental->GetParError(1)) << endl;
        cout << Form("χ²/NDF = %.4f", expFit_accidental->GetChisquare() / expFit_accidental->GetNDF()) << endl;
    } else {
        cout << "Accidental events fit failed with status " << fitStatus_acc << endl;
        // Use default value if fit fails
        tau_accidental = 2.2;
    }
    cout << "----------------------------------------" << endl;

    c->Update();
    plotName = OUTPUT_DIR + "/Accidental_dt.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // New Michel dt plot with double exponential fit
    c->Clear();
    h_dt_michel_double->SetLineWidth(2);
    h_dt_michel_double->SetLineColor(kBlack);
    h_dt_michel_double->GetXaxis()->SetTitle("Time to previous event (Muon) [#mus]");
    h_dt_michel_double->SetTitle("Michel Electron Time Difference (Double Exp Fit)");
    h_dt_michel_double->Draw("HIST");

    TF1* doubleExpFit = new TF1("doubleExpFit", DoubleExpFit, FIT_MIN, FIT_MAX, 4);
    double integral_double = h_dt_michel_double->Integral(h_dt_michel_double->FindBin(FIT_MIN), h_dt_michel_double->FindBin(FIT_MAX));
    double bin_width_double = h_dt_michel_double->GetBinWidth(1);
    double N0_init_double = integral_double * bin_width_double * 0.8; // Assume 80% for N1
    double N2_init_double = integral_double * bin_width_double * 0.2; // Assume 20% for N2
    doubleExpFit->SetParameters(N0_init_double, 2.2, N2_init_double, tau_accidental);
    doubleExpFit->SetParNames("N_{1}", "#tau_{1}", "N_{2}", "#tau_{2}");
    doubleExpFit->SetParLimits(0, 0, integral_double * bin_width_double * 10);
    doubleExpFit->SetParLimits(1, 1.0, 5.0); // Tighter constraint to prevent tau1 < 1 µs
    doubleExpFit->SetParLimits(2, 0, integral_double * bin_width_double * 10);
    // Fix tau2 parameter to the accidental value we obtained
    doubleExpFit->FixParameter(3, tau_accidental);
    doubleExpFit->SetLineColor(kRed);
    doubleExpFit->SetLineWidth(3);

    int fitStatus_double = h_dt_michel_double->Fit(doubleExpFit, "RE+", "SAME", FIT_MIN, FIT_MAX);
    doubleExpFit->Draw("SAME");
    
    gPad->Update();
    TPaveStats *stats_double = (TPaveStats*)h_dt_michel_double->FindObject("stats");
    if (stats_double) {
        stats_double->SetX1NDC(0.6);
        stats_double->SetX2NDC(0.9);
        stats_double->SetY1NDC(0.6);
        stats_double->SetY2NDC(0.9);
        stats_double->SetTextColor(kRed);
        stats_double->Clear();
        stats_double->AddText("Double Exp Fit");
        if (fitStatus_double == 0) {
            stats_double->AddText(Form("#tau_{1} = %.4f #pm %.4f #mus", doubleExpFit->GetParameter(1), doubleExpFit->GetParError(1)));
            stats_double->AddText(Form("N_{1} = %.1f #pm %.1f", doubleExpFit->GetParameter(0), doubleExpFit->GetParError(0)));
            stats_double->AddText(Form("N_{2} = %.1f #pm %.1f", doubleExpFit->GetParameter(2), doubleExpFit->GetParError(2)));
            stats_double->AddText(Form("#tau_{2} = %.4f (fixed)", doubleExpFit->GetParameter(3)));
            stats_double->AddText(Form("#chi^{2}/NDF = %.4f", doubleExpFit->GetChisquare() / doubleExpFit->GetNDF()));
        } else {
            stats_double->AddText(Form("Fit failed (status %d)", fitStatus_double));
        }
        stats_double->Draw();
    }
    
    if (fitStatus_double == 0) {
        double N1 = doubleExpFit->GetParameter(0);
        double N1_err = doubleExpFit->GetParError(0);
        double tau1 = doubleExpFit->GetParameter(1);
        double tau1_err = doubleExpFit->GetParError(1);
        double N2 = doubleExpFit->GetParameter(2);
        double N2_err = doubleExpFit->GetParError(2);
        double tau2 = doubleExpFit->GetParameter(3);
        double chi2 = doubleExpFit->GetChisquare();
        int ndf = doubleExpFit->GetNDF();
        double chi2_ndf = ndf > 0 ? chi2 / ndf : 0;

        cout << "Double Exponential Fit Results (Michel dt, " << FIT_MIN << "-" << FIT_MAX << " µs):\n";
        cout << "Fit Status: " << fitStatus_double << " (0 = success)\n";
        cout << Form("τ₁ = %.4f ± %.4f µs", tau1, tau1_err) << endl;
        cout << Form("N₁ = %.1f ± %.1f", N1, N1_err) << endl;
        cout << Form("N₂ = %.1f ± %.1f", N2, N2_err) << endl;
        cout << Form("τ₂ = %.4f (fixed)", tau2) << endl;
        cout << Form("χ²/NDF = %.4f", chi2_ndf) << endl;
    } else {
        cout << "Double exponential fit failed with status " << fitStatus_double << endl;
    }
    cout << "----------------------------------------" << endl;

    c->Update();
    plotName = OUTPUT_DIR + "/Michel_dt_double_exp.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Compare different exponential fit start times for single exponential
    std::vector<double> fit_starts = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    std::vector<double> taus, tau_errs, chi2ndfs;
    int best_index = -1;
    double min_chi2ndf = 1e9;
    
    for (int i = 0; i < fit_starts.size(); i++) {
        double fit_start = fit_starts[i];
        double fit_end = 16.0;
        
        TF1* expFit_var = new TF1(Form("expFit_var_%.1f", fit_start), ExpFit, fit_start, fit_end, 3);
        
        double C_init_var = 0;
        int bin_12_var = h_dt_michel->FindBin(12.0);
        int bin_16_var = h_dt_michel->FindBin(16.0);
        double min_content_var = 1e9;
        for (int bin = bin_12_var; bin <= bin_16_var; bin++) {
            double content = h_dt_michel->GetBinContent(bin);
            if (content > 0 && content < min_content_var) min_content_var = content;
        }
        if (min_content_var < 1e9) C_init_var = min_content_var;
        else C_init_var = 0.1;
        
        double integral_var = h_dt_michel->Integral(h_dt_michel->FindBin(fit_start), h_dt_michel->FindBin(fit_end));
        double bin_width_var = h_dt_michel->GetBinWidth(1);
        double N0_init_var = (integral_var * bin_width_var - C_init_var * (fit_end - fit_start)) / 2.2;
        if (N0_init_var < 0) N0_init_var = 100;
        
        expFit_var->SetParameters(N0_init_var, 2.2, C_init_var);
        expFit_var->SetParNames("N_{0}", "#tau", "C");
        expFit_var->SetParLimits(0, 0, N0_init_var * 100);
        expFit_var->SetParLimits(1, 0.1, 20.0);
        expFit_var->SetParLimits(2, -C_init_var * 10, C_init_var * 10);
        
        int fitStatus_var = h_dt_michel->Fit(expFit_var, "QRN+", "", fit_start, fit_end);
        
        double tau = expFit_var->GetParameter(1);
        double tau_err = expFit_var->GetParError(1);
        double chi2 = expFit_var->GetChisquare();
        int ndf = expFit_var->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : 999;
        
        taus.push_back(tau);
        tau_errs.push_back(tau_err);
        chi2ndfs.push_back(chi2ndf);
        
        if (chi2ndf < min_chi2ndf && fitStatus_var == 0) {
            min_chi2ndf = chi2ndf;
            best_index = i;
        }
        
        cout << Form("Single Exp Fit Range %.1f-%.1f µs:\n", fit_start, fit_end);
        cout << "Fit Status: " << fitStatus_var << " (0 = success)\n";
        cout << Form("τ = %.4f ± %.4f µs", tau, tau_err) << endl;
        cout << Form("χ²/NDF = %.4f", chi2ndf) << endl;
        cout << "----------------------------------------" << endl;
        
        delete expFit_var;
    }
    
    if (best_index >= 0) {
        cout << Form("Best Single Exp Fit Range: %.1f-16.0 µs\n", fit_starts[best_index]);
        cout << Form("τ = %.4f ± %.4f µs", taus[best_index], tau_errs[best_index]) << endl;
        cout << Form("χ²/NDF = %.4f (minimum)", chi2ndfs[best_index]) << endl;
        cout << "----------------------------------------" << endl;
    }
    
    TCanvas* c_comp = new TCanvas("c_comp", "Single Exp Fit Start Time Comparison", 1200, 800);
    c_comp->SetGrid();
    
    TPad* pad = new TPad("pad", "pad", 0, 0, 1, 1);
    pad->Draw();
    pad->cd();
    
    TGraph* g_chi2 = new TGraph(fit_starts.size(), &fit_starts[0], &chi2ndfs[0]);
    TGraph* g_tau = new TGraph(fit_starts.size(), &fit_starts[0], &taus[0]);
    
    g_chi2->SetTitle("Single Exp Fit Start Time Comparison");
    g_chi2->GetXaxis()->SetTitle("Fit Start Time (#mus)");
    g_chi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
    g_chi2->SetMarkerStyle(20);
    g_chi2->SetMarkerColor(kBlue);
    g_chi2->SetLineColor(kBlue);
    g_chi2->SetLineWidth(2);
    
    g_tau->SetMarkerStyle(22);
    g_tau->SetMarkerColor(kRed);
    g_tau->SetLineColor(kRed);
    g_tau->SetLineWidth(2);
    
    g_chi2->Draw("APL");
    
    pad->Update();
    double ymin = pad->GetUymin();
    double ymax = pad->GetUymax();
    
    double tau_min = *min_element(taus.begin(), taus.end());
    double tau_max = *max_element(taus.begin(), taus.end());
    double scale = (ymax - ymin)/(tau_max - tau_min);
    double offset = ymin - tau_min * scale;
    
    for (int i = 0; i < g_tau->GetN(); i++) {
        double x, y;
        g_tau->GetPoint(i, x, y);
        g_tau->SetPoint(i, x, y * scale + offset);
    }
    
    g_tau->Draw("PL same");
    
    TGaxis* axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                             gPad->GetUxmax(), gPad->GetUymax(),
                             tau_min, tau_max, 510, "+L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetTitle("#tau (#mus)");
    axis->SetTitleColor(kRed);
    axis->Draw();
    
    TLegend* leg_comp = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_comp->AddEntry(g_chi2, "#chi^{2}/ndf", "lp");
    leg_comp->AddEntry(g_tau, "#tau", "lp");
    leg_comp->Draw();
    
    string compPlotName = OUTPUT_DIR + "/SingleExp_FitStartComparison.png";
    c_comp->SaveAs(compPlotName.c_str());
    cout << "Saved single exp comparison plot: " << compPlotName << endl;
    
    delete g_chi2;
    delete g_tau;
    delete leg_comp;
    delete axis;
    delete pad;
    delete c_comp;

    // Compare different double exponential fit start times
    std::vector<double> dbl_taus1, dbl_taus1_err, dbl_chi2ndfs;
    int dbl_best_index = -1;
    double dbl_min_chi2ndf = 1e9;
    
    for (int i = 0; i < fit_starts.size(); i++) {
        double fit_start = fit_starts[i];
        double fit_end = 16.0;
        
        TF1* dblExpFit_var = new TF1(Form("dblExpFit_var_%.1f", fit_start), DoubleExpFit, fit_start, fit_end, 4);
        
        double integral_var = h_dt_michel_double->Integral(h_dt_michel_double->FindBin(fit_start), h_dt_michel_double->FindBin(fit_end));
        double bin_width_var = h_dt_michel_double->GetBinWidth(1);
        
        dblExpFit_var->SetParameters(integral_var * bin_width_var * 0.8, 2.2, 
                                    integral_var * bin_width_var * 0.2, tau_accidental);
        dblExpFit_var->SetParNames("N_{1}", "#tau_{1}", "N_{2}", "#tau_{2}");
        dblExpFit_var->SetParLimits(0, 0, integral_var * bin_width_var * 10);
        dblExpFit_var->SetParLimits(1, 1.0, 5.0); // Tighter constraint for tau1
        dblExpFit_var->SetParLimits(2, 0, integral_var * bin_width_var * 10);
        dblExpFit_var->FixParameter(3, tau_accidental);
        
        int fitStatus_var = h_dt_michel_double->Fit(dblExpFit_var, "QRN+", "", fit_start, fit_end);
        
        double tau1 = dblExpFit_var->GetParameter(1);
        double tau1_err = dblExpFit_var->GetParError(1);
        double chi2 = dblExpFit_var->GetChisquare();
        int ndf = dblExpFit_var->GetNDF();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : 999;
        
        dbl_taus1.push_back(tau1);
        dbl_taus1_err.push_back(tau1_err);
        dbl_chi2ndfs.push_back(chi2ndf);
        
        if (chi2ndf < dbl_min_chi2ndf && fitStatus_var == 0) {
            dbl_min_chi2ndf = chi2ndf;
            dbl_best_index = i;
        }
        
        cout << Form("Double Exp Fit Range %.1f-%.1f µs:\n", fit_start, fit_end);
        cout << "Fit Status: " << fitStatus_var << " (0 = success)\n";
        cout << Form("τ₁ = %.4f ± %.4f µs", tau1, tau1_err) << endl;
        cout << Form("χ²/NDF = %.4f", chi2ndf) << endl;
        cout << "----------------------------------------" << endl;
        
        delete dblExpFit_var;
    }
    
    if (dbl_best_index >= 0) {
        cout << Form("Best Double Exp Fit Range: %.1f-16.0 µs\n", fit_starts[dbl_best_index]);
        cout << Form("τ₁ = %.4f ± %.4f µs", dbl_taus1[dbl_best_index], dbl_taus1_err[dbl_best_index]) << endl;
        cout << Form("χ²/NDF = %.4f (minimum)", dbl_chi2ndfs[dbl_best_index]) << endl;
        cout << "----------------------------------------" << endl;
    }
    
    TCanvas* c_dbl_comp = new TCanvas("c_dbl_comp", "Double Exp Fit Start Time Comparison", 1200, 800);
    c_dbl_comp->SetGrid();
    
    TPad* dbl_pad = new TPad("dbl_pad", "dbl_pad", 0, 0, 1, 1);
    dbl_pad->Draw();
    dbl_pad->cd();
    
    TGraph* g_dbl_chi2 = new TGraph(fit_starts.size(), &fit_starts[0], &dbl_chi2ndfs[0]);
    TGraph* g_dbl_tau1 = new TGraph(fit_starts.size(), &fit_starts[0], &dbl_taus1[0]);
    
    g_dbl_chi2->SetTitle("Double Exp Fit Start Time Comparison");
    g_dbl_chi2->GetXaxis()->SetTitle("Fit Start Time (#mus)");
    g_dbl_chi2->GetYaxis()->SetTitle("#chi^{2}/ndf");
    g_dbl_chi2->SetMarkerStyle(20);
    g_dbl_chi2->SetMarkerColor(kBlue);
    g_dbl_chi2->SetLineColor(kBlue);
    g_dbl_chi2->SetLineWidth(2);
    
    g_dbl_tau1->SetMarkerStyle(22);
    g_dbl_tau1->SetMarkerColor(kRed);
    g_dbl_tau1->SetLineColor(kRed);
    g_dbl_tau1->SetLineWidth(2);
    
    g_dbl_chi2->Draw("APL");
    
    dbl_pad->Update();
    double dbl_ymin = dbl_pad->GetUymin();
    double dbl_ymax = dbl_pad->GetUymax();
    
    double dbl_tau1_min = *min_element(dbl_taus1.begin(), dbl_taus1.end());
    double dbl_tau1_max = *max_element(dbl_taus1.begin(), dbl_taus1.end());
    
    double scale_tau1 = (dbl_ymax - dbl_ymin) / (dbl_tau1_max - dbl_tau1_min);
    double offset_tau1 = dbl_ymin - dbl_tau1_min * scale_tau1;
    
    for (int i = 0; i < g_dbl_tau1->GetN(); i++) {
        double x, y;
        g_dbl_tau1->GetPoint(i, x, y);
        g_dbl_tau1->SetPoint(i, x, y * scale_tau1 + offset_tau1);
    }
    
    g_dbl_tau1->Draw("PL same");
    
    TGaxis* axis_tau1 = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                                   gPad->GetUxmax(), gPad->GetUymax(),
                                   dbl_tau1_min, dbl_tau1_max, 510, "+L");
    axis_tau1->SetLineColor(kRed);
    axis_tau1->SetLabelColor(kRed);
    axis_tau1->SetTitle("#tau_{1} (#mus)");
    axis_tau1->SetTitleColor(kRed);
    axis_tau1->Draw();
    
    TLegend* leg_dbl_comp = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_dbl_comp->AddEntry(g_dbl_chi2, "#chi^{2}/ndf", "lp");
    leg_dbl_comp->AddEntry(g_dbl_tau1, "#tau_{1}", "lp");
    leg_dbl_comp->Draw();
    
    string dblCompPlotName = OUTPUT_DIR + "/DoubleExp_FitStartComparison.png";
    c_dbl_comp->SaveAs(dblCompPlotName.c_str());
    cout << "Saved double exp comparison plot: " << dblCompPlotName << endl;
    
    delete g_dbl_chi2;
    delete g_dbl_tau1;
    delete leg_dbl_comp;
    delete axis_tau1;
    delete dbl_pad;
    delete c_dbl_comp;

    c->Clear();
    h_energy_vs_dt->SetStats(0);
    h_energy_vs_dt->GetXaxis()->SetTitle("dt (#mus)");
    h_energy_vs_dt->Draw("COLZ");
    c->Update();
    plotName = OUTPUT_DIR + "/Michel_Energy_vs_dt.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    c->Clear();
    h_side_vp_muon->SetLineColor(kMagenta);
    h_side_vp_muon->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Side_Veto_Muon.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    c->Clear();
    h_top_vp_muon->SetLineColor(kCyan);
    h_top_vp_muon->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Top_Veto_Muon.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    c->Clear();
    h_trigger_bits->SetLineColor(kGreen);
    h_trigger_bits->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/TriggerBits_Distribution.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;
    
    createVetoPanelPlots(h_veto_panel, OUTPUT_DIR);

    c->Clear();
    h_accidental_energy->SetLineColor(kRed);
    h_accidental_energy->SetLineWidth(2);
    h_accidental_energy->GetXaxis()->SetTitle("Energy (p.e.)");
    h_accidental_energy->SetTitle("Accidental Events Energy Spectrum");
    h_accidental_energy->Draw("HIST");
    
    c->Update();
    plotName = OUTPUT_DIR + "/Accidental_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    delete h_muon_energy;
    delete h_michel_energy;
    delete h_dt_michel;
    delete h_dt_michel_double;
    delete h_energy_vs_dt;
    delete h_side_vp_muon;
    delete h_top_vp_muon;
    delete h_trigger_bits;
    for (int i = 0; i < 10; i++) {
        delete h_veto_panel[i];
    }
    
    delete h_dt_accidental;
    delete h_dt_michel_true;
    delete h_accidental_energy;
    
    delete c;

    cout << "Analysis complete. Results saved in " << OUTPUT_DIR << "/ (*.png)" << endl;
    return 0;
}
