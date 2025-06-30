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
#include <TLatex.h>
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
#include <cstdlib>

using std::cout;
using std::cerr;
using std::endl;
using namespace std;

// Constants
const int N_PMTS = 12;
const int PMT_CHANNEL_MAP[12] = {0,10,7,2,6,3,8,9,11,4,5,1};
const int PULSE_THRESHOLD = 30;
const int BS_UNCERTAINTY = 5;
const int EV61_THRESHOLD = 1200;
const double MUON_ENERGY_THRESHOLD = 50;
const double MICHEL_ENERGY_MIN = 40;
const double MICHEL_ENERGY_MAX = 1000;
const double MICHEL_ENERGY_MAX_DT = 400;
const double MICHEL_DT_MIN = 0.8;
const double MICHEL_DT_MAX = 16.0;
const int ADCSIZE = 45;
const double PEAK_POSITION_RMS_CUT = 2.5;
const double AREA_HEIGHT_RATIO_CUT = 1.2;
const int SATURATION_THRESHOLD_LOW = 10;     // Raw ADC < 10 indicates saturation
const int SATURATION_THRESHOLD_HIGH = 4085;  // Raw ADC > 4085 indicates saturation

// Generate unique output directory with timestamp
string getTimestamp() {
    time_t now = time(nullptr);
    struct tm tstruct;
    char buffer[20];
    tstruct = *localtime(&now);
    strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", &tstruct);
    return string(buffer);
}
string OUTPUT_DIR = "./AnalysisOutput_" + getTimestamp();

const std::vector<double> SIDE_VP_THRESHOLDS = {750, 950, 1200, 1375, 525, 700, 700, 500};
const double TOP_VP_THRESHOLD = 450;
const double FIT_MIN = 1.0;
const double FIT_MAX = 10.0;

// Pulse structures with proper initialization
struct pulse_temp {
    double start = 0;
    double end = 0;
    double peak = 0;
    double energy = 0;
    int peak_position = -1;
    bool is_saturated = false;  // Saturation flag
};

struct pulse {
    double start = 0;
    double end = 0;
    double peak = 0;
    double energy = 0;
    double number = 0;
    bool single = false;
    bool beam = false;
    double trigger = 0;
    double side_vp_energy = 0;
    double top_vp_energy = 0;
    double all_vp_energy = 0;
    double last_muon_time = 0;
    bool is_muon = false;
    bool is_michel = false;
    double peak_position_rms = 0;
    bool is_good_event = false;
    bool is_saturated = false;  // Saturation flag for the event
};

// Fitting functions for SPE calibration
Double_t fitGauss(Double_t *x, Double_t *par) {
    return par[0] * TMath::Gaus(x[0], par[1], par[2]);
}

Double_t six_fit_func(Double_t *x, Double_t *par) {
    return (par[0] * TMath::Gaus(x[0], par[1], par[2]) + 
            par[3] * TMath::Gaus(x[0], par[4], par[5]));
}

Double_t eight_fit_func(Double_t *x, Double_t *par) {
    return (par[0] * TMath::Gaus(x[0], par[1], par[2]) + 
            par[3] * TMath::Gaus(x[0], par[4], par[5]) + 
            par[6] * TMath::Gaus(x[0], 2.0 * par[4], TMath::Sqrt(2.0 * par[5]*par[5] - par[2]*par[2])) + 
            par[7] * TMath::Gaus(x[0], 3.0 * par[4], TMath::Sqrt(3.0 * par[5]*par[5] - 2.0 * par[2]*par[2])));
}

// Exponential fit function for Michel decay time
Double_t ExpFit(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0]/par[1]) + par[2];
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
    return static_cast<double>(most_common);
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

void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    if (data.empty()) {
        mean = 0;
        rms = 0;
        return;
    }
    mean = 0.0;
    for (const auto &value : data) mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data) rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

// Create output directory
bool createOutputDirectory(const string& dirName) {
    struct stat st;
    if (stat(dirName.c_str(), &st)) {
        if (mkdir(dirName.c_str(), 0755)) {
            cerr << "Error: Could not create directory " << dirName << endl;
            return false;
        }
        cout << "Created output directory: " << dirName << endl;
    } else {
        cout << "Output directory already exists: " << dirName << endl;
    }
    return true;
}

// Improved isGoodEvent function
bool isGoodEvent(const std::vector<pulse_temp>& pmt_pulses, const Double_t* mu1, const Double_t* baselineRMS) {
    int countAbove2PE = 0;
    for (int pmt = 0; pmt < N_PMTS; pmt++) {
        if (pmt_pulses[pmt].peak > 0 && pmt_pulses[pmt].peak > 2 * mu1[pmt]) {
            countAbove2PE++;
        }
    }

    if (countAbove2PE >= 3) {
        vector<Double_t> peakPositions;
        for (int pmt = 0; pmt < N_PMTS; pmt++) {
            if (pmt_pulses[pmt].peak > 0) {
                peakPositions.push_back(pmt_pulses[pmt].peak_position);
            }
        }
        
        if (!peakPositions.empty()) {
            Double_t dummyMean;
            Double_t current_rms;
            CalculateMeanAndRMS(peakPositions, dummyMean, current_rms);
            if (current_rms < PEAK_POSITION_RMS_CUT) return true;
        }
    } 
    else {
        int countConditionB = 0;
        for (int pmt = 0; pmt < N_PMTS; pmt++) {
            if (pmt_pulses[pmt].peak > 0) {
                if (pmt_pulses[pmt].peak > 3 * baselineRMS[PMT_CHANNEL_MAP[pmt]] && 
                    (pmt_pulses[pmt].energy / pmt_pulses[pmt].peak) > AREA_HEIGHT_RATIO_CUT) {
                    countConditionB++;
                }
            }
        }

        if (countConditionB >= 3) {
            vector<Double_t> peakPositions;
            for (int pmt = 0; pmt < N_PMTS; pmt++) {
                if (pmt_pulses[pmt].peak > 0) {
                    peakPositions.push_back(pmt_pulses[pmt].peak_position);
                }
            }
            
            if (!peakPositions.empty()) {
                Double_t dummyMean;
                Double_t current_rms;
                CalculateMeanAndRMS(peakPositions, dummyMean, current_rms);
                if (current_rms < PEAK_POSITION_RMS_CUT) return true;
            }
        }
    }

    return false;
}

// Improved SPE calibration function using multi-step fitting
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

    // Create SPE plots directory
    string speDir = OUTPUT_DIR + "/SPE_Fits";
    gSystem->mkdir(speDir.c_str(), kTRUE);

    TH1F *histArea[N_PMTS];
    Long64_t nLEDFlashes[N_PMTS] = {0};
    for (int i = 0; i < N_PMTS; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area", i + 1),
                               Form("PMT %d;ADC Counts;Events", i + 1), 150, -50, 400);
        histArea[i]->SetLineColor(kRed);
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

    // Suppress ROOT fit messages
    Int_t defaultErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

    TCanvas *c = new TCanvas("c", "SPE Fits", 1200, 800);
    c->Divide(4, 3); // Divide into 4x3 grid for 12 PMTs
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    
    for (int i = 0; i < N_PMTS; i++) {
        if (histArea[i]->GetEntries() < 1000) {
            cerr << "Warning: Insufficient data for PMT " << i + 1 << " in " << calibFileName << endl;
            mu1[i] = 0;
            mu1_err[i] = 0;
            delete histArea[i];
            continue;
        }

        c->cd(i+1);
        
        // Step 1: Fit single Gaussian
        TF1 *f1 = new TF1("f1", fitGauss, -50, 50, 3);
        f1->SetParameters(1500, 0, 25);
        f1->SetParNames("A0", "#mu_{0}", "#sigma_{0}");
        histArea[i]->Fit(f1, "Q", "", -50, 50);

        // Step 2: Fit 6-parameter function
        TF1 *f6 = new TF1("f6", six_fit_func, -50, 200, 6);
        f6->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), 
                          1800, 70, 30);
        f6->SetParNames("A0", "#mu_{0}", "#sigma_{0}", "A1", "#mu_{1}", "#sigma_{1}");
        histArea[i]->Fit(f6, "Q", "", -50, 200);

        // Step 3: Fit 8-parameter function
        TF1 *f8 = new TF1("f8", eight_fit_func, -50, 400, 8);
        f8->SetParameters(f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), 
                          f6->GetParameter(3), f6->GetParameter(4), f6->GetParameter(5), 
                          200, 50);
        f8->SetParNames("A0", "#mu_{0}", "#sigma_{0}", "A1", "#mu_{1}", "#sigma_{1}", "A2", "A3");
        f8->SetLineColor(kBlue);
        histArea[i]->Fit(f8, "Q", "", -50, 400);

        // Get results
        mu1[i] = f8->GetParameter(4);
        mu1_err[i] = f8->GetParError(4);

        // Draw
        histArea[i]->Draw();
        f8->Draw("same");

        TLatex tex;
        tex.SetTextFont(42);
        tex.SetTextSize(0.04);
        tex.SetNDC();
        tex.DrawLatex(0.15, 0.85, Form("PMT %d", i+1));
        tex.DrawLatex(0.15, 0.80, Form("mu1 = %.2f #pm %.2f", mu1[i], mu1_err[i]));
        
        delete f1;
        delete f6;
        delete f8;
    }

    // Save combined canvas
    string plotName = OUTPUT_DIR + "/SPE_Fits.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved SPE plot: " << plotName << endl;

    // Restore ROOT error level
    gErrorIgnoreLevel = defaultErrorLevel;
    
    // Clean up histograms
    for (int i = 0; i < N_PMTS; i++) {
        if (histArea[i]) delete histArea[i];
    }
    delete c;
    calibFile->Close();
}

int main(int argc, char *argv[]) {
    // Parse command-line arguments
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <calibration_file> <input_file1> [<input_file2> ...]" << endl;
        return -1;
    }

    string calibFileName = argv[1];
    vector<string> inputFiles;
    for (int i = 2; i < argc; i++) {
        inputFiles.push_back(argv[i]);
    }

    // Generate output directory name
    OUTPUT_DIR = "./AnalysisOutput_" + getTimestamp();
    
    // Create output directory
    if (!createOutputDirectory(OUTPUT_DIR)) {
        return 1;
    }

    cout << "Calibration file: " << calibFileName << endl;
    cout << "Input files:" << endl;
    for (const auto& file : inputFiles) {
        cout << "  " << file << endl;
    }

    // Check if calibration file exists
    if (gSystem->AccessPathName(calibFileName.c_str())) {
        cerr << "Error: Calibration file " << calibFileName << " not found" << endl;
        return -1;
    }

    // Check if at least one input file exists
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

    // Perform SPE calibration with improved method
    Double_t mu1[N_PMTS] = {0};
    Double_t mu1_err[N_PMTS] = {0};
    performCalibration(calibFileName, mu1, mu1_err);

    // Print calibration results
    cout << "SPE Calibration Results (from " << calibFileName << "):\n";
    for (int i = 0; i < N_PMTS; i++) {
        cout << "PMT " << i + 1 << ": mu1 = " << mu1[i] << " ± " << mu1_err[i] << " ADC counts/p.e.\n";
    }

    // Statistics counters
    int total_events = 0;
    int total_good_events = 0;
    int total_muons = 0;
    int total_michels = 0;
    int total_saturated_events = 0;

    // Map to track triggerBits counts
    std::map<int, int> trigger_counts;

    // Define histograms
    TH1D* h_muon_energy = new TH1D("muon_energy", "Muon Energy Distribution (with Michel Electrons);Energy (p.e.);Counts/100 p.e.", 550, -500, 5000);
    TH1D* h_michel_energy = new TH1D("michel_energy", "Michel Electron Energy Distribution;Energy (p.e.);Counts/4 p.e.", 200, 0, 800);
    TH1D* h_dt_michel = new TH1D("DeltaT", "Muon-Michel Time Difference ;Time to Previous event(Muon)(#mus);Counts/0.1 #mus", 160, 0, MICHEL_DT_MAX);
    TH2D* h_energy_vs_dt = new TH2D("energy_vs_dt", "Michel Energy vs Time Difference;dt (#mus);Energy (p.e.)", 160, 0, 1000, 200, 0, 2000);
    TH1D* h_side_vp_muon = new TH1D("side_vp_muon", "Side Veto Energy for Muons;Energy (ADC);Counts", 200, 0, 5000);
    TH1D* h_top_vp_muon = new TH1D("top_vp_muon", "Top Veto Energy for Muons;Energy (ADC);Counts", 200, 0, 1000);
    TH1D* h_trigger_bits = new TH1D("trigger_bits", "Trigger Bits Distribution;Trigger Bits;Counts", 36, 0, 36);
    TH1D* h_peak_position_rms = new TH1D("peak_position_rms", "Peak Position RMS Distribution;RMS (samples);Counts", 100, 0, 10);
    TH1D* h_good_vs_bad = new TH1D("good_vs_bad", "Event Quality;Quality;Counts", 2, 0, 2);
    TH1D* h_saturation = new TH1D("saturation", "Saturation Status;Status;Counts", 2, 0, 2);
    h_good_vs_bad->GetXaxis()->SetBinLabel(1, "Good Events");
    h_good_vs_bad->GetXaxis()->SetBinLabel(2, "Bad Events");
    h_saturation->GetXaxis()->SetBinLabel(1, "Saturated");
    h_saturation->GetXaxis()->SetBinLabel(2, "Not Saturated");

    for (const auto& inputFileName : inputFiles) {
        // Per-file statistics
        int file_events = 0;
        int file_good_events = 0;
        int file_muons = 0;
        int file_michels = 0;
        int file_saturated_events = 0;
        
        // Check if input file exists
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

        // Declaration of leaf types
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

        // Set branch addresses
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
        std::vector<std::pair<double, double>> muon_candidates;

        for (int iEnt = 0; iEnt < numEntries; iEnt++) {
            t->GetEntry(iEnt);
            file_events++;
            total_events++;

            // Fill triggerBits histogram and track counts
            h_trigger_bits->Fill(triggerBits);
            trigger_counts[triggerBits]++;

            // Initialize pulse
            struct pulse p;
            p.start = nsTime / 1000.0; // Convert ns to µs
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
            p.peak_position_rms = 0;
            p.is_good_event = false;
            p.is_saturated = false;

            std::vector<double> all_chan_start, all_chan_end, all_chan_peak, all_chan_energy;
            std::vector<double> side_vp_energy, top_vp_energy;
            std::vector<double> chan_starts_no_outliers;
            std::vector<pulse_temp> pmt_pulses(N_PMTS);
            TH1D h_wf("h_wf", "Waveform", ADCSIZE, 0, ADCSIZE);

            bool pulse_at_end = false;
            int pulse_at_end_count = 0;
            std::vector<double> veto_energies(10, 0); // Channels 12-21

            // First pass: check for saturation in any channel
            bool event_saturated = false;
            for (int iChan = 0; iChan < 23; iChan++) {
                for (int i = 0; i < ADCSIZE; i++) {
                    short rawADC = adcVal[iChan][i];
                    if (rawADC <= SATURATION_THRESHOLD_LOW || rawADC >= SATURATION_THRESHOLD_HIGH) {
                        event_saturated = true;
                        break;
                    }
                }
                if (event_saturated) break;
            }
            p.is_saturated = event_saturated;
            h_saturation->Fill(event_saturated ? 0 : 1);
            if (event_saturated) {
                file_saturated_events++;
                total_saturated_events++;
                // Skip saturated events for further processing
                continue;
            }

            for (int iChan = 0; iChan < 23; iChan++) {
                // Fill waveform histogram with baseline subtraction
                for (int i = 0; i < ADCSIZE; i++) {
                    h_wf.SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);
                }

                // Check beam status (channel 22)
                if (iChan == 22) {
                    double ev61_energy = 0;
                    for (int iBin = 1; iBin <= ADCSIZE; iBin++) {
                        ev61_energy += h_wf.GetBinContent(iBin);
                    }
                    if (ev61_energy > EV61_THRESHOLD) {
                        p.beam = true;
                    }
                }

                // Pulse detection
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
                            pt.start = thresholdBin * 16.0 / 1000.0; // Convert ns to µs
                            pt.peak = peak;
                            pt.end = iBin * 16.0 / 1000.0;
                            pt.peak_position = peakBin;
                            for (int j = peakBin - 1; j >= 1 && h_wf.GetBinContent(j) > BS_UNCERTAINTY; j--) {
                                if (h_wf.GetBinContent(j) > peak * 0.1) {
                                    pt.start = j * 16.0 / 1000.0;
                                }
                                pulseEnergy += h_wf.GetBinContent(j);
                            }
                            pt.energy = pulseEnergy;
                            
                            if (iChan <= 11) {
                                int pmt_idx = -1;
                                for (int k = 0; k < N_PMTS; k++) {
                                    if (PMT_CHANNEL_MAP[k] == iChan) {
                                        pmt_idx = k;
                                        break;
                                    }
                                }
                                if (pmt_idx >= 0) {
                                    pmt_pulses[pmt_idx] = pt;
                                    if (mu1[pmt_idx] > 0) {
                                        pt.energy /= mu1[pmt_idx];
                                        pt.peak /= mu1[pmt_idx];
                                    }
                                    all_chan_start.push_back(pt.start);
                                    all_chan_end.push_back(pt.end);
                                    all_chan_peak.push_back(pt.peak);
                                    all_chan_energy.push_back(pt.energy);
                                    if (pt.energy > 1) p.number += 1;
                                }
                            }
                            peak = 0;
                            peakBin = 0;
                            pulseEnergy = 0;
                            thresholdBin = 0;
                            onPulse = false;
                        }
                    }
                }

                // Store energy for veto panels (ADC)
                if (iChan >= 12 && iChan <= 19) {
                    side_vp_energy.push_back(allPulseEnergy);
                    veto_energies[iChan - 12] = allPulseEnergy;
                } else if (iChan >= 20 && iChan <= 21) {
                    double factor = (iChan == 20) ? 1.07809 : 1.0;
                    top_vp_energy.push_back(allPulseEnergy * factor);
                    veto_energies[iChan - 12] = allPulseEnergy * factor;
                }

                // Check for pulses at waveform end
                if (iChan <= 11 && h_wf.GetBinContent(ADCSIZE) > 100) {
                    pulse_at_end_count++;
                    if (pulse_at_end_count >= 10) pulse_at_end = true;
                }

                h_wf.Reset();
            }

            // Aggregate pulse properties
            if (!all_chan_start.empty()) {
                p.start += mostFrequent(all_chan_start);
            }
            if (!all_chan_end.empty()) {
                p.end += mostFrequent(all_chan_end);
            }
            p.energy = std::accumulate(all_chan_energy.begin(), all_chan_energy.end(), 0.0);
            p.peak = std::accumulate(all_chan_peak.begin(), all_chan_peak.end(), 0.0);
            p.side_vp_energy = std::accumulate(side_vp_energy.begin(), side_vp_energy.end(), 0.0);
            p.top_vp_energy = std::accumulate(top_vp_energy.begin(), top_vp_energy.end(), 0.0);
            p.all_vp_energy = p.side_vp_energy + p.top_vp_energy;

            // Check timing consistency
            for (const auto& start : all_chan_start) {
                if (fabs(start - mostFrequent(all_chan_start)) < 10 * 16.0 / 1000.0) {
                    chan_starts_no_outliers.push_back(start);
                }
            }
            p.single = (variance(chan_starts_no_outliers) < 5 * 16.0 / 1000.0);

            // Calculate peak position RMS for this event
            vector<Double_t> peakPositions;
            for (const auto& pt : pmt_pulses) {
                if (pt.peak > 0) {
                    peakPositions.push_back(pt.peak_position);
                }
            }
            if (!peakPositions.empty()) {
                Double_t dummyMean;
                CalculateMeanAndRMS(peakPositions, dummyMean, p.peak_position_rms);
                h_peak_position_rms->Fill(p.peak_position_rms);
            }

            // Apply afterpulsing/crosstalk cuts
            p.is_good_event = isGoodEvent(pmt_pulses, mu1, baselineRMS);
            h_good_vs_bad->Fill(p.is_good_event ? 0 : 1);

            // Only process muon and Michel identification for good events
            if (p.is_good_event) {
                file_good_events++;
                total_good_events++;

                // Muon detection
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
                    file_muons++;
                    total_muons++;
                    muon_candidates.emplace_back(p.start, p.energy);
                    h_side_vp_muon->Fill(p.side_vp_energy);
                    h_top_vp_muon->Fill(p.top_vp_energy);
                }

                // Michel electron detection
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

                // Define common Michel electron criteria
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

                // Apply additional cut for dt and energy_vs_dt plots
                bool is_michel_for_dt = is_michel_candidate && p.energy <= MICHEL_ENERGY_MAX_DT;

                if (is_michel_candidate) {
                    p.is_michel = true;
                    file_michels++;
                    total_michels++;
                    michel_muon_times.insert(last_muon_time);
                    h_michel_energy->Fill(p.energy);
                }

                if (is_michel_for_dt) {
                    h_dt_michel->Fill(dt);
                    h_energy_vs_dt->Fill(dt, p.energy);
                }
            }

            p.last_muon_time = last_muon_time;
        }

        // Second pass: Fill h_muon_energy for muons associated with Michel electrons
        for (const auto& muon : muon_candidates) {
            if (michel_muon_times.find(muon.first) != michel_muon_times.end()) {
                h_muon_energy->Fill(muon.second);
            }
        }

        // Print stats to console
        cout << "File " << inputFileName << " Statistics:\n";
        cout << "Total Events: " << file_events << "\n";
        cout << "Saturated Events: " << file_saturated_events << "\n";
        cout << "Good Events (passed cuts): " << file_good_events << "\n";
        cout << "Muons Detected: " << file_muons << "\n";
        cout << "Michel Electrons Detected: " << file_michels << "\n";
        cout << "------------------------\n";

        f->Close();
    }

    // Print triggerBits distribution
    cout << "Trigger Bits Distribution (all files):\n";
    for (const auto& pair : trigger_counts) {
        cout << "Trigger " << pair.first << ": " << pair.second << " events\n";
    }
    cout << "------------------------\n";
    cout << "Total saturated events: " << total_saturated_events << endl;

    // Generate analysis plots
    TCanvas *c = new TCanvas("c", "Analysis Plots", 1200, 800);
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    // Muon Energy
    c->Clear();
    h_muon_energy->SetLineColor(kBlue);
    h_muon_energy->Draw();
    c->Update();
    string plotName = OUTPUT_DIR + "/Muon_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Michel Energy
    c->Clear();
    h_michel_energy->SetLineColor(kRed);
    h_michel_energy->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Michel_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Michel dt with exponential fit
    c->Clear();
    h_dt_michel->SetMarkerStyle(20);
    h_dt_michel->SetMarkerSize(1.0);
    h_dt_michel->GetXaxis()->SetTitle("Time to previous event(Muon)#mus");
    h_dt_michel->Draw("PE");

    if (h_dt_michel->GetEntries() > 5) {
        double integral = h_dt_michel->Integral(h_dt_michel->FindBin(FIT_MIN), h_dt_michel->FindBin(FIT_MAX));
        double bin_width = h_dt_michel->GetBinWidth(1);
        double N0_init = integral * bin_width / (FIT_MAX - FIT_MIN);
        double C_init = 0;
        int bin_14 = h_dt_michel->FindBin(14.0);
        int bin_16 = h_dt_michel->FindBin(16.0);
        double min_content = 1e9;
        for (int i = bin_14; i <= bin_16; i++) {
            double content = h_dt_michel->GetBinContent(i);
            if (content > 0 && content < min_content) min_content = content;
        }
        if (min_content < 1e9) C_init = min_content;
        else C_init = 0.1;

        TF1 *expFit = new TF1("expFit", ExpFit, FIT_MIN, FIT_MAX, 3);
        expFit->SetParameters(N0_init, 2.2, C_init);
        expFit->SetParLimits(0, 0, N0_init * 100);
        expFit->SetParLimits(1, 0.1, 20.0);
        expFit->SetParLimits(2, -C_init * 10, C_init * 10);
        expFit->SetParNames("N_{0}", "#tau", "C");
        expFit->SetNpx(1000);

        int fitStatus = h_dt_michel->Fit(expFit, "RE", "", FIT_MIN, FIT_MAX);
        expFit->SetLineColor(kGreen);
        expFit->SetLineWidth(3);
        expFit->Draw("same");
        gPad->Update();

        TPaveStats *stats = (TPaveStats*)h_dt_michel->FindObject("stats");
        if (!stats) {
            stats = new TPaveStats(0.60, 0.60, 0.90, 0.90, "brNDC");
            stats->SetName("stats");
            h_dt_michel->GetListOfFunctions()->Add(stats);
        }
        stats->SetTextColor(kRed);
        stats->SetX1NDC(0.60);
        stats->SetX2NDC(0.90);
        stats->SetY1NDC(0.60);
        stats->SetY2NDC(0.90);
        stats->Clear();
        stats->AddText("DeltaT");
        stats->AddText(Form("#tau = %.4f #pm %.4f #mus", expFit->GetParameter(1), expFit->GetParError(1)));
        stats->AddText(Form("#chi^{2}/NDF = %.4f", expFit->GetChisquare() / expFit->GetNDF()));
        stats->AddText(Form("N_{0} = %.1f #pm %.1f", expFit->GetParameter(0), expFit->GetParError(0)));
        stats->AddText(Form("C = %.1f #pm %.1f", expFit->GetParameter(2), expFit->GetParError(2)));
        stats->Draw();
        gPad->Update();

        double N0 = expFit->GetParameter(0);
        double N0_err = expFit->GetParError(0);
        double tau = expFit->GetParameter(1);
        double tau_err = expFit->GetParError(1);
        double C = expFit->GetParameter(2);
        double C_err = expFit->GetParError(2);
        double chi2 = expFit->GetChisquare();
        int ndf = expFit->GetNDF();
        double chi2_ndf = ndf > 0 ? chi2 / ndf : 0;

        cout << "Exponential Fit Results (Michel dt, 1.5-16 µs):\n";
        cout << Form("Fit Status: %d (0 = success)", fitStatus) << endl;
        cout << Form("N₀ = %.1f ± %.1f", N0, N0_err) << endl;
        cout << Form("τ = %.4f ± %.4f µs", tau, tau_err) << endl;
        cout << Form("C = %.1f ± %.1f", C, C_err) << endl;
        cout << Form("χ² = %.1f", chi2) << endl;
        cout << Form("NDF = %d", ndf) << endl;
        cout << Form("χ²/NDF = %.4f", chi2_ndf) << endl;
        cout << "----------------------------------------" << endl;

        if (fitStatus != 0) {
            cout << "Warning: Exponential fit failed for h_dt_michel (status = " << fitStatus << ")" << endl;
            cout << "Initial Parameters: N0 = " << N0_init << ", τ = 2.2 µs, C = " << C_init << endl;
        }
    }

    c->Update();
    c->Modified();
    c->RedrawAxis();
    plotName = OUTPUT_DIR + "/Michel_dt.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Energy vs dt
    c->Clear();
    h_energy_vs_dt->SetStats(0);
    h_energy_vs_dt->GetXaxis()->SetTitle("dt (#mus)");
    h_energy_vs_dt->Draw("COLZ");
    c->Update();
    plotName = OUTPUT_DIR + "/Michel_Energy_vs_dt.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Side Veto Muon
    c->Clear();
    h_side_vp_muon->SetLineColor(kMagenta);
    h_side_vp_muon->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Side_Veto_Muon.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Top Veto Muon
    c->Clear();
    h_top_vp_muon->SetLineColor(kCyan);
    h_top_vp_muon->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Top_Veto_Muon.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Trigger Bits Distribution
    c->Clear();
    h_trigger_bits->SetLineColor(kGreen);
    h_trigger_bits->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/TriggerBits_Distribution.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Peak Position RMS
    c->Clear();
    h_peak_position_rms->SetLineColor(kBlue);
    h_peak_position_rms->Draw();
    c->Update();
    plotName = OUTPUT_DIR + "/Peak_Position_RMS.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Good vs Bad Events
    c->Clear();
    h_good_vs_bad->SetFillColor(kBlue);
    h_good_vs_bad->Draw("BAR");
    c->Update();
    plotName = OUTPUT_DIR + "/Good_vs_Bad_Events.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Saturation Status
    c->Clear();
    h_saturation->SetFillColor(kOrange);
    h_saturation->Draw("BAR");
    c->Update();
    plotName = OUTPUT_DIR + "/Saturation_Status.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Clean up
    delete h_muon_energy;
    delete h_michel_energy;
    delete h_dt_michel;
    delete h_energy_vs_dt;
    delete h_side_vp_muon;
    delete h_top_vp_muon;
    delete h_trigger_bits;
    delete h_peak_position_rms;
    delete h_good_vs_bad;
    delete h_saturation;
    delete c;

    cout << "Analysis complete. Results saved in " << OUTPUT_DIR << "/ (*.png)" << endl;
    cout << "Total Events Processed: " << total_events << endl;
    cout << "Total Saturated Events: " << total_saturated_events << endl;
    cout << "Total Good Events: " << total_good_events << endl;
    cout << "Total Muons Detected: " << total_muons << endl;
    cout << "Total Michel Electrons Detected: " << total_michels << endl;
    return 0;
}
