#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TF1.h>
#include <TCanvas.h>
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
#include <random>
#include <memory>

using std::cout;
using std::cerr;
using std::endl;
using namespace std;

// Constants
const int N_PMTS = 12;
const int PMT_CHANNEL_MAP[12] = {0,10,7,2,6,3,8,9,11,4,5,1};
const int BS_UNCERTAINTY = 5;
const int EV61_THRESHOLD = 1200;
const double MUON_ENERGY_THRESHOLD = 50;
const double MICHEL_ENERGY_MIN = 40;
const double MICHEL_ENERGY_MAX = 1000;
const double MICHEL_ENERGY_MAX_DT = 500;
const double MICHEL_DT_MIN = 0.8;
const double MICHEL_DT_MAX = 16.0;
const int ADCSIZE = 45;
const double PEAK_POSITION_RMS_CUT = 2.5;
const double AREA_HEIGHT_RATIO_CUT = 1.2;
const int SATURATION_THRESHOLD_LOW = 10;
const int SATURATION_THRESHOLD_HIGH = 4000;

// Filtering parameters
const double CHI2_THRESHOLD = 2.0;
const int MAX_ITERATIONS = 3;

// Problematic time regions for special filtering
const double REGION1_LOW = 4.0;
const double REGION1_HIGH = 4.5;
const double REGION2_LOW = 5.5;
const double REGION2_HIGH = 6.0;
const double CHI2_THRESHOLD_SPECIAL = 1.0;

// Channel-specific trigger thresholds
const std::vector<double> TRIGGER_THRESHOLDS = {
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    50, 50, 50, 50, 50, 50, 50, 50,
    30, 30,
    100
};

// Combined veto thresholds
const std::vector<double> VETO_THRESHOLDS = {
    750, 950, 1200, 1375, 525, 700, 700, 500, 450, 450
};

// Fit ranges
const double FIT_MIN = 1.0;
const double FIT_MAX = 10.0;
const double FIT_MIN_FILTERED = 1.0;
const double FIT_MAX_FILTERED = 10.0;

// Generate unique output directory
string getTimestamp() {
    time_t now = time(nullptr);
    struct tm tstruct;
    char buffer[20];
    tstruct = *localtime(&now);
    strftime(buffer, sizeof(buffer), "%Y%m%d_%H%M%S", &tstruct);
    return string(buffer);
}
string OUTPUT_DIR = "./AnalysisOutput_" + getTimestamp();

// Pulse structures
struct pulse_temp {
    double start = 0;
    double end = 0;
    double peak = 0;
    double energy = 0;
    int peak_position = -1;
    bool is_saturated = false;
    double peak_time_ns = 0; // Added for timing analysis
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
    bool is_saturated = false;
};

// Michel candidate structure
struct MichelCandidate {
    double dt;
    double energy;
    int eventID;
    string fileName;
};

// Fitting functions
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

// Improved SPE calibration function
bool performCalibration(const string &calibFileName, Double_t *mu1, Double_t *mu1_err) {
    TFile *calibFile = TFile::Open(calibFileName.c_str());
    if (!calibFile || calibFile->IsZombie()) {
        cerr << "Error opening calibration file: " << calibFileName << endl;
        return false;
    }

    TTree *calibTree = (TTree*)calibFile->Get("tree");
    if (!calibTree) {
        cerr << "Error accessing tree in calibration file" << endl;
        calibFile->Close();
        delete calibFile;
        return false;
    }

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

    Int_t defaultErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kError;

    TCanvas *c = new TCanvas("c", "SPE Fits", 1200, 800);
    c->Divide(4, 3);
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
        
        TF1 *f1 = new TF1("f1", fitGauss, -50, 50, 3);
        f1->SetParameters(1500, 0, 25);
        f1->SetParNames("A0", "#mu_{0}", "#sigma_{0}");
        histArea[i]->Fit(f1, "Q", "", -50, 50);

        TF1 *f6 = new TF1("f6", six_fit_func, -50, 200, 6);
        f6->SetParameters(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), 
                          1800, 70, 30);
        f6->SetParNames("A0", "#mu_{0}", "#sigma_{0}", "A1", "#mu_{1}", "#sigma_{1}");
        histArea[i]->Fit(f6, "Q", "", -50, 200);

        TF1 *f8 = new TF1("f8", eight_fit_func, -50, 400, 8);
        f8->SetParameters(f6->GetParameter(0), f6->GetParameter(1), f6->GetParameter(2), 
                          f6->GetParameter(3), f6->GetParameter(4), f6->GetParameter(5), 
                          200, 50);
        f8->SetParNames("A0", "#mu_{0}", "#sigma_{0}", "A1", "#mu_{1}", "#sigma_{1}", "A2", "A3");
        f8->SetLineColor(kBlue);
        histArea[i]->Fit(f8, "Q", "", -50, 400);

        mu1[i] = f8->GetParameter(4);
        mu1_err[i] = f8->GetParError(4);

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

    string plotName = OUTPUT_DIR + "/SPE_Fits.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved SPE plot: " << plotName << endl;

    gErrorIgnoreLevel = defaultErrorLevel;
    
    for (int i = 0; i < N_PMTS; i++) {
        if (histArea[i]) delete histArea[i];
    }
    delete c;
    calibFile->Close();
    delete calibFile;
    return true;
}

void CreateStatsBox(TH1* hist, TF1* fitFunc = nullptr, 
                   const char* title = "DeltaT",
                   int textColor = kRed,
                   double x1 = 0.60, double x2 = 0.90,
                   double y1 = 0.60, double y2 = 0.90) {
    
    if (!gPad || !hist) return;
    
    // First remove any existing stats box
    hist->SetStats(0);
    gPad->Modified();
    gPad->Update();
    
    // Create new stats box
    TPaveStats* stats = new TPaveStats(x1, y1, x2, y2, "brNDC");
    stats->SetName("custom_stats");
    stats->SetBorderSize(1);
    stats->SetFillStyle(1001);
    stats->SetFillColor(0); // Transparent
    stats->SetTextFont(42);
    stats->SetTextSize(0.030);
    stats->SetTextColor(textColor);
    stats->SetTextAlign(12); // Left-bottom aligned
    
    // Add title
    stats->AddText(title);
    
    // Add histogram statistics
    stats->AddText(Form("Entries = %d", (int)hist->GetEntries()));
    stats->AddText(Form("Mean = %.3f #mus", hist->GetMean()));
    stats->AddText(Form("Std Dev = %.3f #mus", hist->GetRMS()));
    
    // Add fit results if available
    if (fitFunc) {
        stats->AddText(Form("#chi^{2}/NDF = %.4f", 
                           fitFunc->GetChisquare() / fitFunc->GetNDF()));
        stats->AddText(Form("N_{0} = %.1f #pm %.1f", 
                           fitFunc->GetParameter(0), 
                           fitFunc->GetParError(0)));
       stats->AddText(Form("#tau = %.4f #pm %.4f #mus", 
                           fitFunc->GetParameter(1), 
                           fitFunc->GetParError(1)));                    
    }
    
    // Add to histogram and draw
    hist->GetListOfFunctions()->Add(stats);
    stats->Draw();
    gPad->Modified();
    gPad->Update();
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

    OUTPUT_DIR = "./AnalysisOutput_" + getTimestamp();
    if (!createOutputDirectory(OUTPUT_DIR)) {
        return 1;
    }

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
    if (!performCalibration(calibFileName, mu1, mu1_err)) {
        cerr << "Calibration failed. Exiting..." << endl;
        return 1;
    }

    cout << "SPE Calibration Results (from " << calibFileName << "):\n";
    for (int i = 0; i < N_PMTS; i++) {
        cout << "PMT " << i + 1 << ": mu1 = " << mu1[i] << " ± " << mu1_err[i] << " ADC counts/p.e.\n";
    }

    int total_events = 0;
    int total_good_events = 0;
    int total_muons = 0;
    int total_michels = 0;
    int total_saturated_events = 0;
    std::map<int, int> trigger_counts;

    // Define histograms
    TH1D* h_muon_energy = new TH1D("muon_energy", "Muon Energy Distribution (with Michel Electrons);Energy (p.e.);Counts/100 p.e.", 550, -500, 5000);
    TH1D* h_michel_energy = new TH1D("michel_energy", "Michel Electron Energy Distribution;Energy (p.e.);Counts/4 p.e.", 200, 0, 800);
    TH1D* h_michel_energy_filtered = new TH1D("michel_energyFiltered", "Michel Electron Energy Distribution;Energy (p.e.);Counts/4 p.e.", 200, 0, 800);
    TH1D* h_dt_michel = new TH1D("DeltaTInital", "Initial Muon-Michel Time Difference;Time to Previous event(Muon)(#mus);Counts/0.25 #mus", 64, 0, MICHEL_DT_MAX);
    TH1D* h_dt_michel_filtered = new TH1D("DeltaT", "Muon-Michel Time Difference;Time to Previous event(Muon)(#mus);Counts/0.25 #mus", 64, 0, MICHEL_DT_MAX);
    TH2D* h_energy_vs_dt = new TH2D("energy_vs_dt", "Michel Energy vs Time Difference;dt (#mus);Energy (p.e.)", 160, 0, 1000, 200, 0, 2000);
    TH1D* h_side_vp_muon = new TH1D("side_vp_muon", "Side Veto Energy for Muons;Energy (ADC);Counts", 200, 0, 5000);
    TH1D* h_top_vp_muon = new TH1D("top_vp_muon", "Top Veto Energy for Muons;Energy (ADC);Counts", 200, 0, 1000);
    TH1D* h_trigger_bits = new TH1D("trigger_bits", "Trigger Bits Distribution;Trigger Bits;Counts", 36, 0, 36);
    TH1D* h_peak_position_rms = new TH1D("peak_position_rms", "Peak Position RMS Distribution;RMS (samples);Counts", 100, 0, 10);
    TH1D* h_good_vs_bad = new TH1D("good_vs_bad", "Event Quality;Quality;Counts", 2, 0, 2);
    TH1D* h_saturation = new TH1D("saturation", "Saturation Status;Status;Counts", 2, 0, 2);
    TH1D* h_special_regions = new TH1D("h_special_regions", "Events in Problem Regions;Time (#mus);Counts", 
                                      (REGION1_HIGH-REGION1_LOW)*40, REGION1_LOW, REGION2_HIGH);
    h_good_vs_bad->GetXaxis()->SetBinLabel(1, "Good Events");
    h_good_vs_bad->GetXaxis()->SetBinLabel(2, "Bad Events");
    h_saturation->GetXaxis()->SetBinLabel(1, "Saturated");
    h_saturation->GetXaxis()->SetBinLabel(2, "Not Saturated");

    // NEW: Histograms for PMT timing analysis
    TH1D* h_dt_neither = new TH1D("h_dt_neither", "#Delta t between PMTs (neither saturated); #Delta t (ns); Counts", 200, -100, 100);
    TH1D* h_dt_one_sat = new TH1D("h_dt_one_sat", "#Delta t between PMTs (one saturated); #Delta t (ns); Counts", 200, -100, 100);
    TH1D* h_dt_both_sat = new TH1D("h_dt_both_sat", "#Delta t between PMTs (both saturated); #Delta t (ns); Counts", 200, -100, 100);
    TH2D* h_dt_vs_amplitude = new TH2D("h_dt_vs_amplitude", "#Delta t vs Amplitude; Amplitude (ADC); #Delta t (ns)", 100, 0, 5000, 100, -50, 50);

    // Global tracking
    double last_muon_time = -1000.0;
    std::map<double, double> global_muons;
    std::set<double> global_muon_times_with_michel;
    std::vector<MichelCandidate> michel_candidates;

    for (const auto& inputFileName : inputFiles) {
        int file_events = 0;
        int file_good_events = 0;
        int file_muons = 0;
        int file_michels = 0;
        int file_saturated_events = 0;

        std::unique_ptr<TFile> f(TFile::Open(inputFileName.c_str()));
        if (!f || f->IsZombie()) {
            cout << "Error: Cannot open file: " << inputFileName << ". Skipping..." << endl;
            continue;
        }

        TTree* t = dynamic_cast<TTree*>(f->Get("tree"));
        if (!t) {
            cout << "Error: Cannot find tree in file: " << inputFileName << endl;
            continue;
        }

        // Verify branches
        if (!t->GetBranch("eventID") || !t->GetBranch("nSamples") || !t->GetBranch("adcVal") ||
            !t->GetBranch("baselineMean") || !t->GetBranch("baselineRMS") || !t->GetBranch("pulseH") ||
            !t->GetBranch("peakPosition") || !t->GetBranch("area") || !t->GetBranch("nsTime") ||
            !t->GetBranch("triggerBits")) {
            cerr << "Error: Missing required branches in tree" << endl;
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

        for (int iEnt = 0; iEnt < numEntries; iEnt++) {
            t->GetEntry(iEnt);
            file_events++;
            total_events++;

            h_trigger_bits->Fill(triggerBits);
            trigger_counts[triggerBits]++;

            struct pulse p;
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
            std::vector<double> veto_energies(10, 0);

            // Check for saturation at event level
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
            }

            // Set event time using nsTime and verify triggerBits
            p.start = nsTime / 1000.0;
            p.end = p.start;
            bool trigger_found = false;
            for (int iChan = 0; iChan < 23; iChan++) {
                if (triggerBits & (1 << iChan)) {
                    for (int i = 0; i < ADCSIZE; i++) {
                        double iBinContent = adcVal[iChan][i] - baselineMean[iChan];
                        if (iBinContent >= TRIGGER_THRESHOLDS[iChan]) {
                            trigger_found = true;
                            break;
                        }
                    }
                }
            }
            if (!trigger_found) {
                cout << "Warning: No channel with set trigger bit exceeded threshold in event " << eventID << endl;
            }

            // Process pulses
            for (int iChan = 0; iChan < 23; iChan++) {
                if (pulseH[iChan] < TRIGGER_THRESHOLDS[iChan]) continue;

                for (int i = 0; i < ADCSIZE; i++) {
                    h_wf.SetBinContent(i + 1, adcVal[iChan][i] - baselineMean[iChan]);
                }

                if (iChan == 22) {
                    double ev61_energy = 0;
                    for (int iBin = 1; iBin <= ADCSIZE; iBin++) {
                        double iBinContent = h_wf.GetBinContent(iBin);
                        if (iBinContent >= TRIGGER_THRESHOLDS[iChan]) {
                            ev61_energy += iBinContent;
                        }
                    }
                    if (ev61_energy > EV61_THRESHOLD) {
                        p.beam = true;
                    }
                }

                bool onPulse = false;
                int thresholdBin = 0, peakBin = 0;
                double peak = 0, pulseEnergy = 0;
                double allPulseEnergy = 0;
                bool pmt_saturated = false; // Per-PMT saturation flag

                for (int iBin = 1; iBin <= ADCSIZE; iBin++) {
                    double iBinContent = h_wf.GetBinContent(iBin);
                    
                    // Check for saturation in this PMT
                    if (iBinContent <= SATURATION_THRESHOLD_LOW || iBinContent >= SATURATION_THRESHOLD_HIGH) {
                        pmt_saturated = true;
                    }
                    
                    if (iBinContent >= TRIGGER_THRESHOLDS[iChan] && !onPulse) {
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
                            // Calculate tpeak and adjust start time
                            double tpeak = (peakBin - thresholdBin) * 16.0; // in ns
                            pt.start = (nsTime + tpeak) / 1000.0; // in microseconds
                            pt.peak = peak;
                            pt.end = iBin * 16.0 / 1000.0;
                            pt.peak_position = peakBin;
                            pt.peak_time_ns = nsTime + (peakBin - 1) * 16.0; // Absolute peak time in ns
                            pt.is_saturated = pmt_saturated;
                            
                            for (int j = peakBin - 1; j >= 1 && h_wf.GetBinContent(j) > BS_UNCERTAINTY; j--) {
                                if (h_wf.GetBinContent(j) > peak * 0.1) {
                                    pt.start = j * 16.0 / 1000.0; // Adjust start if earlier significant content
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
                            if (iChan >= 12 && iChan <= 21) {
                                for (int j = thresholdBin; j <= iBin; j++) {
                                    allPulseEnergy += h_wf.GetBinContent(j);
                                }
                            }
                            peak = 0;
                            peakBin = 0;
                            pulseEnergy = 0;
                            thresholdBin = 0;
                            onPulse = false;
                            pmt_saturated = false;
                        }
                    }
                }

                if (iChan >= 12 && iChan <= 19) {
                    side_vp_energy.push_back(allPulseEnergy);
                    veto_energies[iChan - 12] = allPulseEnergy;
                } else if (iChan >= 20 && iChan <= 21) {
                    double factor = (iChan == 20) ? 1.07809 : 1.0;
                    top_vp_energy.push_back(allPulseEnergy * factor);
                    veto_energies[iChan - 12] = allPulseEnergy * factor;
                }

                if (iChan <= 11 && h_wf.GetBinContent(ADCSIZE) > 100) {
                    pulse_at_end_count++;
                    if (pulse_at_end_count >= 10) pulse_at_end = true;
                }

                h_wf.Reset();
            }

            if (!all_chan_end.empty()) {
                p.end = p.start + mostFrequent(all_chan_end);
            }
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

            p.is_good_event = isGoodEvent(pmt_pulses, mu1, baselineRMS);
            h_good_vs_bad->Fill(p.is_good_event ? 0 : 1);

            // NEW: PMT timing analysis for saturated events
            // Collect active PMTs with significant pulses (pulseH > threshold)
            std::vector<int> active_pmts;
            for (int i = 0; i < N_PMTS; i++) {
                if (pmt_pulses[i].peak > TRIGGER_THRESHOLDS[i]) {
                    active_pmts.push_back(i);
                }
            }
            
            // Only analyze events with at least two active PMTs
            if (active_pmts.size() >= 2) {
                // Consider all pairs of active PMTs
                for (size_t i_idx = 0; i_idx < active_pmts.size(); i_idx++) {
                    for (size_t j_idx = i_idx+1; j_idx < active_pmts.size(); j_idx++) {
                        int i = active_pmts[i_idx];
                        int j = active_pmts[j_idx];
                        
                        double t1 = pmt_pulses[i].peak_time_ns;
                        double t2 = pmt_pulses[j].peak_time_ns;
                        double dt = t1 - t2;
                        
                        bool sat_i = pmt_pulses[i].is_saturated;
                        bool sat_j = pmt_pulses[j].is_saturated;
                        
                        // Fill appropriate histogram based on saturation status
                        if (!sat_i && !sat_j) {
                            h_dt_neither->Fill(dt);
                        } else if (sat_i && sat_j) {
                            h_dt_both_sat->Fill(dt);
                        } else {
                            h_dt_one_sat->Fill(dt);
                        }
                        
                        // Fill 2D histogram for correlation
                        h_dt_vs_amplitude->Fill(pmt_pulses[i].peak, dt);
                        h_dt_vs_amplitude->Fill(pmt_pulses[j].peak, -dt);
                    }
                }
            }

            if (p.is_good_event) {
                file_good_events++;
                total_good_events++;

                bool veto_hit = false;
                for (int i = 0; i < 10; i++) {
                    if (veto_energies[i] > VETO_THRESHOLDS[i]) {
                        veto_hit = true;
                        break;
                    }
                }

                if ((p.energy > MUON_ENERGY_THRESHOLD && veto_hit) ||
                    (pulse_at_end && p.energy > MUON_ENERGY_THRESHOLD / 2 && veto_hit)) {
                    p.is_muon = true;
                    last_muon_time = p.start;
                    global_muons[p.start] = p.energy;
                    file_muons++;
                    total_muons++;
                    h_side_vp_muon->Fill(p.side_vp_energy);
                    h_top_vp_muon->Fill(p.top_vp_energy);
                }

                double dt = p.start - last_muon_time;
                bool veto_low = true;
                for (int i = 0; i < 10; i++) {
                    if (veto_energies[i] > VETO_THRESHOLDS[i]) {
                        veto_low = false;
                        break;
                    }
                }

                bool is_michel_candidate = p.energy >= MICHEL_ENERGY_MIN &&
                                          p.energy <= MICHEL_ENERGY_MAX &&
                                          dt >= MICHEL_DT_MIN &&
                                          dt <= MICHEL_DT_MAX &&
                                          p.number >= 10 &&
                                          veto_low &&
                                          p.trigger != 1 &&
                                          p.trigger != 4 &&
                                          p.trigger != 8 &&
                                          p.trigger != 16;

                bool is_michel_for_dt = is_michel_candidate && p.energy <= MICHEL_ENERGY_MAX_DT;

                if (is_michel_candidate) {
                    p.is_michel = true;
                    file_michels++;
                    total_michels++;
                    global_muon_times_with_michel.insert(last_muon_time);
                    h_michel_energy->Fill(p.energy);
                }

                if (is_michel_for_dt) {
                    MichelCandidate candidate;
                    candidate.dt = dt;
                    candidate.energy = p.energy;
                    candidate.eventID = eventID;
                    candidate.fileName = inputFileName;
                    michel_candidates.push_back(candidate);
                    h_dt_michel->Fill(dt);
                    h_energy_vs_dt->Fill(dt, p.energy);
                    
                    // Fill filtered energy spectrum if in good time window
                    if (dt >= MICHEL_DT_MIN && dt <= MICHEL_DT_MAX) {
                        h_michel_energy_filtered->Fill(p.energy);
                    }
                    
                    if ((dt >= REGION1_LOW && dt <= REGION1_HIGH) || 
                        (dt >= REGION2_LOW && dt <= REGION2_HIGH)) {
                        h_special_regions->Fill(dt);
                    }
                }
            }

            p.last_muon_time = last_muon_time;
        }

        cout << "File " << inputFileName << " Statistics:\n";
        cout << "Total Events: " << file_events << "\n";
        cout << "Saturated Events: " << file_saturated_events << "\n";
        cout << "Good Events (passed cuts): " << file_good_events << "\n";
        cout << "Muons Detected: " << file_muons << "\n";
        cout << "Michel Electrons Detected: " << file_michels << "\n";
        cout << "------------------------\n";
    }

    // Fill muon energy histogram
    for (const auto& muon : global_muons) {
        if (global_muon_times_with_michel.find(muon.first) != global_muon_times_with_michel.end()) {
            h_muon_energy->Fill(muon.second);
        }
    }

    // Print triggerBits distribution
    cout << "Trigger Bits Distribution (all files):\n";
    for (const auto& pair : trigger_counts) {
        cout << "Trigger " << pair.first << ": " << pair.second << " events\n";
    }
    cout << "Total saturated events: " << total_saturated_events << endl;

    // Fill h_dt_michel_filtered with events from 0.8-16 µs
    h_dt_michel_filtered->Reset();
    h_michel_energy_filtered->Reset();  // Reset the filtered energy spectrum too
    std::ofstream logFile((OUTPUT_DIR + "/excluded_events.txt").c_str());
    for (const auto& candidate : michel_candidates) {
        if (candidate.dt >= MICHEL_DT_MIN && candidate.dt <= MICHEL_DT_MAX) {
            h_dt_michel_filtered->Fill(candidate.dt);
            h_michel_energy_filtered->Fill(candidate.energy);  // Fill filtered energy
        } else {
            logFile << "Excluded (outside plotting range): File = " << candidate.fileName 
                    << ", EventID = " << candidate.eventID 
                    << ", dt = " << candidate.dt 
                    << ", energy = " << candidate.energy << endl;
        }
    }

    // Initial fit to identify high chi-square bins in 1-10 µs
    TCanvas *c = new TCanvas("c", "Analysis Plots", 1200, 800);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::vector<MichelCandidate> filtered_candidates = michel_candidates;
    double filtered_chi2 = 0, filtered_ndf = 0, filtered_n0 = 0, filtered_n0err = 0, filtered_tau = 0, filtered_tauerr = 0;
    
    // Iterative filtering loop
    int iteration = 0;
    double prev_chi2ndf = 1e9;
    bool continue_filtering = true;
    
    while (continue_filtering && iteration < MAX_ITERATIONS) {
        iteration++;
        cout << "\n=== Iteration " << iteration << " of event filtering ===" << endl;
        
        if (h_dt_michel_filtered->GetEntries() < 10) {
            cout << "Warning: Insufficient entries (" << h_dt_michel_filtered->GetEntries() 
                 << ") for exponential fit. Skipping iteration." << endl;
            break;
        }
        
        double integral = h_dt_michel_filtered->Integral(h_dt_michel_filtered->FindBin(FIT_MIN_FILTERED), 
                                                        h_dt_michel_filtered->FindBin(FIT_MAX_FILTERED));
        double bin_width = h_dt_michel_filtered->GetBinWidth(1);
        double N0_init = integral * bin_width / (FIT_MAX_FILTERED - FIT_MIN_FILTERED);
        double C_init = 0;
        int bin_14 = h_dt_michel_filtered->FindBin(14.0);
        int bin_16 = h_dt_michel_filtered->FindBin(16.0);
        double min_content = 1e9;
        for (int i = bin_14; i <= bin_16; i++) {
            double content = h_dt_michel_filtered->GetBinContent(i);
            if (content > 0 && content < min_content) min_content = content;
        }
        if (min_content < 1e9) C_init = min_content;
        else C_init = 0.1;

        TH1D *h_log = (TH1D*)h_dt_michel_filtered->Clone("h_log");
        for (int i = 1; i <= h_log->GetNbinsX(); i++) {
            double content = h_log->GetBinContent(i);
            double error = h_log->GetBinError(i);
            if (content > C_init && error > 0) {
                h_log->SetBinContent(i, log(content - C_init));
                h_log->SetBinError(i, error / (content - C_init));
            } else {
                h_log->SetBinContent(i, 0);
                h_log->SetBinError(i, 0);
            }
        }
        TF1 *linearFit = new TF1("linearFit", "[0] - x/[1]", FIT_MIN_FILTERED, FIT_MAX_FILTERED);
        linearFit->SetParameters(log(N0_init), 2.2);
        h_log->Fit(linearFit, "Q", "", FIT_MIN_FILTERED, FIT_MAX_FILTERED);
        double tau_init = linearFit->GetParameter(1);
        delete linearFit;
        delete h_log;

        TF1 *expFit = new TF1("expFit", ExpFit, FIT_MIN_FILTERED, FIT_MAX_FILTERED, 3);
        expFit->SetParameters(N0_init, tau_init, C_init);
        expFit->SetParLimits(0, 0, N0_init * 100);
        expFit->SetParLimits(1, 0.1, 20.0);
        expFit->SetParLimits(2, -C_init * 10, C_init * 10);
        expFit->SetParNames("N_{0}", "#tau", "C");
        expFit->SetNpx(1000);

        h_dt_michel_filtered->Fit(expFit, "RE", "", FIT_MIN_FILTERED, FIT_MAX_FILTERED);

        // Calculate chi-square contributions in 1-10 µs
        TH1D *h_chi2_contrib_initial = new TH1D("h_chi2_contrib_initial", 
                                               "Initial Chi-square Contribution;Time (#mus);(Data - Fit)^2/#sigma^2", 
                                               h_dt_michel_filtered->GetNbinsX(), 
                                               h_dt_michel_filtered->GetXaxis()->GetXmin(), 
                                               h_dt_michel_filtered->GetXaxis()->GetXmax());
        double chi2_initial = 0;
        int ndf_initial = 0;
        std::vector<int> bad_bins;
        std::vector<double> expected_counts;
        for (int i = h_dt_michel_filtered->FindBin(FIT_MIN_FILTERED); i <= h_dt_michel_filtered->FindBin(FIT_MAX_FILTERED); i++) {
            double data = h_dt_michel_filtered->GetBinContent(i);
            double center = h_dt_michel_filtered->GetBinCenter(i);
            double fit_val = expFit->Eval(center);
            double error = h_dt_michel_filtered->GetBinError(i);
            expected_counts.push_back(fit_val);
            if (error > 0) {
                double chi2_contrib = pow(data - fit_val, 2) / (error * error);
                h_chi2_contrib_initial->SetBinContent(i, chi2_contrib);
                chi2_initial += chi2_contrib;
                ndf_initial++;
                
                double chi2_threshold_used = CHI2_THRESHOLD;
                if ((center >= REGION1_LOW && center <= REGION1_HIGH) ||
                    (center >= REGION2_LOW && center <= REGION2_HIGH)) {
                    chi2_threshold_used = CHI2_THRESHOLD_SPECIAL;
                }
                
                if (chi2_contrib > chi2_threshold_used) {
                    bad_bins.push_back(i);
                }
            }
        }
        ndf_initial -= 3;

        TCanvas *c_chi2_initial = new TCanvas("c_chi2_initial", "Initial Chi-square Contributions", 1200, 800);
        h_chi2_contrib_initial->SetLineColor(kRed);
        h_chi2_contrib_initial->Draw();
        string chi2PlotName = OUTPUT_DIR + "/Chi2_Contributions_initial_iter" + std::to_string(iteration) + ".png";
        c_chi2_initial->SaveAs(chi2PlotName.c_str());
        cout << "Saved plot: " << chi2PlotName << endl;
        delete h_chi2_contrib_initial;
        delete c_chi2_initial;

        cout << "Fit Results (Iteration " << iteration << "):\n";
        cout << "χ² = " << chi2_initial << endl;
        cout << "NDF = " << ndf_initial << endl;
        cout << "χ²/NDF = " << chi2_initial / ndf_initial << endl;
        cout << "Bad bins (chi2 contribution > " << CHI2_THRESHOLD << " or > " << CHI2_THRESHOLD_SPECIAL << " in special regions): " << bad_bins.size() << endl;

        // Filter events based on chi-square contributions in 1-10 µs
        std::random_device rd;
        std::mt19937 gen(rd());
        int events_removed = 0;
        
        for (int bin : bad_bins) {
            double bin_center = h_dt_michel_filtered->GetBinCenter(bin);
            double bin_width = h_dt_michel_filtered->GetBinWidth(1);
            double bin_low = bin_center - bin_width / 2;
            double bin_high = bin_center + bin_width / 2;
            double observed = h_dt_michel_filtered->GetBinContent(bin);
            double expected = expFit->Eval(bin_center);
            int events_to_remove = static_cast<int>(observed - expected);
            if (events_to_remove <= 0) continue;

            std::vector<size_t> bin_candidates;
            for (size_t i = 0; i < filtered_candidates.size(); i++) {
                if (filtered_candidates[i].dt >= bin_low && filtered_candidates[i].dt < bin_high) {
                    bin_candidates.push_back(i);
                }
            }

            if (bin_candidates.size() > static_cast<size_t>(events_to_remove)) {
                std::shuffle(bin_candidates.begin(), bin_candidates.end(), gen);
                for (int i = 0; i < events_to_remove; i++) {
                    size_t idx = bin_candidates[i];
                    if ((bin_center >= REGION1_LOW && bin_center <= REGION1_HIGH) ||
                        (bin_center >= REGION2_LOW && bin_center <= REGION2_HIGH)) {
                        logFile << "SPECIAL REGION REMOVAL: Iteration " << iteration 
                                << ", File = " << filtered_candidates[idx].fileName 
                                << ", EventID = " << filtered_candidates[idx].eventID 
                                << ", dt = " << filtered_candidates[idx].dt 
                                << ", energy = " << filtered_candidates[idx].energy 
                                << ", bin = " << bin << endl;
                    } else {
                        logFile << "Excluded (high chi2): Iteration " << iteration 
                                << ", File = " << filtered_candidates[idx].fileName 
                                << ", EventID = " << filtered_candidates[idx].eventID 
                                << ", dt = " << filtered_candidates[idx].dt 
                                << ", energy = " << filtered_candidates[idx].energy 
                                << ", bin = " << bin << endl;
                    }
                    filtered_candidates.erase(filtered_candidates.begin() + idx);
                    for (size_t j = i; j < bin_candidates.size(); j++) {
                        if (bin_candidates[j] > idx) bin_candidates[j]--;
                    }
                }
                events_removed += events_to_remove;
            }
        }
        
        cout << "Removed " << events_removed << " events in iteration " << iteration << endl;
        delete expFit;

        // Refill histogram with remaining events
        h_dt_michel_filtered->Reset();
        h_michel_energy_filtered->Reset();  // Reset the filtered energy spectrum too
        for (const auto& candidate : filtered_candidates) {
            if (candidate.dt >= MICHEL_DT_MIN && candidate.dt <= MICHEL_DT_MAX) {
                h_dt_michel_filtered->Fill(candidate.dt);
                h_michel_energy_filtered->Fill(candidate.energy);  // Fill filtered energy
            }
        }
        
        // Check if we should continue filtering
        double current_chi2ndf = (ndf_initial > 0) ? chi2_initial / ndf_initial : 0;
        double improvement = prev_chi2ndf - current_chi2ndf;
        
        cout << "Current χ²/NDF: " << current_chi2ndf << " (Previous: " 
             << prev_chi2ndf << ", Improvement: " << improvement << ")" << endl;
             
        if (improvement > 0.5 || current_chi2ndf > 2.0) {
            prev_chi2ndf = current_chi2ndf;
        } else {
            continue_filtering = false;
            cout << "Stopping iterations - insufficient improvement" << endl;
        }
    }

    logFile.close();
    cout << "Final entries in h_dt_michel_filtered: " << h_dt_michel_filtered->GetEntries() << endl;

    // ======================================================================
    // NEW: PMT Timing Analysis - Overlay Plot
    // ======================================================================
    cout << "\n=== PMT Timing Analysis ===" << endl;
    cout << "Neither saturated entries: " << h_dt_neither->GetEntries() << endl;
    cout << "One saturated entries: " << h_dt_one_sat->GetEntries() << endl;
    cout << "Both saturated entries: " << h_dt_both_sat->GetEntries() << endl;
    
    // Create overlay plot
    TCanvas* c_timing = new TCanvas("c_timing", "PMT Timing Comparison", 800, 600);
    
    // Set styles
    h_dt_neither->SetLineColor(kBlue);
    h_dt_neither->SetLineWidth(2);
    h_dt_neither->SetTitle("#Delta t between PMT pairs; #Delta t (ns); Counts");
    
    h_dt_one_sat->SetLineColor(kRed);
    h_dt_one_sat->SetLineWidth(2);
    h_dt_one_sat->SetLineStyle(2);
    
    h_dt_both_sat->SetLineColor(kGreen);
    h_dt_both_sat->SetLineWidth(2);
    h_dt_both_sat->SetLineStyle(3);
    
    // Find maximum for scaling
    double max_val = max(h_dt_neither->GetMaximum(), max(h_dt_one_sat->GetMaximum(), h_dt_both_sat->GetMaximum()));
    h_dt_neither->SetMaximum(max_val * 1.2);
    
    // Draw
    h_dt_neither->Draw("HIST");
    h_dt_one_sat->Draw("HIST SAME");
    h_dt_both_sat->Draw("HIST SAME");
    
    // Add legend
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.85);
    leg->AddEntry(h_dt_neither, "Neither saturated", "l");
    leg->AddEntry(h_dt_one_sat, "One saturated", "l");
    leg->AddEntry(h_dt_both_sat, "Both saturated", "l");
    leg->SetBorderSize(0);
    leg->Draw();
    
    // Add stats box
    TPaveText* stats = new TPaveText(0.7, 0.5, 0.9, 0.7, "NDC");
    stats->SetBorderSize(1);
    stats->SetFillColor(0);
    stats->AddText(Form("Neither: %d", (int)h_dt_neither->GetEntries()));
    stats->AddText(Form("One sat: %d", (int)h_dt_one_sat->GetEntries()));
    stats->AddText(Form("Both sat: %d", (int)h_dt_both_sat->GetEntries()));
    stats->AddText(Form("RMS (neither): %.1f ns", h_dt_neither->GetRMS()));
    stats->AddText(Form("RMS (one sat): %.1f ns", h_dt_one_sat->GetRMS()));
    stats->AddText(Form("RMS (both): %.1f ns", h_dt_both_sat->GetRMS()));
    stats->Draw();
    
    // Save plot
    string timingPlotName = OUTPUT_DIR + "/PMT_Timing_Overlay.png";
    c_timing->SaveAs(timingPlotName.c_str());
    cout << "Saved PMT timing overlay plot: " << timingPlotName << endl;
    
    // Plot Δt vs amplitude
    TCanvas* c_dt_amp = new TCanvas("c_dt_amp", "Delta t vs Amplitude", 800, 600);
    h_dt_vs_amplitude->SetStats(0);
    h_dt_vs_amplitude->Draw("COLZ");
    
    string dtAmpPlotName = OUTPUT_DIR + "/DeltaT_vs_Amplitude.png";
    c_dt_amp->SaveAs(dtAmpPlotName.c_str());
    cout << "Saved Δt vs amplitude plot: " << dtAmpPlotName << endl;
    
    // ======================================================================
    // Main plotting section
    // ======================================================================
    
    // Muon Energy
    c->Clear();
    h_muon_energy->SetLineColor(kBlue);
    h_muon_energy->Draw();
    CreateStatsBox(h_muon_energy, nullptr, "Muon Energy", kRed);
    c->Update();
    string plotName = OUTPUT_DIR + "/Muon_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Michel Energy
    c->Clear();
    h_michel_energy->SetLineColor(kRed);
    h_michel_energy->Draw();
    CreateStatsBox(h_michel_energy, nullptr, "Michel Energy", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Michel_Energy.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Filtered Michel Energy
    c->Clear();
    h_michel_energy_filtered->SetLineColor(kRed);
    h_michel_energy_filtered->SetLineWidth(2);
    h_michel_energy_filtered->Draw();
    CreateStatsBox(h_michel_energy_filtered, nullptr, "Filtered Michel Energy", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Michel_Energy_filtered.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Michel dt with exponential fit
    c->Clear();
    h_dt_michel->SetMarkerStyle(20);
    h_dt_michel->SetMarkerSize(1.0);
    h_dt_michel->GetXaxis()->SetTitle("Time to previous event (Muon) (#mus)");
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

        TH1D *h_log = (TH1D*)h_dt_michel->Clone("h_log");
        for (int i = 1; i <= h_log->GetNbinsX(); i++) {
            double content = h_log->GetBinContent(i);
            double error = h_log->GetBinError(i);
            if (content > C_init && error > 0) {
                h_log->SetBinContent(i, log(content - C_init));
                h_log->SetBinError(i, error / (content - C_init));
            } else {
                h_log->SetBinContent(i, 0);
                h_log->SetBinError(i, 0);
            }
        }
        TF1 *linearFit = new TF1("linearFit", "[0] - x/[1]", FIT_MIN, FIT_MAX);
        linearFit->SetParameters(log(N0_init), 2.2);
        h_log->Fit(linearFit, "Q", "", FIT_MIN, FIT_MAX);
        double tau_init = linearFit->GetParameter(1);
        delete linearFit;
        delete h_log;

        TF1 *expFit = new TF1("expFit", ExpFit, FIT_MIN, FIT_MAX, 3);
        expFit->SetParameters(N0_init, tau_init, C_init);
        expFit->SetParLimits(0, 0, N0_init * 100);
        expFit->SetParLimits(1, 0.1, 20.0);
        expFit->SetParLimits(2, -C_init * 10, C_init * 10);
        expFit->SetParNames("N_{0}", "#tau", "C");
        expFit->SetNpx(1000);

        h_dt_michel->Fit(expFit, "RE", "", FIT_MIN, FIT_MAX);
        expFit->SetLineColor(kGreen);
        expFit->SetLineWidth(3);
        expFit->Draw("same");

        // Create the stats box with fit parameters
        CreateStatsBox(h_dt_michel, expFit, "DeltaT", kRed);

        cout << "Exponential Fit Results (Michel dt, " << FIT_MIN << "-" << FIT_MAX << " µs):\n";
        cout << Form("N₀ = %.1f ± %.1f", expFit->GetParameter(0), expFit->GetParError(0)) << endl;
        cout << Form("τ = %.4f ± %.4f µs", expFit->GetParameter(1), expFit->GetParError(1)) << endl;
        cout << Form("C = %.1f ± %.1f", expFit->GetParameter(2), expFit->GetParError(2)) << endl;
        cout << Form("χ²/NDF = %.4f", expFit->GetChisquare() / expFit->GetNDF()) << endl;
        cout << "----------------------------------------" << endl;

        delete expFit;
    } else {
        cout << "Warning: Insufficient entries (" << h_dt_michel->GetEntries() << ") for exponential fit." << endl;
        CreateStatsBox(h_dt_michel, nullptr, "DeltaT", kRed);
    }

    c->Update();
    c->Modified();
    c->RedrawAxis();
    plotName = OUTPUT_DIR + "/Michel_dt.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Filtered dt fit after event removal
    c->Clear();
    h_dt_michel_filtered->SetMarkerStyle(20);
    h_dt_michel_filtered->SetMarkerSize(1.0);
    h_dt_michel_filtered->GetXaxis()->SetTitle("Time to previous event (Muon) (#mus)");
    h_dt_michel_filtered->GetXaxis()->SetRangeUser(0, MICHEL_DT_MAX);
    h_dt_michel_filtered->Draw("PE");

    if (h_dt_michel_filtered->GetEntries() > 5) {
        double integral = h_dt_michel_filtered->Integral(h_dt_michel_filtered->FindBin(FIT_MIN_FILTERED), 
                                                        h_dt_michel_filtered->FindBin(FIT_MAX_FILTERED));
        double bin_width = h_dt_michel_filtered->GetBinWidth(1);
        double N0_init = integral * bin_width / (FIT_MAX_FILTERED - FIT_MIN_FILTERED);
        double C_init = 0;
        int bin_14 = h_dt_michel_filtered->FindBin(14.0);
        int bin_16 = h_dt_michel_filtered->FindBin(16.0);
        double min_content = 1e9;
        for (int i = bin_14; i <= bin_16; i++) {
            double content = h_dt_michel_filtered->GetBinContent(i);
            if (content > 0 && content < min_content) min_content = content;
        }
        if (min_content < 1e9) C_init = min_content;
        else C_init = 0.1;

        TH1D *h_log = (TH1D*)h_dt_michel_filtered->Clone("h_log");
        for (int i = 1; i <= h_log->GetNbinsX(); i++) {
            double content = h_log->GetBinContent(i);
            double error = h_log->GetBinError(i);
            if (content > C_init && error > 0) {
                h_log->SetBinContent(i, log(content - C_init));
                h_log->SetBinError(i, error / (content - C_init));
            } else {
                h_log->SetBinContent(i, 0);
                h_log->SetBinError(i, 0);
            }
        }
        TF1 *linearFit = new TF1("linearFit", "[0] - x/[1]", FIT_MIN_FILTERED, FIT_MAX_FILTERED);
        linearFit->SetParameters(log(N0_init), 2.2);
        h_log->Fit(linearFit, "Q", "", FIT_MIN_FILTERED, FIT_MAX_FILTERED);
        double tau_init = linearFit->GetParameter(1);
        delete linearFit;
        delete h_log;

        TF1 *expFit = new TF1("expFit", ExpFit, FIT_MIN_FILTERED, FIT_MAX_FILTERED, 3);
        expFit->SetParameters(N0_init, tau_init, C_init);
        expFit->SetParLimits(0, 0, N0_init * 100);
        expFit->SetParLimits(1, 0.1, 20.0);
        expFit->SetParLimits(2, -C_init * 10, C_init * 10);
        expFit->SetParNames("N_{0}", "#tau", "C");
        expFit->SetNpx(1000);

        h_dt_michel_filtered->Fit(expFit, "RE", "", FIT_MIN_FILTERED, FIT_MAX_FILTERED);
        expFit->SetLineColor(kGreen);
        expFit->SetLineWidth(3);
        expFit->Draw("same");

        // Create the stats box with fit parameters
        CreateStatsBox(h_dt_michel_filtered, expFit, " DeltaT", kRed);

        cout << "Single-Exponential Fit Results (Filtered Michel dt, " << FIT_MIN_FILTERED << "-" << FIT_MAX_FILTERED << " µs):\n";
        cout << Form("N₀ = %.1f ± %.1f", expFit->GetParameter(0), expFit->GetParError(0)) << endl;
        cout << Form("τ = %.4f ± %.4f µs", expFit->GetParameter(1), expFit->GetParError(1)) << endl;
        cout << Form("C = %.1f ± %.1f", expFit->GetParameter(2), expFit->GetParError(2)) << endl;
        cout << Form("χ²/NDF = %.4f", expFit->GetChisquare() / expFit->GetNDF()) << endl;
        cout << "----------------------------------------" << endl;

        delete expFit;
    } else {
        cout << "Warning: Insufficient entries (" << h_dt_michel_filtered->GetEntries() << ") for filtered exponential fit." << endl;
        CreateStatsBox(h_dt_michel_filtered, nullptr, "Filtered DeltaT", kRed);
    }

    c->Update();
    c->Modified();
    c->RedrawAxis();
    plotName = OUTPUT_DIR + "/Michel_dt_filtered.png";
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
    CreateStatsBox(h_side_vp_muon, nullptr, "Side Veto Muon", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Side_Veto_Muon.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Top Veto Muon
    c->Clear();
    h_top_vp_muon->SetLineColor(kCyan);
    h_top_vp_muon->Draw();
    CreateStatsBox(h_top_vp_muon, nullptr, "Top Veto Muon", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Top_Veto_Muon.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Trigger Bits Distribution
    c->Clear();
    h_trigger_bits->SetLineColor(kGreen);
    h_trigger_bits->Draw();
    CreateStatsBox(h_trigger_bits, nullptr, "Trigger Bits", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/TriggerBits_Distribution.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Peak Position RMS
    c->Clear();
    h_peak_position_rms->SetLineColor(kBlue);
    h_peak_position_rms->Draw();
    CreateStatsBox(h_peak_position_rms, nullptr, "Peak Position RMS", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Peak_Position_RMS.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Good vs Bad Events
    c->Clear();
    h_good_vs_bad->SetFillColor(kBlue);
    h_good_vs_bad->Draw("BAR");
    CreateStatsBox(h_good_vs_bad, nullptr, "Event Quality", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Good_vs_Bad_Events.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Saturation Status
    c->Clear();
    h_saturation->SetFillColor(kOrange);
    h_saturation->Draw("BAR");
    CreateStatsBox(h_saturation, nullptr, "Saturation Status", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Saturation_Status.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Special Regions Distribution
    c->Clear();
    h_special_regions->SetLineColor(kMagenta);
    h_special_regions->Draw();
    CreateStatsBox(h_special_regions, nullptr, "Special Regions", kRed);
    c->Update();
    plotName = OUTPUT_DIR + "/Special_Regions.png";
    c->SaveAs(plotName.c_str());
    cout << "Saved plot: " << plotName << endl;

    // Save histograms to ROOT file
    TFile *outFile = new TFile((OUTPUT_DIR + "/analysis_output.root").c_str(), "RECREATE");
    h_muon_energy->Write();
    h_michel_energy->Write();
    h_michel_energy_filtered->Write();
    h_dt_michel->Write();
    h_dt_michel_filtered->Write();
    h_energy_vs_dt->Write();
    h_side_vp_muon->Write();
    h_top_vp_muon->Write();
    h_trigger_bits->Write();
    h_peak_position_rms->Write();
    h_good_vs_bad->Write();
    h_saturation->Write();
    h_special_regions->Write();
    
    // Write timing analysis histograms
    h_dt_neither->Write();
    h_dt_one_sat->Write();
    h_dt_both_sat->Write();
    h_dt_vs_amplitude->Write();
    
    outFile->Close();
    delete outFile;

    // Clean up
    delete h_muon_energy;
    delete h_michel_energy;
    delete h_michel_energy_filtered;
    delete h_dt_michel;
    delete h_dt_michel_filtered;
    delete h_energy_vs_dt;
    delete h_side_vp_muon;
    delete h_top_vp_muon;
    delete h_trigger_bits;
    delete h_peak_position_rms;
    delete h_good_vs_bad;
    delete h_saturation;
    delete h_special_regions;
    
    // Delete timing analysis histograms
    delete h_dt_neither;
    delete h_dt_one_sat;
    delete h_dt_both_sat;
    delete h_dt_vs_amplitude;
    
    delete c;
    delete c_timing;
    delete c_dt_amp;

    cout << "Analysis complete. Results saved in " << OUTPUT_DIR << "/ (*.png, analysis_output.root)" << endl;
    cout << "Total Events Processed: " << total_events << endl;
    cout << "Total Saturated Events: " << total_saturated_events << endl;
    cout << "Total Good Events: " << total_good_events << endl;
    cout << "Total Muons Detected: " << total_muons << endl;
    cout << "Total Michel Electrons Detected: " << total_michels << endl;
    
    // Print timing analysis summary
    cout << "\n===== PMT Timing Reliability Summary =====" << endl;
    cout << "Neither saturated: RMS = " << h_dt_neither->GetRMS() << " ns" << endl;
    cout << "One saturated: RMS = " << h_dt_one_sat->GetRMS() << " ns" << endl;
    cout << "Both saturated: RMS = " << h_dt_both_sat->GetRMS() << " ns" << endl;
    
    return 0;
}
