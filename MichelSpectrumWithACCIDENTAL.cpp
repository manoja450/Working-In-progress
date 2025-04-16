//Muon and Michel Electron Analysis Code with Separate Canvas for Accidental Fit
// This code analyses cosmic muon events and their subsequent Michel electron decay
// in a particle detector. It includes PMT calibration using LED data, muon event selection
// with SiPM shower rejection, Michel electron identification, muon lifetime measurement,
// energy spectrum analysis, multiple exponential fits, an overlayed summary plot, and now
// an additional exponential accidental background model.
// A separate canvas is created to visualize the accidental fit in the 10–20 µs region.

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TChain.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>
#include <ctime>
#include <cstdlib>
#include <TLine.h>
#include <TLatex.h>

using namespace std;

// Constants for detector channels and thresholds.
const int N_PMTS = 12;
const int N_SIPMS = 10;
const int PMT_CHANNEL_MAP[N_PMTS] = {0,10,7,2,6,3,8,9,11,4,5,1};
const int SIPM_CHANNEL_MAP[N_SIPMS] = {12,13,14,15,16,17,18,19,20,21};
const double SIPM_THRESHOLDS[N_SIPMS] = {800,800,1100,1200,550,600,650,450,600,650};
const double PMT_THRESHOLDS[N_PMTS] = {4800,6000,5000,6000,6000,4707,4500,3000,2000,5000,4500,4800};

// SPE fitting function (4-Gaussian model)
Double_t SPEfit(Double_t *x, Double_t *par) {
    Double_t term1 = par[0] * exp(-0.5 * pow((x[0]-par[1])/par[2], 2));
    Double_t term2 = par[3] * exp(-0.5 * pow((x[0]-par[4])/par[5], 2));
    Double_t term3 = par[6] * exp(-0.5 * pow((x[0]-sqrt(2)*par[4]) / sqrt(2*pow(par[5],2)-pow(par[2],2)), 2));
    Double_t term4 = par[7] * exp(-0.5 * pow((x[0]-sqrt(3)*par[4]) / sqrt(3*pow(par[5],2)-2*pow(par[2],2)), 2));
    return term1 + term2 + term3 + term4;
}

// Muon decay exponential function: f(x) = A * exp(-x/tau)
Double_t DecayFit(Double_t *x, Double_t *par) {
    return par[0] * exp(-x[0]/par[1]);
}

// Calculate mean and RMS of a vector.
void CalculateMeanAndRMS(const vector<Double_t> &data, Double_t &mean, Double_t &rms) {
    mean = 0.0;
    for (const auto &value : data)
        mean += value;
    mean /= data.size();
    
    rms = 0.0;
    for (const auto &value : data)
        rms += pow(value - mean, 2);
    rms = sqrt(rms / data.size());
}

// Generate a unique output directory name using current time and a random number.
string generateOutputDirName() {
    time_t now = time(0);
    tm *ltm = localtime(&now);
    int randomNum = rand() % 10000;
    return "MichelAnalysis_" + to_string(1900 + ltm->tm_year) +
           to_string(1 + ltm->tm_mon) + to_string(ltm->tm_mday) + "_" +
           to_string(ltm->tm_hour) + to_string(ltm->tm_min) + to_string(ltm->tm_sec) +
           "_" + to_string(randomNum);
}

// PMT calibration using LED data (triggerBits == 16).
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
    TH1F *histArea[N_PMTS];
    Long64_t nLEDFlashes[N_PMTS] = {0};
    for (int i = 0; i < N_PMTS; i++) {
        histArea[i] = new TH1F(Form("PMT%d_Area", i+1),
                               Form("PMT %d;ADC Counts;Events", i+1), 150, -50, 400);
    }
    Int_t triggerBits;
    Double_t area[23];
    calibTree->SetBranchAddress("triggerBits", &triggerBits);
    calibTree->SetBranchAddress("area", area);
    Long64_t nEntries = calibTree->GetEntries();
    cout << "Processing " << nEntries << " calibration events..." << endl;
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        calibTree->GetEntry(entry);
        if (triggerBits != 16)
            continue;
        for (int pmt = 0; pmt < N_PMTS; pmt++) {
            histArea[pmt]->Fill(area[PMT_CHANNEL_MAP[pmt]]);
            nLEDFlashes[pmt]++;
        }
    }
    for (int i = 0; i < N_PMTS; i++) {
        if (histArea[i]->GetEntries() < 1000) {
            cerr << "Warning: Insufficient data for PMT " << i+1 << endl;
            mu1[i] = 0;
            mu1_err[i] = 0;
            continue;
        }
        TF1 *fitFunc = new TF1("fitFunc", SPEfit, -50, 400, 8);
        Double_t histMean = histArea[i]->GetMean();
        Double_t histRMS = histArea[i]->GetRMS();
        fitFunc->SetParameters(1000, histMean - histRMS, histRMS/2,
                                1000, histMean, histRMS,
                                500, 200);
        histArea[i]->Fit(fitFunc, "Q0", "", -50, 400);
        mu1[i] = fitFunc->GetParameter(4);
        Double_t sigma_mu1 = fitFunc->GetParError(4);
        Double_t sigma1 = fitFunc->GetParameter(5);
        mu1_err[i] = sqrt(pow(sigma_mu1, 2) + pow(sigma1/sqrt(nLEDFlashes[i]), 2));
        delete fitFunc;
        delete histArea[i];
    }
    cout << "\nPMT Calibration Results (1PE peak positions):\n";
    cout << "PMT#  HardwareCh  mu1 [ADC]  Error [ADC]  N_events\n";
    cout << "--------------------------------------------------\n";
    for (int i = 0; i < N_PMTS; i++) {
        printf("PMT%02d     %2d       %6.2f      %5.2f      %6lld\n", 
               i+1, PMT_CHANNEL_MAP[i], mu1[i], mu1_err[i], nLEDFlashes[i]);
    }
    cout << endl;
    calibFile->Close();
}

// Main event selection function.
bool passesMainEventSelection(const Double_t *pulseH, const Double_t *baselineRMS,
                              const Double_t *area, const Int_t *peakPosition,
                              const Double_t *mu1) {
    int countAbove2PE = 0;
    for (int pmt = 0; pmt < N_PMTS; pmt++) {
        if (pulseH[PMT_CHANNEL_MAP[pmt]] > 2 * mu1[pmt])
            countAbove2PE++;
    }
    if (countAbove2PE >= 3) {
        vector<Double_t> peakPositions;
        for (int pmt = 0; pmt < N_PMTS; pmt++)
            peakPositions.push_back(peakPosition[PMT_CHANNEL_MAP[pmt]]);
        Double_t dummyMean, currentRMS;
        CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
        if (currentRMS < 2.5)
            return true;
    }
    else {
        int countConditionB = 0;
        for (int pmt = 0; pmt < N_PMTS; pmt++) {
            int ch = PMT_CHANNEL_MAP[pmt];
            if (pulseH[ch] > 3 * baselineRMS[ch] && (area[ch] / pulseH[ch]) > 1.2)
                countConditionB++;
        }
        if (countConditionB >= 3) {
            vector<Double_t> peakPositions;
            for (int pmt = 0; pmt < N_PMTS; pmt++)
                peakPositions.push_back(peakPosition[PMT_CHANNEL_MAP[pmt]]);
            Double_t dummyMean, currentRMS;
            CalculateMeanAndRMS(peakPositions, dummyMean, currentRMS);
            if (currentRMS < 2.5)
                return true;
        }
    }
    return false;
}

// A custom composite function to fit the Michel decay plus accidental background.
// par[0] : N0 (Michel amplitude)
// par[1] : τ_Michel (Michel lifetime)
// par[2] : A_accidental (accidental amplitude)
// tau_accidental : obtained from the accidental fit and fixed.
Double_t CombinedDecayAccidental(Double_t *x, Double_t *par) {
    // tau_accidental is defined as a global static variable after the accidental fit.
    extern double gTauAccidental;
    return par[0] * exp(-x[0]/par[1]) + par[2] * exp(-x[0]/gTauAccidental);
}
// Global variable to store tau_accidental for use in the fit function.
double gTauAccidental = 0;

// Main analysis function.
void analyzeMuonMichel(TChain *analysisChain, const Double_t *mu1, const string &outputDir) {
    gErrorIgnoreLevel = kError;
    Double_t area[23], pulseH[23], baselineRMS[23];
    Int_t triggerBits, peakPosition[23];
    Long64_t nsTime;
    analysisChain->SetBranchAddress("triggerBits", &triggerBits);
    analysisChain->SetBranchAddress("area", area);
    analysisChain->SetBranchAddress("pulseH", pulseH);
    analysisChain->SetBranchAddress("baselineRMS", baselineRMS);
    analysisChain->SetBranchAddress("peakPosition", peakPosition);
    analysisChain->SetBranchAddress("nsTime", &nsTime);
    mkdir(outputDir.c_str(), 0777);
    
    // Create histograms.
    TH1F *histDeltaT = new TH1F("DeltaT",
        "Muon-Michel Time Difference;Time to previous event(Muon)(#mus);Counts/0.1 #mus",
        100, 1.0, 10.0);
    TH1F *histMichelSpectrum = new TH1F("MichelSpectrum",
        "Michel Electron Energy;Photoelectrons;Events/10 p.e.", 100, 0, 1000);
    TH1F *histMuonPMTHits = new TH1F("MuonPMTHits",
        "PMT Multiplicity;Number of PMTs hit;Events", 12, 0, 12);
    TH1F *histSiPMMultiplicity = new TH1F("SiPMMultiplicity",
        "SiPM Multiplicity;Number of SiPMs hit;Events", 10, 0, 10);
    TH1F *histTriggerBits = new TH1F("TriggerBits",
        "Trigger Bits Distribution;Trigger Bits;Events", 64, 0, 64);
    
    // Adjust axis title offsets and margins for better display.
    histDeltaT->GetXaxis()->SetTitleOffset(1.1);
    histDeltaT->GetYaxis()->SetTitleOffset(1.4);
    histMichelSpectrum->GetXaxis()->SetTitleOffset(1.1);
    histMichelSpectrum->GetYaxis()->SetTitleOffset(1.4);
    histMuonPMTHits->GetXaxis()->SetTitleOffset(1.1);
    histMuonPMTHits->GetYaxis()->SetTitleOffset(1.4);
    histSiPMMultiplicity->GetXaxis()->SetTitleOffset(1.1);
    histSiPMMultiplicity->GetYaxis()->SetTitleOffset(1.4);
    histTriggerBits->GetXaxis()->SetTitleOffset(1.1);
    histTriggerBits->GetYaxis()->SetTitleOffset(1.4);
    
    int muonCount = 0, michelCount = 0, totalTrigger34Events = 0;
    // Beam event handling: record time of last beam event (10 µs = 10,000 ns)
    Long64_t lastBeamTime = -1000000;  
    Long64_t nEntries = analysisChain->GetEntries();
    cout << "Analyzing " << nEntries << " events..." << endl;
    
    // Create a histogram for the accidental time window (10-20 µs)
    TH1F *histAccidental = new TH1F("AccidentalDeltaT",
        "Accidental Time Spectrum;Time (#mus);Counts/0.1 #mus", 100, 10.0, 20.0);
    
    for (Long64_t entry = 0; entry < nEntries; entry++) {
        analysisChain->GetEntry(entry);
        
        // --- Beam Event Handling: record time of beam events (triggerBits==1)
        if (triggerBits == 1) {
            lastBeamTime = nsTime;
            continue; // Skip processing beam event.
        }
        // Reject any event occurring within 10 µs after a beam event.
        if ((nsTime - lastBeamTime) < 10000)
            continue;
        // ---------------------------
        
        histTriggerBits->Fill(triggerBits);
        if (triggerBits != 34 && triggerBits != 2)
            continue;
        totalTrigger34Events++;
        if (!passesMainEventSelection(pulseH, baselineRMS, area, peakPosition, mu1))
            continue;
        int sipmHitCount = 0;
        for (int i = 0; i < N_SIPMS; i++) {
            int ch = SIPM_CHANNEL_MAP[i];
            if (area[ch] >= SIPM_THRESHOLDS[i])
                sipmHitCount++;
        }
        histSiPMMultiplicity->Fill(sipmHitCount);
        if (sipmHitCount > 3)
            continue;
        int pmtHitCount = 0;
        for (int pmt = 0; pmt < N_PMTS; pmt++) {
            if (area[PMT_CHANNEL_MAP[pmt]] >= PMT_THRESHOLDS[pmt])
                pmtHitCount++;
        }
        histMuonPMTHits->Fill(pmtHitCount);
        muonCount++;
        Long64_t muonTime = nsTime;
        // Look for the Michel electron candidate in subsequent events.
        for (Long64_t nextEntry = entry + 1; nextEntry < nEntries; nextEntry++) {
            analysisChain->GetEntry(nextEntry);
            
            // --- Beam Window Rejection for Michel Candidate ---
            if (triggerBits == 1) {
                lastBeamTime = nsTime;
                continue;
            }
            if ((nsTime - lastBeamTime) < 10000)
                continue;
            // ----------------------------------------------------
            
            if (triggerBits != 34 && triggerBits != 2)
                continue;
            double deltaT = (nsTime - muonTime) * 1e-3; // convert ns to µs
            
            // Fill the accidental histogram if within 10-20 µs region
            if(deltaT >= 10 && deltaT <= 20)
                histAccidental->Fill(deltaT);
            
            if (deltaT > 10)
                break;
            if (deltaT < 1)
                continue;
            if (!passesMainEventSelection(pulseH, baselineRMS, area, peakPosition, mu1))
                continue;
            int michelPMTCount = 0;
            double michelEnergy = 0.0;
            for (int pmt = 0; pmt < N_PMTS; pmt++) {
                double pe = area[PMT_CHANNEL_MAP[pmt]] / mu1[pmt];
                if (pe >= 2) {
                    michelPMTCount++;
                    michelEnergy += pe;
                }
            }
            if (michelPMTCount >= 11) {
                michelCount++;
                histDeltaT->Fill(deltaT);
                histMichelSpectrum->Fill(michelEnergy);
                break;
            }
        }
    }
    cout << "\nANALYSIS SUMMARY:" << endl;
    cout << "Total events processed: " << nEntries << endl;
    cout << "TriggerBits==34 events: " << totalTrigger34Events << endl;
    cout << "Muon candidates: " << muonCount << endl;
    cout << "Michel electrons: " << michelCount << endl;
    cout << "Michel fraction: " << (muonCount > 0 ? 100.0 * michelCount / muonCount : 0) << "%" << endl;
    
    //===============================================
    // Fit the Accidental Region (10-20 µs) to obtain τ_accidental.
    //===============================================
    TF1 *accidentalFit = new TF1("accidentalFit", "[0]*exp(-x/[1])", 10, 20);
    accidentalFit->SetParameters(histAccidental->GetMaximum(), 5.0); // initial guess: A, τ_accidental ~5 µs
    histAccidental->Fit(accidentalFit, "RQ");
    gTauAccidental = accidentalFit->GetParameter(1);
    cout << Form("Accidental Fit: τ_accidental = %.4f µs", gTauAccidental) << endl;
    
    //----- Draw the accidental fit on a separate canvas -----//
    TCanvas *cAccidental = new TCanvas("cAccidental", "Accidental Background Fit", 800, 600);
    histAccidental->SetMarkerStyle(20);
    histAccidental->Draw("PE");
    accidentalFit->SetLineColor(kRed);
    accidentalFit->Draw("same");
    cAccidental->SaveAs((outputDir + "/accidental_fit.png").c_str());
    //---------------------------------------------------------//
    
    //===============================================
    // Multiple Exponential Fits on histDeltaT.
    //===============================================
    vector< pair<double,double> > fitRanges = {
         make_pair(1.0, 10.0),
         make_pair(1.5, 10.0),
         make_pair(2.0, 10.0),
         make_pair(2.5, 10.0)
    };
    vector<double> taus, tau_errors, chi2_ndf;
    vector<TF1*> fitFunctions;
    
    TCanvas *cFits = new TCanvas("cFits", "Time Difference Fits", 1200, 800);
    cFits->SetLeftMargin(0.15);
    cFits->SetRightMargin(0.10);
    cFits->SetBottomMargin(0.15);
    cFits->SetTopMargin(0.10);
    cFits->Divide(2,2);
    
    for (size_t i = 0; i < fitRanges.size(); i++) {
         double fitMin = fitRanges[i].first;
         double fitMax = fitRanges[i].second;
         TF1 *expFit = new TF1(Form("expFit_%zu", i), DecayFit, fitMin, fitMax, 2);
         expFit->SetParameters(histDeltaT->GetBinContent(histDeltaT->GetMaximumBin()), 2.2);
         expFit->SetParLimits(0, 0, histDeltaT->GetBinContent(histDeltaT->GetMaximumBin()) * 2);
         expFit->SetParLimits(1, 0.1, 10);
         histDeltaT->Fit(expFit, "RQ", "", fitMin, fitMax);
         
         double tau = expFit->GetParameter(1);
         double tauErr = expFit->GetParError(1);
         double chi2 = expFit->GetChisquare();
         double ndf = expFit->GetNDF();
         taus.push_back(tau);
         tau_errors.push_back(tauErr);
         chi2_ndf.push_back(chi2 / ndf);
         fitFunctions.push_back(expFit);
         
         cout << Form("Fit Range: %.1f - %.1f µs", fitMin, fitMax) << endl;
         cout << Form("#tau = %.4f ± %.4f µs", tau, tauErr) << endl;
         cout << Form("#chi^{2}/NDF = %.4f", chi2/ndf) << endl;
         cout << "----------------------------------------" << endl;
         
         cFits->cd(i+1);
         histDeltaT->SetMarkerStyle(20);
         histDeltaT->Draw("PE");
         expFit->Draw("same");
         TPaveText *pt = new TPaveText(0.60, 0.60, 0.80, 0.75, "NDC");
         pt->SetFillColor(0);
         pt->SetTextSize(0.035);
         pt->SetBorderSize(0);
         pt->AddText(Form("Fit Range: %.1f-%.1f #mus", fitMin, fitMax));
         pt->AddText(Form("#tau = %.4f #pm %.4f #mus", tau, tauErr));
         pt->AddText(Form("#chi^{2}/ndf = %.4f",chi2/ndf));
         pt->Draw("same");
    }
    cFits->SaveAs((outputDir + "/time_difference_fits.png").c_str());
    
    //===============================================
    // Composite Fit: Michel signal plus accidental background.
    //===============================================
    TF1 *compositeFit = new TF1("compositeFit", CombinedDecayAccidental, 2.5, 10.0, 3);
    compositeFit->SetParameters(histDeltaT->GetBinContent(histDeltaT->GetMaximumBin()),
                                2.2,
                                accidentalFit->GetParameter(0)); 
    histDeltaT->Fit(compositeFit, "RQ", "", 2.5, 10.0);
    
    double N0 = compositeFit->GetParameter(0);
    double tau_Michel = compositeFit->GetParameter(1);
    double A_acc = compositeFit->GetParameter(2);
    double chi2_ndf_comp = compositeFit->GetChisquare() / compositeFit->GetNDF();
    
    cout << "\nComposite Fit (2.5-10 µs):" << endl;
    cout << Form("  N0           = %.4f", N0) << endl;
    cout << Form("  τ_Michel     = %.4f µs", tau_Michel) << endl;
    cout << Form("  A_accidental = %.4f", A_acc) << endl;
    cout << Form("  τ_accidental fixed = %.4f µs", gTauAccidental) << endl;
    cout << Form("  χ²/ndf       = %.4f", chi2_ndf_comp) << endl;
    
    //===============================================
    // Overlayed Summary Plot (Same Y Axis)
    //===============================================
    TGraphErrors *gTau = new TGraphErrors();
    TGraph *gChi2 = new TGraph();
    int nPoints = taus.size();
    double commonYMin = taus[0], commonYMax = taus[0];
    for (int i = 0; i < nPoints; i++) {
         double x_value = fitRanges[i].first;
         gTau->SetPoint(i, x_value, taus[i]);
         gTau->SetPointError(i, 0, tau_errors[i]);
         gChi2->SetPoint(i, x_value, chi2_ndf[i]);
         if (taus[i] < commonYMin)
              commonYMin = taus[i];
         if (taus[i] > commonYMax)
              commonYMax = taus[i];
         if (chi2_ndf[i] < commonYMin)
              commonYMin = chi2_ndf[i];
         if (chi2_ndf[i] > commonYMax)
              commonYMax = chi2_ndf[i];
    }
    commonYMin *= 0.9;
    commonYMax *= 1.1;
    
    gTau->SetMarkerStyle(20);
    gTau->SetMarkerColor(kBlue);
    gTau->SetLineColor(kBlue);
    gChi2->SetMarkerStyle(21);
    gChi2->SetMarkerColor(kRed);
    gChi2->SetLineColor(kRed);
    
    TCanvas *cSummary = new TCanvas("cSummary", "Fit Summary Comparison", 800, 600);
    cSummary->SetLeftMargin(0.15);
    cSummary->SetRightMargin(0.10);
    cSummary->SetBottomMargin(0.15);
    cSummary->SetTopMargin(0.10);
    
    gTau->SetTitle("#tau and #chi^{2}/ndf vs. Fit Lower Bound;Fit Lower Bound (#mus);Value");
    gTau->Draw("ALP");
    gTau->GetYaxis()->SetRangeUser(commonYMin, commonYMax);
    gChi2->Draw("LP SAME");
    
    TLegend *leg = new TLegend(0.70, 0.75, 0.90, 0.90);
    leg->AddEntry(gTau, "#tau(#mus)", "lp");
    leg->AddEntry(gChi2, "#chi^{2}/ndf", "lp");
    leg->SetTextSize(0.035);
    leg->Draw();
    
    cSummary->SaveAs((outputDir + "/fit_comparison.png").c_str());
    
    //===============================================
    // Other Histogram Plots with Adjusted Margins and Offsets
    //===============================================
    auto plotHistogram = [&outputDir](TH1F* hist, const string& name,
                                      const string& title, bool addCutLine = false, double cutValue = 0) {
        TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1200, 800);
        c->SetLeftMargin(0.15);
        c->SetRightMargin(0.10);
        c->SetBottomMargin(0.15);
        c->SetTopMargin(0.10);
        hist->GetXaxis()->SetTitleOffset(1.1);
        hist->GetYaxis()->SetTitleOffset(1.3);
        hist->GetXaxis()->SetLabelSize(0.045);
        hist->GetYaxis()->SetLabelSize(0.045);
        // For PMT multiplicity plot, set axis properties to encourage scientific notation.
        if (name == "muon_pmt_hits") {
            hist->GetYaxis()->SetNoExponent(false);
            hist->GetYaxis()->SetMaxDigits(3);
            hist->GetYaxis()->SetNoExponent(kTRUE);
            // For compatibility with older ROOT versions
            hist->GetYaxis()->SetNoExponent(kFALSE);   // Let ROOT decide if scientific format is needed
            hist->GetYaxis()->SetTitleOffset(1.4);       // Move the axis title a bit if overlapping
        }
        hist->Draw();
        c->SaveAs((outputDir + "/" + name + ".png").c_str());
        delete c;
    };
    
    plotHistogram(histMichelSpectrum, "michel_spectrum", "Michel Electron Energy Spectrum");
    plotHistogram(histMuonPMTHits, "muon_pmt_hits", "PMT Multiplicity");
    plotHistogram(histSiPMMultiplicity, "sipm_multiplicity", "SiPM Multiplicity Distribution", true, 2.5);
    plotHistogram(histTriggerBits, "trigger_bits", "Trigger Bits Distribution");
    
    // Save histograms to a ROOT file.
    TFile *outFile = new TFile((outputDir + "/results.root").c_str(), "RECREATE");
    histDeltaT->Write();
    histMichelSpectrum->Write();
    histMuonPMTHits->Write();
    histSiPMMultiplicity->Write();
    histTriggerBits->Write();
    outFile->Close();
    
    // Cleanup.
    delete histDeltaT;
    delete histMichelSpectrum;
    delete histMuonPMTHits;
    delete histSiPMMultiplicity;
    delete histTriggerBits;
    delete histAccidental;
    delete cFits;
    delete cSummary;
    delete outFile;
}

// Main function.
int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file1.root> [input_file2.root ...]" << endl;
        cerr << "Note: First file is used for calibration (triggerBits==16), others for analysis" << endl;
        return 1;
    }
    string outputDir = generateOutputDirName();
    cout << "Creating output directory: " << outputDir << endl;
    mkdir(outputDir.c_str(), 0777);
    Double_t mu1[N_PMTS], mu1_err[N_PMTS];
    performCalibration(argv[1], mu1, mu1_err);
    TChain *analysisChain = new TChain("tree");
    for (int i = 1; i < argc; i++) {
        analysisChain->Add(argv[i]);
        cout << "Added file: " << argv[i] << endl;
    }
    analyzeMuonMichel(analysisChain, mu1, outputDir);
    delete analysisChain;
    cout << "Analysis complete. Results in: " << outputDir << endl;
    return 0;
}
