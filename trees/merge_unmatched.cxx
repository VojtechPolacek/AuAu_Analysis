#include <TFile.h>
#include <TH1D.h>
#include <TH2.h>
#include <TTree.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <iostream>
#include <TString.h>
#include <vector>
#include <TLegend.h>

std::vector<TH1D*> h_EEC_unmatched_vec;

void merge_unmatched(){
    h_EEC_unmatched_vec.clear(); // clear in case of re-running in same session

    const char* radii[]        = {"R0.2", "R0.3", "R0.4"};
    const char* types[]        = {"full", "charged"};
    const char* centralities[] = {"CENT_0_10", "MID_20_40", "PERI_60_80"};

    std::vector<TString> files;
    FILE* inlist = fopen("/gpfs/mnt/gpfs01/star/pwg/polacvo1/Analysis_AuAu_EEC/trees/inlist.txt", "r");
    if (!inlist) {
        std::cerr << "Failed to open inlist.txt" << std::endl;
        return;
    }

    char line[1024];
    while (fgets(line, sizeof(line), inlist)) {
        TString fname(line);
        fname.Remove(TString::kTrailing, '\n');
        fname.Remove(TString::kTrailing, '\r');
        if (fname.Length() > 0) files.push_back(fname);
    }
    fclose(inlist);
    std::cout << "Found " << files.size() << " files to merge." << std::endl;

    // Initialize vector with nullptr for all 18 combinations
    for (const auto& r : radii)
        for (const auto& t : types)
            for (const auto& c : centralities)
                h_EEC_unmatched_vec.push_back(nullptr);

    int nProcessed = 0;
    int ErrorCount = 0;
    int TotalHistograms = 0;

    for (const auto& fname : files) {
        TFile* f = TFile::Open(fname);
        if (!f || f->IsZombie()) {
            std::cerr << "Skipping bad file: " << fname << std::endl;
            continue;
        }

        int idx = 0;
        for (const auto& r : radii) {
            for (const auto& t : types) {
                for (const auto& c : centralities) {
                    TString path = Form("%s/%s/%s/hEEC_unmatched_all", r, t, c);
                    TH1D* h = (TH1D*)f->Get(path);

                    if (!h || h->GetEntries() == 0) {
                        idx++;
                        ErrorCount++;
                        TotalHistograms++;
                        continue;
                    }

                    if (h_EEC_unmatched_vec[idx] == nullptr) {
                        // unique name per combination to avoid ROOT name collisions
                        TString cloneName = Form("hEEC_unmatched_all_%s_%s_%s", r, t, c);
                        h_EEC_unmatched_vec[idx] = (TH1D*)h->Clone(cloneName);
                        h_EEC_unmatched_vec[idx]->SetDirectory(0);
                    } else {
                        h_EEC_unmatched_vec[idx]->Add(h);
                    }

                    idx++;
                    TotalHistograms++;
                }
            }
        }
        f->Close();
        nProcessed++;
        if (nProcessed % 100 == 0) {
            float ratio = (static_cast<float>(ErrorCount) / TotalHistograms) * 100;
            std::cout << "Processed " << nProcessed << "/" << files.size() << " files..." << std::endl;
            std::cout << "Current error count: " << ErrorCount << std::endl;
            std::cout << "Total histograms processed so far: " << TotalHistograms << std::endl;
            std::cout << "Current error ratio: " << ratio << " %" << std::endl;
            std::cout << "-----------------------------------------" << std::endl;
        }
    }

    // Write output file
    TFile* outFile = TFile::Open("/gpfs/mnt/gpfs01/star/pwg/polacvo1/Analysis_AuAu_EEC/trees/unmatched_merged.root", "RECREATE");

    int idx = 0;
    for (const auto& r : radii) {
        for (const auto& t : types) {
            for (const auto& c : centralities) {
                if (h_EEC_unmatched_vec[idx]) {
                    // Create directories level by level
                    if (!outFile->GetDirectory(r)) outFile->mkdir(r);
                    TDirectory* rdir = (TDirectory*)outFile->Get(r);

                    if (!rdir->GetDirectory(t)) rdir->mkdir(t);
                    TDirectory* tdir = (TDirectory*)rdir->Get(t);

                    if (!tdir->GetDirectory(c)) tdir->mkdir(c);
                    TDirectory* cdir = (TDirectory*)tdir->Get(c);

                    cdir->cd();
                    int bytes = h_EEC_unmatched_vec[idx]->Write();
                    if (bytes == 0)
                        std::cerr << "Warning: Write() returned 0 for " << r << "/" << t << "/" << c << std::endl;
                    else
                        std::cout << "Written " << bytes << " bytes for " << r << "/" << t << "/" << c << std::endl;

                    outFile->cd();
                } else {
                    std::cerr << "No valid histogram to write for " << r << "/" << t << "/" << c << std::endl;
                }
                idx++;
            }
        }
    }

    outFile->Close();
    std::cout << "Done. Processed " << nProcessed << " files -> unmatched_merged.root" << std::endl;
    std::cout << "Total histograms processed: " << TotalHistograms << std::endl;
    std::cout << "Total errors encountered: " << ErrorCount << std::endl;
}