// This code is to fill the histograms for PID purpose, 1st we need to pair the pion and then fill the histograms

#include "Iff2017.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TVector3.h>
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"
#include <fstream>
#include <iostream>
#include <TFile.h>

// Need to include histograms for ToF only, TPC only and ToF or TPC
// Need to find out the ToF and TPC yield

using namespace std;
void Iff2017::Loop()
{
    // gROOT->Reset();

    const int nBin = 5;
    const int nEtaBin = 9;
    const int ncharge = 2;
    const int nDir = 2;
    const char *etaDir[nDir] = {"Gt", "Lt"};
    const char *charge[ncharge] = {"Pos", "Neg"};

    double pT[6] = {2.60, 4.628, 5.643, 6.979, 9.265, 25};                                                      // pT bin range for A_UT Vs Minv
    double M[6] = {0.20, 0.4848, 0.6565, 0.8431, 1.1710, 4.000};                                                // Minv bin for A_UT Vs pT
    double eta_range[10] = {-1.200, -0.7243, -0.5539, -0.3410, -0.0225, 0.2833, 0.4852, 0.6346, 0.7695, 1.200}; // Eta bin for A_UT Vs Eta

    const char *det[3] = {"", "ToF", "TPCorToF"};

    const char *ToFCut[3] = {"TOF_Nocut", "TOF_Cut", "DontUse"};

    for (int j = 0; j < ncharge; j++)
    {
        for (int k = 0; k < nDir; k++)
        {
            for (int i = 0; i < nBin; i++)
            {
                h_dEdxVsp_pt[j][k][i] = new TH2D(Form("dEdx_Vs_p%s%s_ptBin%i", charge[j], etaDir[k], i), "dEdx_vs_p;log10(p)[GeV/c];log10(dEdx)[GeV/cm]", 100, -0.35, 1.55, 100, -0.120, 1.30);
                h_dEdxVsp_M[j][k][i] = new TH2D(Form("dEdx_Vs_p%s%s_MBin%i", charge[j], etaDir[k], i), "dEdx_vs_p;log10(p)[GeV/c];log10(dEdx)[GeV/cm]", 100, -0.35, 1.55, 100, -0.120, 1.30);
                for (int ndet = 0; ndet < 3; ndet++)
                {
                    // cout << "Hello There"<< endl;
                    hsigmaPionM[ndet][j][k][i] = new TH1D(Form("hsigmaPion%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hsigmaKaonM[ndet][j][k][i] = new TH1D(Form("hsigmaKaon%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hsigmaProtonM[ndet][j][k][i] = new TH1D(Form("hsigmaProton%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hsigmaElectronM[ndet][j][k][i] = new TH1D(Form("hsigmaElectron%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hmsqrM[ndet][j][k][i] = new TH1D(Form("hmsqr%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 2);

                    hsigmaPionPt[ndet][j][k][i] = new TH1D(Form("hsigmaPion%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hsigmaKaonPt[ndet][j][k][i] = new TH1D(Form("hsigmaKaon%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hsigmaProtonPt[ndet][j][k][i] = new TH1D(Form("hsigmaProton%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);
                    hsigmaElectronPt[ndet][j][k][i] = new TH1D(Form("hsigmaElectron%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -15, 15);

                    hmsqrPt[ndet][j][k][i] = new TH1D(Form("hmsqr%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 2);

                    hPivsKpt[ndet][j][k][i] = new TH2D(Form("hPivsK%s%s_ptBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                    hPivsPpt[ndet][j][k][i] = new TH2D(Form("hPivsP%s%s_ptBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                    hPivsEpt[ndet][j][k][i] = new TH2D(Form("hPivsE%s%s_ptBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -3, 10);

                    hPivsKm[ndet][j][k][i] = new TH2D(Form("hPivsK%s%s_mBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                    hPivsPm[ndet][j][k][i] = new TH2D(Form("hPivsP%s%s_mBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                    hPivsEm[ndet][j][k][i] = new TH2D(Form("hPivsE%s%s_mBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -3, 10);

                    hPvsKpt[ndet][j][k][i] = new TH2D(Form("hPvsK%s%s_ptBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                    hPvsKm[ndet][j][k][i] = new TH2D(Form("hPvsK%s%s_mBin%i%s", charge[j], etaDir[k], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);

                    hpM[ndet][j][k][i] = new TH1D(Form("hp%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), " ", 100, 0, 25);
                    hptM[ndet][j][k][i] = new TH1D(Form("hpt%s_%s_Mbin%i%s", charge[j], etaDir[k], i, det[ndet]), " ", 100, 0, 25);
                    hpP[ndet][j][k][i] = new TH1D(Form("hp%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), " ", 100, 0, 25);
                    hptP[ndet][j][k][i] = new TH1D(Form("hpt%s_%s_Ptbin%i%s", charge[j], etaDir[k], i, det[ndet]), " ", 100, 0, 25);

                    hsigmaPion_TpcVsToF_pt[ndet][j][k][i] = new TH2D(Form("hsigmaPion_TpcVsToF%s%s_ptBin%i%s", charge[j], etaDir[k], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
                    hsigmaKaon_TpcVsToF_pt[ndet][j][k][i] = new TH2D(Form("hsigmaKaon_TpcVsToF%s%s_ptBin%i%s", charge[j], etaDir[k], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
                    hsigmaProton_TpcVsToF_pt[ndet][j][k][i] = new TH2D(Form("hsigmaProton_TpcVsToF%s%s_ptBin%i%s", charge[j], etaDir[k], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);

                    hsigmaPion_TpcVsToF_M[ndet][j][k][i] = new TH2D(Form("hsigmaPion_TpcVsToF%s%s_MBin%i%s", charge[j], etaDir[k], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
                    hsigmaKaon_TpcVsToF_M[ndet][j][k][i] = new TH2D(Form("hsigmaKaon_TpcVsToF%s%s_MBin%i%s", charge[j], etaDir[k], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
                    hsigmaProton_TpcVsToF_M[ndet][j][k][i] = new TH2D(Form("hsigmaProton_TpcVsToF%s%s_MBin%i%s", charge[j], etaDir[k], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);

                } // nBin loop
            }     // nDir loop
        }         // ncharge loop
    }             // TPC, TOF, TPC or ToF loop
    for (int ndet = 0; ndet < 3; ndet++)
    {
        for (int j = 0; j < ncharge; j++)
        {
            for (int i = 0; i < nEtaBin; i++)
            {
                hsigmaPionEta[ndet][j][i] = new TH1D(Form("hsigmaPion%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -15, 15);
                hsigmaKaonEta[ndet][j][i] = new TH1D(Form("hsigmaKaon%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -15, 15);
                hsigmaElectronEta[ndet][j][i] = new TH1D(Form("hsigmaElectron%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -15, 15);
                hsigmaProtonEta[ndet][j][i] = new TH1D(Form("hsigmaProton%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -15, 15);
                hmsqrEta[ndet][j][i] = new TH1D(Form("hmsqr%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -1, 2);

                hPivsKeta[ndet][j][i] = new TH2D(Form("hPivsK%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                hPivsPeta[ndet][j][i] = new TH2D(Form("hPivsP%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);
                hPivsEeta[ndet][j][i] = new TH2D(Form("hPivsE%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -1, 1, 100, -3, 10);
                hPvsKeta[ndet][j][i] = new TH2D(Form("hPvsK%s_etaBin%i%s", charge[j], i, det[ndet]), "", 100, -1, 1, 100, -5, 5);

                hsigmaPion_TpcVsToF_eta[ndet][j][i] = new TH2D(Form("hsigmaPion_TpcVsToF%s_etaBin%i%s", charge[j], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
                hsigmaKaon_TpcVsToF_eta[ndet][j][i] = new TH2D(Form("hsigmaKaon_TpcVsToF%s_etaBin%i%s", charge[j], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
                hsigmaProton_TpcVsToF_eta[ndet][j][i] = new TH2D(Form("hsigmaProton_TpcVsToF%s_etaBin%i%s", charge[j], i, ToFCut[ndet]), "", 100, -5, 10, 100, -5, 10);
            }
        }
    }
    cout << "There is No error Eta Loop" << endl;
    h_dEdx_p = new TH2D("dEdx_Vs_p", "dEdx_Vs_p; log10(p)[GeV/c]; log10(dEdx)[GeV/cm]", 100, -0.35, 1.55, 100, -0.120, 1.30);
    hnSigmaPionTPC = new TH1D("nsigmaPionTPC", "nsigmaPionTPC", 100, -5, 100);
    hnSigmaPionToF = new TH1D("nsigmaPionToF", "nsigmaPionToF", 100, -5, 100);

    double pT_pos, pT_neg;
    double p1x, p2x, p1y, p2y, p1z, p2z, psx, psy, psz, ps, cone, R, Rx, Ry, Rz, R1, Rx1, Ry1, Rz1, p_tr, p_tr2, ptr, ptr2;
    double p1, p2, E1, E2, Minv, pT_pair, pT_min_pair, eta_pair, fitPts_min_pair;
    double cosPhiRB, sinPhiRB, cosPhiSB, sinPhiSB, cosPhiRY, sinPhiRY, cosPhiSY, sinPhiSY;
    double PhiR_cosB, PhiR_sinB, PhiRB, PhiS_cosB, PhiS_sinB, PhiSB, PhiRSB, PhiR_cosY, PhiR_sinY, PhiRY, PhiS_cosY, PhiS_sinY, PhiSY, PhiRSY, phiSB, phiRB, phiDiffB, phiRSB;
    int spin51 = 0, spin53 = 0, spin83 = 0, spin85 = 0, blueUp = 0, blueDown = 0, yellowUp = 0, yellowDown = 0;
    double msqr1, msqr2;
    double pi = 3.14159265359;
    //  double CONE_CUT_MIN = 0.05;
    double m_pion = 0.1396; // GeV
    int fillnum;
    unsigned int evTime;

    int paircount = 0;
    int ntrack1 = 0, ntrack2 = 0;
    if (fChain == 0)
        return;
    Long64_t nentries = fChain->GetEntries();
    cout << "Muji" << endl;
    cout << "Entries:\t " << nentries << endl;

    // Event Loop

    cout << fChain->GetEntries() << " \t Event Num" << endl;
    for (Long64_t jentry = 0; jentry < nentries; jentry++)
    {
        // cout << jentry << " \t Event Num" << endl;
        // for(Long64_t jentry=0;jentry<50;jentry++){

        fChain->GetEntry(jentry);
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        // cout << fspinconfig << "\t spin config\t "<< endl;
        // cout << fVZ << "\t fVZ config\t "<< endl;

        if (fabs(fVZ) > 90)
            continue;
        if (fspinconfig != 51 && fspinconfig != 53 && fspinconfig != 83 && fspinconfig != 85)
            continue;

        // Trigger Selection
        int trig_JP0VPDMB30 = 0;
        int trig_JP1VPDMB30 = 0;
        int trig_JP0VPDMBLt30 = 0;
        int trig_JP0VPDMBGt30 = 0;
        int trig_JP1VPDMBLt30 = 0;
        int trig_JP1VPDMBGt30 = 0;
        int trig_JP2 = 0;
        int trig_BHT3 = 0;
        int trig_BHT1VPDMB30 = 0;
        int trig_BHT2BBCMB = 0;
        int trig_JPsiHTTP = 0;
        int trig_VPDMB30 = 0;

        for (int trig = 0; trig < ftrigger->size(); trig++)
        {

            if (ftrigger->at(trig) == 570401)
            {
                trig_JP2 = 1;
            } // JP2 trigger Flag

            if (ftrigger->at(trig) == 570403)
            {
                trig_JP1VPDMB30 = 1;
            } // JP1 Trigger

            if (ftrigger->at(trig) == 570404)
            {
                trig_JP0VPDMB30 = 1;
            } // JP0 Trigger

            if (ftrigger->at(trig) == 570201)
            {
                trig_BHT3 = 1;
            } // BHT3 trigger Flag
            if (ftrigger->at(trig) == 570215)
            {
                trig_BHT2BBCMB = 1;
            } // BHT2 Trigger

            if (ftrigger->at(trig) == 570214)
            {
                trig_BHT1VPDMB30 = 1;
            } // BHT1 Trigger

            if (ftrigger->at(trig) == 570229 || ftrigger->at(trig) == 570219 || ftrigger->at(trig) == 570209)
            {
                trig_JPsiHTTP = 1;
            } // JPsiHTTP trigger
            if (ftrigger->at(trig) == 570001)
            {
                trig_VPDMB30 = 1;
            } // VPDMB30 Trigger

        } // frigger Loop

        if (trig_JP2 != 1 && trig_JP0VPDMB30 != 1 && trig_JP1VPDMB30 != 1 && trig_BHT3 != 1 && trig_BHT1VPDMB30 != 1 && trig_BHT2BBCMB != 1 && trig_JPsiHTTP != 1 && trig_VPDMB30 != 1)
            continue;

        double msqrPlus = 0;
        double msqrMinus = 0;
        // Track Loop
        for (int tr = 0; tr < fmaxpar; tr++)
        {
            if (ffitPtsPoss[tr] <= 0)
                continue;
            ntrack1++;
            double fitPtsRatio_tr = static_cast<double>(ffitPts[tr]) / (static_cast<double>(ffitPtsPoss[tr]));

            if (fp[tr] > 0.0 && fdca[tr] < 1. && ffitPts[tr] > 15 && fhitsdedx[tr] > 20 && feta[tr] <= 1. && feta[tr] >= -1. && fitPtsRatio_tr > .51)

            {
                double log_fp = log10(fp[tr]);
                double log_fdEdx = log10(fdEdx[tr]);
                // h_dEdx_p->Fill(fp[tr], fdEdx[tr]);
                h_dEdx_p->Fill(log_fp, log_fdEdx);
            }

            if (fp[tr] > 2 && fnSigmaPion[tr] >= -15. && fnSigmaPion[tr] <= 15. && fdca[tr] < 1. && ffitPts[tr] > 15 && fhitsdedx[tr] > 20 && feta[tr] <= 1. && feta[tr] >= -1. && fitPtsRatio_tr > .51)
            {

                double inverse_pion_beta1 = sqrt(pow(m_pion, 2) / pow(fp[tr], 2) + 1);
                if (fBetaToF[tr] != -999)
                    msqr1 = fp[tr] * fp[tr] * (1 / (fBetaToF[tr] * fBetaToF[tr]) - 1);

                // dummy histograms just to check ToF yeild
                hnSigmaPionTPC->Fill(fnSigmaPion[tr]);
                if (fnSigmaPionToF[tr] != -999)
                {
                    hnSigmaPionToF->Fill(fnSigmaPionToF[tr]);
                }
                for (int tr2 = tr + 1; tr2 < fmaxpar; tr2++)
                {
                    if (ffitPtsPoss[tr2] <= 0)
                        continue;
                    ntrack2++;
                    double fitPtsRatio_tr2 = (static_cast<double>(ffitPts[tr2])) / (static_cast<double>(ffitPtsPoss[tr2]));
                    if (fp[tr2] > 2 && fnSigmaPion[tr2] >= -15. && fnSigmaPion[tr2] <= 15. && fdca[tr2] < 1. && ffitPts[tr2] > 15 && fhitsdedx[tr2] > 20 && feta[tr2] <= 1. && feta[tr2] >= -1. && fitPtsRatio_tr2 > .51)
                    {

                        double inverse_pion_beta2 = sqrt(pow(m_pion, 2) / pow(fp[tr2], 2) + 1);
                        if (fBetaToF[tr2] != -999)
                            msqr2 = fp[tr2] * fp[tr2] * (1 / (fBetaToF[tr2] * fBetaToF[tr2]) - 1);

                        if (fcharge[tr] != fcharge[tr2])
                        {
                            double phiDiff = fphi[tr] - fphi[tr2];
                            if (phiDiff > pi)
                                phiDiff -= (2 * pi);
                            if (phiDiff < ((-1) * pi))
                                phiDiff += (2 * pi);
                            if (phiDiff > pi || phiDiff < -pi)
                                continue;
                            cone = sqrt(pow(feta[tr] - feta[tr2], 2) + pow(phiDiff, 2)); // cone cut < 0.7
                            if (cone <= 0.7)
                            {
                                paircount++;
                                if (fcharge[tr] > 0)
                                {
                                    p1x = fpT[tr] * cos(fphi[tr]);
                                    p2x = fpT[tr2] * cos(fphi[tr2]);
                                    p1y = fpT[tr] * sin(fphi[tr]);
                                    p2y = fpT[tr2] * sin(fphi[tr2]);
                                    p1z = fpT[tr] * sinh(feta[tr]);
                                    // pz= pcos(theta),pT=psin(theta), pz=pTcot(theta),sinh(eta)=((e^(eta)-e^(-eta))/2)=cot(theta);
                                    p2z = fpT[tr2] * sinh(feta[tr2]);
                                    p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
                                    p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);

                                    p_tr = fpT[tr] * cosh(feta[tr]);
                                    p_tr2 = fpT[tr2] * cosh(feta[tr2]);

                                    ptr = fp[tr];
                                    ptr2 = fp[tr2];
                                }
                                if (fcharge[tr] < 0)
                                {
                                    p1x = fpT[tr2] * cos(fphi[tr2]);
                                    p2x = fpT[tr] * cos(fphi[tr]);
                                    p1y = fpT[tr2] * sin(fphi[tr2]);
                                    p2y = fpT[tr] * sin(fphi[tr]);
                                    p1z = fpT[tr2] * sinh(feta[tr2]);
                                    p2z = fpT[tr] * sinh(feta[tr]);
                                    p1 = sqrt(p1x * p1x + p1y * p1y + p1z * p1z);
                                    p2 = sqrt(p2x * p2x + p2y * p2y + p2z * p2z);
                                    // cross check momentum ....
                                    p_tr = fpT[tr2] * cosh(feta[tr2]);
                                    p_tr2 = fpT[tr] * cosh(feta[tr]);

                                    ptr = fp[tr2];
                                    ptr2 = fp[tr];
                                }
                                // cout << "p_tr: "<< p_tr <<",  "<< p1<<", "<<ptr<< endl;
                                // cout << "p_tr2: "<< p_tr2 <<",  "<< p2<<", "<<ptr2<< endl;
                                // momentum all good !!
                                if (fpT[tr] > fpT[tr2])
                                    pT_min_pair = fpT[tr2];
                                if (fpT[tr] < fpT[tr2])
                                    pT_min_pair = fpT[tr];
                                if (ffitPts[tr] > ffitPts[tr2])
                                    fitPts_min_pair = ffitPts[tr2];
                                if (ffitPts[tr] < ffitPts[tr2])
                                    fitPts_min_pair = ffitPts[tr];

                                // Components of sum vector
                                psx = p1x + p2x;
                                psy = p1y + p2y;
                                psz = p1z + p2z;
                                // Sum Vector
                                ps = sqrt(psx * psx + psy * psy + psz * psz);
                                // Relatve momentum of dihadron system
                                // Charge ordering is important. In R = (Ph1-Ph2)*.5 ,I want  Ph1 to be positive and Ph2 to be negative always.
                                double Rx1, Rz1, Ry1;
                                if (fcharge[tr] > 0)
                                {
                                    Rx1 = (p1x - p2x);
                                    Ry1 = (p1y - p2y);
                                    Rz1 = (p1z - p2z);
                                }
                                else if (fcharge[tr2] > 0)
                                {
                                    Rx1 = (p2x - p1x);
                                    Ry1 = (p2y - p1y);
                                    Rz1 = (p2z - p1z);
                                }
                                Rx = (p1x - p2x);
                                Ry = (p1y - p2y);
                                Rz = (p1z - p2z);
                                R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz); // R and R1 are same
                                R1 = sqrt(Rx1 * Rx1 + Ry1 * Ry1 + Rz1 * Rz1);
                                // cout << "R1: "<< R1 << ", R: "<< R << endl; //Exact same output
                                // calculate M,pt,eta
                                E1 = sqrt(m_pion * m_pion + p1 * p1);
                                E2 = sqrt(m_pion * m_pion + p2 * p2);
                                Minv = sqrt(2 * m_pion * m_pion + 2 * (E1 * E2 - p1x * p2x - p1y * p2y - p1z * p2z));
                                pT_pair = sqrt((p1x + p2x) * (p1x + p2x) + (p1y + p2y) * (p1y + p2y)); // pair is vector sum of two track and pT_pair is the transverse part of the resultant ve    ctor of two track momentum vectors
                                eta_pair = TMath::ASinH((p1z + p2z) / pT_pair);

                                double sigmaPionP = 0;
                                double sigmaKaonP = 0;
                                double sigmaProtonP = 0;
                                double sigmaElectronP = 0;
                                double sigmaPionN = 0;
                                double sigmaProtonN = 0;
                                double sigmaElectronN = 0;
                                double sigmaKaonN = 0;

                                double sigmaPionToFP = 0;
                                double sigmaKaonToFP = 0;
                                double sigmaProtonToFP = 0;
                                double sigmaPionToFN = 0;
                                double sigmaProtonToFN = 0;
                                double sigmaKaonToFN = 0;

                                double pP = 0;
                                double ptP = 0;
                                double pN = 0;
                                double ptN = 0;
                                double log_dEdxP = 0;
                                double log_dEdxN = 0;
                                double log_pP = 0;
                                double log_pN = 0;

                                if (fcharge[tr] > 0 || fcharge[tr2] > 0)
                                {
                                    if (fcharge[tr] > 0)
                                    {
                                        pP = fp[tr];
                                        ptP = fpT[tr];
                                        sigmaPionP = fnSigmaPion[tr];
                                        sigmaKaonP = fnSigmaKaon[tr];
                                        sigmaProtonP = fnSigmaProton[tr];
                                        sigmaElectronP = fnSigmaElectron[tr];

                                        sigmaPionToFP = fnSigmaPionToF[tr];
                                        sigmaKaonToFP = fnSigmaKaonToF[tr];
                                        sigmaProtonToFP = fnSigmaProtonToF[tr];

                                        msqrPlus = msqr1;
                                        log_dEdxP = log10(fdEdx[tr]);
                                        log_pP = log10(fp[tr]);
                                    }
                                    if (fcharge[tr2] > 0)
                                    {
                                        pP = fp[tr2];
                                        ptP = fpT[tr2];
                                        sigmaPionP = fnSigmaPion[tr2];
                                        sigmaKaonP = fnSigmaKaon[tr2];
                                        sigmaProtonP = fnSigmaProton[tr2];
                                        sigmaElectronP = fnSigmaElectron[tr2];

                                        sigmaPionToFP = fnSigmaPionToF[tr2];
                                        sigmaKaonToFP = fnSigmaKaonToF[tr2];
                                        sigmaProtonToFP = fnSigmaProtonToF[tr2];

                                        msqrPlus = msqr2;
                                        log_dEdxP = log10(fdEdx[tr2]);
                                        log_pP = log10(fp[tr2]);
                                    }
                                    // cout << sigmaPionP<< "\t sigmaPionP\t "<< endl;

                                    for (int l = 0; l < 5; l++)
                                    {
                                        if (pT_pair >= pT[l] && pT_pair < pT[l + 1] && eta_pair > 0)
                                        {

                                            h_dEdxVsp_pt[0][0][l]->Fill(log_pP, log_dEdxP);
                                            // Selecting only Tracks that have ToF info
                                            // if (sigmaPionToFP != -999)
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            {
                                                hsigmaPion_TpcVsToF_pt[0][0][0][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_pt[0][0][0][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_pt[0][0][0][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // ToF_noCut
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            {
                                                hsigmaPion_TpcVsToF_pt[1][0][0][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_pt[1][0][0][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_pt[1][0][0][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // ToF_Cut

                                            // TPC
                                            hmsqrPt[0][0][0][l]->Fill(msqrPlus);
                                            hsigmaPionPt[0][0][0][l]->Fill(sigmaPionP);
                                            hsigmaProtonPt[0][0][0][l]->Fill(sigmaProtonP);
                                            hsigmaKaonPt[0][0][0][l]->Fill(sigmaKaonP);
                                            hsigmaElectronPt[0][0][0][l]->Fill(sigmaElectronP);
                                            hPivsKpt[0][0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                            hPivsPpt[0][0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                            hPivsEpt[0][0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                            hPvsKpt[0][0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                            hpP[0][0][0][l]->Fill(pP);
                                            hptP[0][0][0][l]->Fill(ptP);

                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2) || sigmaPionToFP == -999)
                                            {
                                                // TPC or ToF
                                                hmsqrPt[2][0][0][l]->Fill(msqrPlus);
                                                hsigmaPionPt[2][0][0][l]->Fill(sigmaPionP);
                                                hsigmaProtonPt[2][0][0][l]->Fill(sigmaProtonP);
                                                hsigmaKaonPt[2][0][0][l]->Fill(sigmaKaonP);
                                                hsigmaElectronPt[2][0][0][l]->Fill(sigmaElectronP);
                                                hPivsKpt[2][0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPpt[2][0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEpt[2][0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKpt[2][0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpP[2][0][0][l]->Fill(pP);
                                                hptP[2][0][0][l]->Fill(ptP);
                                            } // TPC or ToF control
                                            // ToF only
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                hmsqrPt[1][0][0][l]->Fill(msqrPlus);
                                                hsigmaPionPt[1][0][0][l]->Fill(sigmaPionP);
                                                hsigmaProtonPt[1][0][0][l]->Fill(sigmaProtonP);
                                                hsigmaKaonPt[1][0][0][l]->Fill(sigmaKaonP);
                                                hsigmaElectronPt[1][0][0][l]->Fill(sigmaElectronP);
                                                hPivsKpt[1][0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPpt[1][0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEpt[1][0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKpt[1][0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpP[1][0][0][l]->Fill(pP);
                                                hptP[1][0][0][l]->Fill(ptP);
                                            } // ToF Control

                                        } // contol statement for pT bining and eta_pair>0
                                        if (Minv >= M[l] && Minv < M[l + 1] && eta_pair > 0)
                                        {
                                            h_dEdxVsp_M[0][0][l]->Fill(log_pP, log_dEdxP);
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            // if (sigmaPionToFP != -999)
                                            {
                                                hsigmaPion_TpcVsToF_M[0][0][0][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_M[0][0][0][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_M[0][0][0][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // TOF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_M[1][0][0][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_M[1][0][0][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_M[1][0][0][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // ToF_Cut

                                            hmsqrM[0][0][0][l]->Fill(msqrPlus);
                                            hsigmaPionM[0][0][0][l]->Fill(sigmaPionP);
                                            hsigmaProtonM[0][0][0][l]->Fill(sigmaProtonP);
                                            hsigmaKaonM[0][0][0][l]->Fill(sigmaKaonP);
                                            hsigmaElectronM[0][0][0][l]->Fill(sigmaElectronP);
                                            hPivsKm[0][0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                            hPivsPm[0][0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                            hPivsEm[0][0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                            hPvsKm[0][0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                            hpM[0][0][0][l]->Fill(pP);
                                            hptM[0][0][0][l]->Fill(ptP);
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2) || sigmaPionToFP == -999)
                                            {
                                                // TPC or TOF
                                                hmsqrM[2][0][0][l]->Fill(msqrPlus);
                                                hsigmaPionM[2][0][0][l]->Fill(sigmaPionP);
                                                hsigmaProtonM[2][0][0][l]->Fill(sigmaProtonP);
                                                hsigmaKaonM[2][0][0][l]->Fill(sigmaKaonP);
                                                hsigmaElectronM[2][0][0][l]->Fill(sigmaElectronP);
                                                hPivsKm[2][0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPm[2][0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEm[2][0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKm[2][0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpM[2][0][0][l]->Fill(pP);
                                                hptM[2][0][0][l]->Fill(ptP);
                                            } // TPC or ToF
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                // TOF only
                                                hmsqrM[1][0][0][l]->Fill(msqrPlus);
                                                hsigmaPionM[1][0][0][l]->Fill(sigmaPionP);
                                                hsigmaProtonM[1][0][0][l]->Fill(sigmaProtonP);
                                                hsigmaKaonM[1][0][0][l]->Fill(sigmaKaonP);
                                                hsigmaElectronM[1][0][0][l]->Fill(sigmaElectronP);
                                                hPivsKm[1][0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPm[1][0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEm[1][0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKm[1][0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpM[1][0][0][l]->Fill(pP);
                                                hptM[1][0][0][l]->Fill(ptP);
                                            } // ToF only

                                        } // control statement for Minv bining and eta_pair>0

                                        if (pT_pair >= pT[l] && pT_pair < pT[l + 1] && eta_pair < 0)
                                        {
                                            h_dEdxVsp_pt[0][1][l]->Fill(log_pP, log_dEdxP);

                                            // if (sigmaPionToFP != -999)
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            {
                                                hsigmaPion_TpcVsToF_pt[0][0][1][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_pt[0][0][1][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_pt[0][0][1][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // ToF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {

                                                hsigmaPion_TpcVsToF_pt[1][0][1][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_pt[1][0][1][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_pt[1][0][1][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // ToF_Cut

                                            // TPC
                                            hmsqrPt[0][0][1][l]->Fill(msqrPlus);
                                            hsigmaPionPt[0][0][1][l]->Fill(sigmaPionP);
                                            hsigmaProtonPt[0][0][1][l]->Fill(sigmaProtonP);
                                            hsigmaKaonPt[0][0][1][l]->Fill(sigmaKaonP);
                                            hsigmaElectronPt[0][0][1][l]->Fill(sigmaElectronP);
                                            hPivsKpt[0][0][1][l]->Fill(sigmaKaonP, sigmaPionP);
                                            hPivsPpt[0][0][1][l]->Fill(sigmaProtonP, sigmaPionP);
                                            hPivsEpt[0][0][1][l]->Fill(sigmaElectronP, sigmaPionP);
                                            hPvsKpt[0][0][1][l]->Fill(sigmaKaonP, sigmaProtonP);

                                            hpP[0][0][1][l]->Fill(pP);
                                            hptP[0][0][1][l]->Fill(ptP);

                                            // if ((sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2) || sigmaPionToFP == -999)
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            {
                                                // TPC or TOF
                                                hmsqrPt[2][0][1][l]->Fill(msqrPlus);
                                                hsigmaPionPt[2][0][1][l]->Fill(sigmaPionP);
                                                hsigmaProtonPt[2][0][1][l]->Fill(sigmaProtonP);
                                                hsigmaKaonPt[2][0][1][l]->Fill(sigmaKaonP);
                                                hsigmaElectronPt[2][0][1][l]->Fill(sigmaElectronP);
                                                hPivsKpt[2][0][1][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPpt[2][0][1][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEpt[2][0][1][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKpt[2][0][1][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpP[2][0][1][l]->Fill(pP);
                                                hptP[2][0][1][l]->Fill(ptP);
                                            } // TPC or TOF
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //  if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                // TOF only
                                                hmsqrPt[1][0][1][l]->Fill(msqrPlus);
                                                hsigmaPionPt[1][0][1][l]->Fill(sigmaPionP);
                                                hsigmaProtonPt[1][0][1][l]->Fill(sigmaProtonP);
                                                hsigmaKaonPt[1][0][1][l]->Fill(sigmaKaonP);
                                                hsigmaElectronPt[1][0][1][l]->Fill(sigmaElectronP);
                                                hPivsKpt[1][0][1][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPpt[1][0][1][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEpt[1][0][1][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKpt[1][0][1][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpP[1][0][1][l]->Fill(pP);
                                                hptP[1][0][1][l]->Fill(ptP);
                                            } // ToF only

                                        } // control statement for pT binning with eta_pair<0

                                        if (Minv >= M[l] && Minv < M[l + 1] && eta_pair < 0)
                                        {

                                            h_dEdxVsp_M[0][1][l]->Fill(log_pP, log_dEdxP);
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            // if (sigmaPionToFP != -999)
                                            {
                                                hsigmaPion_TpcVsToF_M[0][0][1][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_M[0][0][1][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_M[0][0][1][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // TOF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_M[1][0][1][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_M[1][0][1][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_M[1][0][1][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // TOF_Cut

                                            // TPC
                                            hmsqrM[0][0][1][l]->Fill(msqrPlus);
                                            hsigmaPionM[0][0][1][l]->Fill(sigmaPionP);
                                            hsigmaProtonM[0][0][1][l]->Fill(sigmaProtonP);
                                            hsigmaKaonM[0][0][1][l]->Fill(sigmaKaonP);
                                            hsigmaElectronM[0][0][1][l]->Fill(sigmaElectronP);
                                            hPivsKm[0][0][1][l]->Fill(sigmaKaonP, sigmaPionP);
                                            hPivsPm[0][0][1][l]->Fill(sigmaProtonP, sigmaPionP);
                                            hPivsEm[0][0][1][l]->Fill(sigmaElectronP, sigmaPionP);
                                            hPvsKm[0][0][1][l]->Fill(sigmaKaonP, sigmaProtonP);
                                            hpM[0][0][1][l]->Fill(pP);
                                            hptM[0][0][1][l]->Fill(ptP);

                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2) || sigmaPionToFP == -999)
                                            {
                                                hmsqrM[2][0][1][l]->Fill(msqrPlus);
                                                hsigmaPionM[2][0][1][l]->Fill(sigmaPionP);
                                                hsigmaProtonM[2][0][1][l]->Fill(sigmaProtonP);
                                                hsigmaKaonM[2][0][1][l]->Fill(sigmaKaonP);
                                                hsigmaElectronM[2][0][1][l]->Fill(sigmaElectronP);
                                                hPivsKm[2][0][1][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPm[2][0][1][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEm[2][0][1][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKm[2][0][1][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpM[2][0][1][l]->Fill(pP);
                                                hptM[2][0][1][l]->Fill(ptP);
                                            } // TPC or TOF
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //  if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                hmsqrM[1][0][1][l]->Fill(msqrPlus);
                                                hsigmaPionM[1][0][1][l]->Fill(sigmaPionP);
                                                hsigmaProtonM[1][0][1][l]->Fill(sigmaProtonP);
                                                hsigmaKaonM[1][0][1][l]->Fill(sigmaKaonP);
                                                hsigmaElectronM[1][0][1][l]->Fill(sigmaElectronP);
                                                hPivsKm[1][0][1][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPm[1][0][1][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEm[1][0][1][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKm[1][0][1][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hpM[1][0][1][l]->Fill(pP);
                                                hptM[1][0][1][l]->Fill(ptP);
                                            } // ToF
                                        }     // contorl statment for Minv bining with eta_pair<0

                                    } // 5 pT or Minv Bin loop

                                    for (int l = 0; l < 9; l++)
                                    {
                                        if (eta_pair >= eta_range[l] && eta_pair < eta_range[l + 1])
                                        {
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            // if (sigmaPionToFP != -999)
                                            {
                                                hsigmaPion_TpcVsToF_eta[0][0][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_eta[0][0][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_eta[0][0][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // TOF_Nocut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_eta[1][0][l]->Fill(sigmaPionToFP, sigmaPionP);
                                                hsigmaKaon_TpcVsToF_eta[1][0][l]->Fill(sigmaKaonToFP, sigmaKaonP);
                                                hsigmaProton_TpcVsToF_eta[1][0][l]->Fill(sigmaProtonToFP, sigmaProtonP);
                                            } // ToF_Cut

                                            // TPC
                                            hmsqrEta[0][0][l]->Fill(msqrPlus);
                                            hPivsKeta[0][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                            hPivsPeta[0][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                            hPivsEeta[0][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                            hPvsKeta[0][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                            hsigmaPionEta[0][0][l]->Fill(sigmaPionP);
                                            hsigmaKaonEta[0][0][l]->Fill(sigmaKaonP);
                                            hsigmaProtonEta[0][0][l]->Fill(sigmaProtonP);
                                            hsigmaElectronEta[0][0][l]->Fill(sigmaElectronP);

                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2) || sigmaPionToFP == -999)
                                            {

                                                hmsqrEta[2][0][l]->Fill(msqrPlus);
                                                hPivsKeta[2][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPeta[2][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEeta[2][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKeta[2][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hsigmaPionEta[2][0][l]->Fill(sigmaPionP);
                                                hsigmaKaonEta[2][0][l]->Fill(sigmaKaonP);
                                                hsigmaProtonEta[2][0][l]->Fill(sigmaProtonP);
                                                hsigmaElectronEta[2][0][l]->Fill(sigmaElectronP);
                                            } // ToF or TPC
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFP != -999 && fabs(sigmaPionToFP) < 2)
                                            {
                                                hmsqrEta[1][0][l]->Fill(msqrPlus);
                                                hPivsKeta[1][0][l]->Fill(sigmaKaonP, sigmaPionP);
                                                hPivsPeta[1][0][l]->Fill(sigmaProtonP, sigmaPionP);
                                                hPivsEeta[1][0][l]->Fill(sigmaElectronP, sigmaPionP);
                                                hPvsKeta[1][0][l]->Fill(sigmaKaonP, sigmaProtonP);
                                                hsigmaPionEta[1][0][l]->Fill(sigmaPionP);
                                                hsigmaKaonEta[1][0][l]->Fill(sigmaKaonP);
                                                hsigmaProtonEta[1][0][l]->Fill(sigmaProtonP);
                                                hsigmaElectronEta[1][0][l]->Fill(sigmaElectronP);
                                            } // ToF only

                                        } // eta_pair window control
                                    }     // 9 loop for Eta Bin

                                } // both track +ve charge

                                if (fcharge[tr] < 0 || fcharge[tr2] < 0)
                                {

                                    if (fcharge[tr] < 0)
                                    {
                                        pN = fp[tr];
                                        ptN = fpT[tr];
                                        sigmaPionN = fnSigmaPion[tr];
                                        sigmaKaonN = fnSigmaKaon[tr];
                                        sigmaProtonN = fnSigmaProton[tr];
                                        sigmaElectronN = fnSigmaElectron[tr];

                                        sigmaPionToFN = fnSigmaPionToF[tr];
                                        sigmaKaonToFN = fnSigmaKaonToF[tr];
                                        sigmaProtonToFN = fnSigmaProtonToF[tr];

                                        msqrMinus = msqr1;
                                        log_dEdxN = log10(fdEdx[tr]);
                                        log_pN = log10(fp[tr]);
                                    }
                                    if (fcharge[tr2] < 0)
                                    {
                                        pN = fp[tr2];
                                        ptN = fpT[tr2];
                                        sigmaPionN = fnSigmaPion[tr2];
                                        sigmaKaonN = fnSigmaKaon[tr2];
                                        sigmaProtonN = fnSigmaProton[tr2];
                                        sigmaElectronN = fnSigmaElectron[tr2];

                                        sigmaPionToFN = fnSigmaPionToF[tr2];
                                        sigmaKaonToFN = fnSigmaKaonToF[tr2];
                                        sigmaProtonToFN = fnSigmaProtonToF[tr2];

                                        msqrMinus = msqr2;
                                        log_dEdxN = log10(fdEdx[tr2]);
                                        log_pN = log10(fp[tr2]);
                                    }

                                    for (int l = 0; l < 5; l++)
                                    {
                                        if (pT_pair >= pT[l] && pT_pair < pT[l + 1] && eta_pair > 0)
                                        {
                                            h_dEdxVsp_pt[1][0][l]->Fill(log_pN, log_dEdxN);
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            // if (sigmaPionToFN != -999)
                                            {
                                                hsigmaPion_TpcVsToF_pt[0][1][0][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_pt[0][1][0][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_pt[0][1][0][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // ToF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_pt[1][1][0][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_pt[1][1][0][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_pt[1][1][0][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_Cut

                                            // TPC
                                            hmsqrPt[0][1][0][l]->Fill(msqrMinus);
                                            hsigmaPionPt[0][1][0][l]->Fill(sigmaPionN);
                                            hsigmaProtonPt[0][1][0][l]->Fill(sigmaProtonN);
                                            hsigmaKaonPt[0][1][0][l]->Fill(sigmaKaonN);
                                            hsigmaElectronPt[0][1][0][l]->Fill(sigmaElectronN);
                                            hPivsKpt[0][1][0][l]->Fill(sigmaKaonN, sigmaPionN);
                                            hPivsPpt[0][1][0][l]->Fill(sigmaProtonN, sigmaPionN);
                                            hPivsEpt[0][1][0][l]->Fill(sigmaElectronN, sigmaPionN);
                                            hPvsKpt[0][1][0][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            hpP[0][1][0][l]->Fill(pN);
                                            hptP[0][1][0][l]->Fill(ptN);

                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            //    if ((sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2) || sigmaPionToFN == -999)
                                            {
                                                hmsqrPt[2][1][0][l]->Fill(msqrMinus);
                                                hsigmaPionPt[2][1][0][l]->Fill(sigmaPionN);
                                                hsigmaProtonPt[2][1][0][l]->Fill(sigmaProtonN);
                                                hsigmaKaonPt[2][1][0][l]->Fill(sigmaKaonN);
                                                hsigmaElectronPt[2][1][0][l]->Fill(sigmaElectronN);
                                                hPivsKpt[2][1][0][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPpt[2][1][0][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEpt[2][1][0][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKpt[2][1][0][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            } // TPC or TOF
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //  if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hmsqrPt[1][1][0][l]->Fill(msqrMinus);
                                                hsigmaPionPt[1][1][0][l]->Fill(sigmaPionN);
                                                hsigmaProtonPt[1][1][0][l]->Fill(sigmaProtonN);
                                                hsigmaKaonPt[1][1][0][l]->Fill(sigmaKaonN);
                                                hsigmaElectronPt[1][1][0][l]->Fill(sigmaElectronN);
                                                hPivsKpt[1][1][0][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPpt[1][1][0][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEpt[1][1][0][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKpt[1][1][0][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            } // ToF only

                                        } // control statment for pT binning and eta_pair>0

                                        if (Minv >= M[l] && Minv < M[l + 1] && eta_pair > 0)
                                        {
                                            h_dEdxVsp_M[1][0][l]->Fill(log_pN, log_dEdxN);

                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            // if (sigmaPionToFN != -999)
                                            {
                                                hsigmaPion_TpcVsToF_M[0][1][0][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_M[0][1][0][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_M[0][1][0][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //   if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_M[1][1][0][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_M[1][1][0][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_M[1][1][0][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_Cut
                                              // TPC
                                            hmsqrM[0][1][0][l]->Fill(msqrMinus);
                                            hsigmaPionM[0][1][0][l]->Fill(sigmaPionN);
                                            hsigmaProtonM[0][1][0][l]->Fill(sigmaProtonN);
                                            hsigmaKaonM[0][1][0][l]->Fill(sigmaKaonN);
                                            hsigmaElectronM[0][1][0][l]->Fill(sigmaElectronN);
                                            hPivsKm[0][1][0][l]->Fill(sigmaKaonN, sigmaPionN);
                                            hPivsPm[0][1][0][l]->Fill(sigmaProtonN, sigmaPionN);
                                            hPivsEm[0][1][0][l]->Fill(sigmaElectronN, sigmaPionN);
                                            hPvsKm[0][1][0][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            hpM[0][1][0][l]->Fill(pN);
                                            hptM[0][1][0][l]->Fill(ptN);
                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            //   if ((sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2) || sigmaPionToFN == -999)
                                            {
                                                hmsqrM[2][1][0][l]->Fill(msqrMinus);
                                                hsigmaPionM[2][1][0][l]->Fill(sigmaPionN);
                                                hsigmaProtonM[2][1][0][l]->Fill(sigmaProtonN);
                                                hsigmaKaonM[2][1][0][l]->Fill(sigmaKaonN);
                                                hsigmaElectronM[2][1][0][l]->Fill(sigmaElectronN);
                                                hPivsKm[2][1][0][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPm[2][1][0][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEm[2][1][0][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKm[2][1][0][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            } // TPC or TOF
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //   if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hmsqrM[1][1][0][l]->Fill(msqrMinus);
                                                hsigmaPionM[1][1][0][l]->Fill(sigmaPionN);
                                                hsigmaProtonM[1][1][0][l]->Fill(sigmaProtonN);
                                                hsigmaKaonM[1][1][0][l]->Fill(sigmaKaonN);
                                                hsigmaElectronM[1][1][0][l]->Fill(sigmaElectronN);
                                                hPivsKm[1][1][0][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPm[1][1][0][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEm[1][1][0][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKm[1][1][0][l]->Fill(sigmaKaonN, sigmaProtonN);

                                            } // TOF Only

                                        } // control statment for Minv bining and eta_pair>0

                                        if (pT_pair >= pT[l] && pT_pair < pT[l + 1] && eta_pair < 0)
                                        {
                                            h_dEdxVsp_pt[1][1][l]->Fill(log_pN, log_dEdxN);
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            // if (sigmaPionToFN != -999)
                                            {
                                                hsigmaPion_TpcVsToF_pt[0][1][1][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_pt[0][1][1][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_pt[0][1][1][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_No Cut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //    if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_pt[1][1][1][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_pt[1][1][1][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_pt[1][1][1][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_Cut

                                            hmsqrPt[0][1][1][l]->Fill(msqrMinus);
                                            hsigmaPionPt[0][1][1][l]->Fill(sigmaPionN);
                                            hsigmaProtonPt[0][1][1][l]->Fill(sigmaProtonN);
                                            hsigmaKaonPt[0][1][1][l]->Fill(sigmaKaonN);
                                            hsigmaElectronPt[0][1][1][l]->Fill(sigmaElectronN);
                                            hPivsKpt[0][1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                            hPivsPpt[0][1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                            hPivsEpt[0][1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                            hPvsKpt[0][1][1][l]->Fill(sigmaKaonN, sigmaProtonN);

                                            hpP[0][1][1][l]->Fill(pN);
                                            hptP[0][1][1][l]->Fill(ptN);
                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2) || sigmaPionToFN == -999)
                                            {

                                                hmsqrPt[2][1][1][l]->Fill(msqrMinus);
                                                hsigmaPionPt[2][1][1][l]->Fill(sigmaPionN);
                                                hsigmaProtonPt[2][1][1][l]->Fill(sigmaProtonN);
                                                hsigmaKaonPt[2][1][1][l]->Fill(sigmaKaonN);
                                                hsigmaElectronPt[2][1][1][l]->Fill(sigmaElectronN);
                                                hPivsKpt[2][1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPpt[2][1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEpt[2][1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKpt[2][1][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                                hpP[2][1][1][l]->Fill(pN);
                                                hptP[2][1][1][l]->Fill(ptN);
                                            } // TPC or TOF
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //  if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hmsqrPt[1][1][1][l]->Fill(msqrMinus);
                                                hsigmaPionPt[1][1][1][l]->Fill(sigmaPionN);
                                                hsigmaProtonPt[1][1][1][l]->Fill(sigmaProtonN);
                                                hsigmaKaonPt[1][1][1][l]->Fill(sigmaKaonN);
                                                hsigmaElectronPt[1][1][1][l]->Fill(sigmaElectronN);
                                                hPivsKpt[1][1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPpt[1][1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEpt[1][1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKpt[1][1][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                                hpP[1][1][1][l]->Fill(pN);
                                                hptP[1][1][1][l]->Fill(ptN);
                                            } // ToF only

                                        } // control statement for pT_bining and eta_pair<0

                                        if (Minv >= M[l] && Minv < M[l + 1] && eta_pair < 0)
                                        {
                                            h_dEdxVsp_M[1][1][l]->Fill(log_pN, log_dEdxN);

                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            //  if (sigmaPionToFN != -999)
                                            {
                                                hsigmaPion_TpcVsToF_M[0][1][1][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_M[0][1][1][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_M[0][1][1][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // ToF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            //     if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_M[1][1][1][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_M[1][1][1][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_M[1][1][1][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_Cut
                                            // TPC
                                            hmsqrM[0][1][1][l]->Fill(msqrMinus);
                                            hsigmaPionM[0][1][1][l]->Fill(sigmaPionN);
                                            hsigmaProtonM[0][1][1][l]->Fill(sigmaProtonN);
                                            hsigmaKaonM[0][1][1][l]->Fill(sigmaKaonN);
                                            hsigmaElectronM[0][1][1][l]->Fill(sigmaElectronN);
                                            hPivsKm[0][1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                            hPivsPm[0][1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                            hPivsEm[0][1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                            hPvsKm[0][1][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            hpM[0][1][1][l]->Fill(pN);
                                            hptM[0][1][1][l]->Fill(ptN);

                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2) || sigmaPionToFN == -999)
                                            {
                                                hmsqrM[2][1][1][l]->Fill(msqrMinus);
                                                hsigmaPionM[2][1][1][l]->Fill(sigmaPionN);
                                                hsigmaProtonM[2][1][1][l]->Fill(sigmaProtonN);
                                                hsigmaKaonM[2][1][1][l]->Fill(sigmaKaonN);
                                                hsigmaElectronM[2][1][1][l]->Fill(sigmaElectronN);
                                                hPivsKm[2][1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPm[2][1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEm[2][1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKm[2][1][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                                hpM[2][1][1][l]->Fill(pN);
                                                hptM[2][1][1][l]->Fill(ptN);
                                            } // ToF or TPC

                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hmsqrM[1][1][1][l]->Fill(msqrMinus);
                                                hsigmaPionM[1][1][1][l]->Fill(sigmaPionN);
                                                hsigmaProtonM[1][1][1][l]->Fill(sigmaProtonN);
                                                hsigmaKaonM[1][1][1][l]->Fill(sigmaKaonN);
                                                hsigmaElectronM[1][1][1][l]->Fill(sigmaElectronN);
                                                hPivsKm[1][1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPm[1][1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEm[1][1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKm[1][1][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                                hpM[1][1][1][l]->Fill(pN);
                                                hptM[1][1][1][l]->Fill(ptN);
                                            } // ToF only

                                        } // control statement for Minv bining with eta_pair<0
                                    }     // 5 pT or Minv bin loop

                                    for (int l = 0; l < 9; l++)
                                    {
                                        if (eta_pair >= eta_range[l] && eta_pair < eta_range[l + 1])
                                        {
                                            if (fnSigmaPionToF[tr] != -999 && fnSigmaPionToF[tr2] != -999)
                                            //    if (sigmaPionToFN != -999)
                                            {
                                                hsigmaPion_TpcVsToF_eta[0][1][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_eta[0][1][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_eta[0][1][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // ToF_NoCut
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hsigmaPion_TpcVsToF_eta[1][1][l]->Fill(sigmaPionToFN, sigmaPionN);
                                                hsigmaKaon_TpcVsToF_eta[1][1][l]->Fill(sigmaKaonToFN, sigmaKaonN);
                                                hsigmaProton_TpcVsToF_eta[1][1][l]->Fill(sigmaProtonToFN, sigmaProtonN);
                                            } // TOF_NoCut
                                            // TPC
                                            hmsqrEta[0][1][l]->Fill(msqrMinus);
                                            hsigmaPionEta[0][1][l]->Fill(sigmaPionN);
                                            hsigmaKaonEta[0][1][l]->Fill(sigmaKaonN);
                                            hsigmaProtonEta[0][1][l]->Fill(sigmaProtonN);
                                            hsigmaElectronEta[0][1][l]->Fill(sigmaElectronN);
                                            hPivsKeta[0][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                            hPivsPeta[0][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                            hPivsEeta[0][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                            hPvsKeta[0][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            // TPC or TOF
                                            if ((((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) || fnSigmaPionToF[tr] == -999)) && (((fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2) || fnSigmaPionToF[tr2] == -999)))
                                            // if ((sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2) || sigmaPionToFN == -999)
                                            {

                                                hmsqrEta[2][1][l]->Fill(msqrMinus);
                                                hsigmaPionEta[2][1][l]->Fill(sigmaPionN);
                                                hsigmaKaonEta[2][1][l]->Fill(sigmaKaonN);
                                                hsigmaProtonEta[2][1][l]->Fill(sigmaProtonN);
                                                hsigmaElectronEta[2][1][l]->Fill(sigmaElectronN);
                                                hPivsKeta[2][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPeta[2][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEeta[2][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKeta[2][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            } // TOF or TPC
                                            if ((fnSigmaPionToF[tr] != -999 && fabs(fnSigmaPionToF[tr]) < 2) && (fnSigmaPionToF[tr2] != -999 && fabs(fnSigmaPionToF[tr2]) < 2))
                                            // if (sigmaPionToFN != -999 && fabs(sigmaPionToFN) < 2)
                                            {
                                                hmsqrEta[1][1][l]->Fill(msqrMinus);
                                                hsigmaPionEta[1][1][l]->Fill(sigmaPionN);
                                                hsigmaKaonEta[1][1][l]->Fill(sigmaKaonN);
                                                hsigmaProtonEta[1][1][l]->Fill(sigmaProtonN);
                                                hsigmaElectronEta[1][1][l]->Fill(sigmaElectronN);
                                                hPivsKeta[1][1][l]->Fill(sigmaKaonN, sigmaPionN);
                                                hPivsPeta[1][1][l]->Fill(sigmaProtonN, sigmaPionN);
                                                hPivsEeta[1][1][l]->Fill(sigmaElectronN, sigmaPionN);
                                                hPvsKeta[1][1][l]->Fill(sigmaKaonN, sigmaProtonN);
                                            } // ToF only

                                        } // control statement for eta_bining
                                    }     // 9 loop for nine Eta bin
                                }         // both track -ve charge

                            } // cone cut
                        }     // opposite charge control statement
                    }         // Track 2 quality cut
                }             // Track 2 Loop

            } // Track 1 qualtiy cut
            // cout << "Hello" << endl;
        } // Track 1 Loop

        // cout << jentry << "jentry" << endl;
    } // Event Loop
    cout << "Pion paring Loop" << endl;

} // Iff2017:Loop()

Iff2017::Iff2017(char *ifile)
{
    fChain = new TChain("ftree");
    fChain->Add(ifile);
    Init();
}

Iff2017::~Iff2017()
{
    if (!fChain)
        return;
    delete fChain->GetCurrentFile();
}

Int_t Iff2017::GetEntry(Long64_t entry)
{
    if (!fChain)
        return 0;
    return fChain->GetEntry(entry);
}
Long64_t Iff2017::LoadTree(Long64_t entry)
{
    if (!fChain)
        return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0)
        return centry;
    if (!fChain->InheritsFrom(TChain::Class()))
        return centry;
    TChain *chain = (TChain *)fChain;
    if (chain->GetTreeNumber() != fCurrent)
    {
        fCurrent = chain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void Iff2017::Init()
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers

    fCurrent = -1;
    fChain->SetMakeClass(1);
    // fChain->SetBranchAddress("frefmult", &frefmult, &b_frefmult);
    fChain->SetBranchAddress("fmaxpar", &fmaxpar, &b_fmaxpar);
    fChain->SetBranchAddress("ffillNum", &ffillNum, &b_ffillNum);
    fChain->SetBranchAddress("frunNum", &frunNum, &b_frunNum);
    fChain->SetBranchAddress("fspinconfig", &fspinconfig, &b_fspinconfig);
    fChain->SetBranchAddress("ftrigger", &ftrigger, &b_ftrigger);
    fChain->SetBranchAddress("fVZ", &fVZ, &b_fVZ);
    fChain->SetBranchAddress("fevTime", &fevTime, &b_fevTime);
    // fChain->SetBranchAddress("fverRank", &fverRank, &b_fverRank);
    fChain->SetBranchAddress("fpT", fpT, &b_fpT);
    fChain->SetBranchAddress("fp", fp, &b_fp);
    fChain->SetBranchAddress("feta", feta, &b_feta);
    fChain->SetBranchAddress("fphi", fphi, &b_fphi);
    fChain->SetBranchAddress("fcharge", fcharge, &b_fcharge);
    fChain->SetBranchAddress("fnSigmaPion", fnSigmaPion, &b_fnSigmaPion);
    fChain->SetBranchAddress("fnSigmaKaon", fnSigmaKaon, &b_fnSigmaKaon);
    fChain->SetBranchAddress("fnSigmaProton", fnSigmaProton, &b_fnSigmaProton);
    fChain->SetBranchAddress("fnSigmaElectron", fnSigmaElectron, &b_fnSigmaElectron);

    fChain->SetBranchAddress("fnSigmaPionToF", fnSigmaPionToF, &b_fnSigmaPionToF);
    fChain->SetBranchAddress("fnSigmaKaonToF", fnSigmaKaonToF, &b_fnSigmaKaonToF);
    fChain->SetBranchAddress("fnSigmaProtonToF", fnSigmaProtonToF, &b_fnSigmaProtonToF);
    fChain->SetBranchAddress("fnSigmaElectronToF", fnSigmaElectronToF, &b_fnSigmaElectronToF);

    fChain->SetBranchAddress("fdEdx", fdEdx, &b_fdEdx);
    fChain->SetBranchAddress("fdca", fdca, &b_fdca);
    fChain->SetBranchAddress("ffitPts", ffitPts, &b_ffitPts);
    fChain->SetBranchAddress("ffitPtsPoss", ffitPtsPoss, &b_ffitPtsPoss);
    fChain->SetBranchAddress("fhitsdedx", fhitsdedx, &b_fhitsdedx);
    fChain->SetBranchAddress("fBetaToF", fBetaToF, &b_fBetaToF);
    fChain->SetBranchAddress("fvzVpd", &fvzVpd);
    Notify();
} // Init() ended

Bool_t Iff2017::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.
    return kTRUE;
}

void Iff2017::Finish(char *ofile)
{
    TFile *fout = new TFile(ofile, "recreate");
    for (int ndet = 0; ndet < 3; ndet++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                for (int i = 0; i < 5; i++)
                {

                    hmsqrPt[ndet][j][k][i]->SetDirectory(fout);
                    hmsqrM[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaPionM[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaKaonM[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaProtonM[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaElectronM[ndet][j][k][i]->SetDirectory(fout);

                    hsigmaPionPt[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaKaonPt[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaProtonPt[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaElectronPt[ndet][j][k][i]->SetDirectory(fout);

                    hsigmaPion_TpcVsToF_pt[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaKaon_TpcVsToF_pt[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaProton_TpcVsToF_pt[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaPion_TpcVsToF_M[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaKaon_TpcVsToF_M[ndet][j][k][i]->SetDirectory(fout);
                    hsigmaProton_TpcVsToF_M[ndet][j][k][i]->SetDirectory(fout);

                    hPivsKpt[ndet][j][k][i]->SetDirectory(fout);
                    hPvsKpt[ndet][j][k][i]->SetDirectory(fout);
                    hPivsPpt[ndet][j][k][i]->SetDirectory(fout);
                    hPivsEpt[ndet][j][k][i]->SetDirectory(fout);
                    hPivsKm[ndet][j][k][i]->SetDirectory(fout);
                    hPvsKm[ndet][j][k][i]->SetDirectory(fout);
                    hPivsPm[ndet][j][k][i]->SetDirectory(fout);
                    hPivsEm[ndet][j][k][i]->SetDirectory(fout);

                    hpM[ndet][j][k][i]->SetDirectory(fout);
                    hptM[ndet][j][k][i]->SetDirectory(fout);
                    hpP[ndet][j][k][i]->SetDirectory(fout);
                    hptP[ndet][j][k][i]->SetDirectory(fout);
                    h_dEdxVsp_pt[j][k][i]->SetDirectory(fout);
                    h_dEdxVsp_M[j][k][i]->SetDirectory(fout);

                    hmsqrPt[ndet][j][k][i]->Write();
                    hmsqrM[ndet][j][k][i]->Write();

                    hsigmaPionM[ndet][j][k][i]->Write();
                    hsigmaKaonM[ndet][j][k][i]->Write();
                    hsigmaProtonM[ndet][j][k][i]->Write();
                    hsigmaElectronM[ndet][j][k][i]->Write();

                    hsigmaPionPt[ndet][j][k][i]->Write();
                    hsigmaKaonPt[ndet][j][k][i]->Write();
                    hsigmaProtonPt[ndet][j][k][i]->Write();
                    hsigmaElectronPt[ndet][j][k][i]->Write();

                    hsigmaPion_TpcVsToF_pt[ndet][j][k][i]->Write();
                    hsigmaKaon_TpcVsToF_pt[ndet][j][k][i]->Write();
                    hsigmaProton_TpcVsToF_pt[ndet][j][k][i]->Write();
                    hsigmaPion_TpcVsToF_M[ndet][j][k][i]->Write();
                    hsigmaKaon_TpcVsToF_M[ndet][j][k][i]->Write();
                    hsigmaProton_TpcVsToF_M[ndet][j][k][i]->Write();

                    hPvsKpt[ndet][j][k][i]->Write();
                    hPivsKpt[ndet][j][k][i]->Write();
                    hPivsPpt[ndet][j][k][i]->Write();
                    hPivsEpt[ndet][j][k][i]->Write();
                    hPivsKm[ndet][j][k][i]->Write();
                    hPvsKm[ndet][j][k][i]->Write();
                    hPivsPm[ndet][j][k][i]->Write();
                    hPivsEm[ndet][j][k][i]->Write();

                    hpP[ndet][j][k][i]->Write();
                    hptP[ndet][j][k][i]->Write();
                    hpM[ndet][j][k][i]->Write();
                    hptM[ndet][j][k][i]->Write();
                    h_dEdxVsp_pt[j][k][i]->Write();
                    h_dEdxVsp_M[j][k][i]->Write();
                }
            }
        }
    }
    for (int ndet = 0; ndet < 3; ndet++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int i = 0; i < 9; i++)
            {
                hmsqrEta[ndet][j][i]->SetDirectory(fout);
                hsigmaPionEta[ndet][j][i]->SetDirectory(fout);
                hsigmaKaonEta[ndet][j][i]->SetDirectory(fout);
                hsigmaProtonEta[ndet][j][i]->SetDirectory(fout);
                hsigmaElectronEta[ndet][j][i]->SetDirectory(fout);
                hPvsKeta[ndet][j][i]->SetDirectory(fout);
                hPivsKeta[ndet][j][i]->SetDirectory(fout);
                hPivsPeta[ndet][j][i]->SetDirectory(fout);
                hPivsEeta[ndet][j][i]->SetDirectory(fout);

                hsigmaPion_TpcVsToF_eta[ndet][j][i]->SetDirectory(fout);
                hsigmaKaon_TpcVsToF_eta[ndet][j][i]->SetDirectory(fout);
                hsigmaProton_TpcVsToF_eta[ndet][j][i]->SetDirectory(fout);

                hmsqrEta[ndet][j][i]->Write();
                hsigmaPionEta[ndet][j][i]->Write();
                hsigmaKaonEta[ndet][j][i]->Write();
                hsigmaProtonEta[ndet][j][i]->Write();
                hsigmaElectronEta[ndet][j][i]->Write();
                hPvsKeta[ndet][j][i]->Write();
                hPivsKeta[ndet][j][i]->Write();
                hPivsPeta[ndet][j][i]->Write();
                hPivsEeta[ndet][j][i]->Write();

                hsigmaPion_TpcVsToF_eta[ndet][j][i]->Write();
                hsigmaKaon_TpcVsToF_eta[ndet][j][i]->Write();
                hsigmaProton_TpcVsToF_eta[ndet][j][i]->Write();
            }
        }
    } // TPC, ToF , TPC or ToF loop
    hnSigmaPionTPC->SetDirectory(fout);
    hnSigmaPionTPC->Write();
    hnSigmaPionToF->SetDirectory(fout);
    hnSigmaPionToF->Write();
    h_dEdx_p->SetDirectory(fout);
    h_dEdx_p->Write();
    fout->Close();
} // Finish ended;

void Iff2017::Show(Long64_t entry)
{
    // eventTime-> Write();
    // pions->Write();
    //  Print contents of entry.
    //  If entry is not specified, print current entry
    if (!fChain)
        return;
    fChain->Show(entry);
}

Int_t Iff2017::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
ClassImp(Iff2017);
