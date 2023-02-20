

/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 22 17:04:53 2010 by ROOT version 5.22/00
// from TTree ftree/LRC
// found on file: 00E5DE07942C320D8D35F5F8AD40C261_7total.root
//////////////////////////////////////////////////////////

#ifndef Iff2017_h
#define Iff2017_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include <iostream>
#include <fstream>

using namespace std;

class Iff2017
{
public:
    TChain *fChain; //! pointer to the analyzed TTree or TChain
    Int_t fCurrent; //! current Tree number in a TChain

    // Declaration of leaf types
    Int_t frefmult;
    Int_t fmaxpar;
    Int_t ffillNum;
    Int_t frunNum;
    Int_t fspinconfig;
    // Int_t           ftrigger;
    Double_t fVZ;
    Double_t fvzVpd; // Included by Babu 11/22/2019
    unsigned int fevTime;
    Double_t fverRank;                 // included by Babu  11/22/2019
    Double_t fpT[2422];                //[fmaxpar]
    Double_t fp[2422];                 //[fmaxpar]
    Double_t feta[2422];               //[fmaxpar]
    Double_t fphi[2422];               //[fmaxpar]
    Short_t fcharge[2422];             //[fmaxpar]
    Double_t fnSigmaPion[2422];        //[fmaxpar]
    Double_t fnSigmaKaon[2422];        //[fmaxpar]
    Double_t fnSigmaProton[2422];      //[fmaxpar]
    Double_t fnSigmaElectron[2422];    //[fmaxpar]
    Double_t fnSigmaPionToF[2422];     //[fmaxpar]
    Double_t fnSigmaKaonToF[2422];     //[fmaxpar]
    Double_t fnSigmaProtonToF[2422];   //[fmaxpar]
    Double_t fnSigmaElectronToF[2422]; //[fmaxpar]
    Double_t fdEdx[2422];              //[fmaxpar]
    Double_t fdca[2422];               //[fmaxpar]
    UShort_t ffitPts[2422];            //[fmaxpar]
    UShort_t ffitPtsPoss[2422];        //[fmaxpar]
    UShort_t fhitsdedx[2422];          //[fmaxpar]
    Double_t fBetaToF[2422];           //[fmaxpar]
    vector<int> *ftrigger;
    int trigId;

    // List of branches
    TBranch *b_frefmult;    //!
    TBranch *b_fmaxpar;     //!
    TBranch *b_ffillNum;    //!
    TBranch *b_frunNum;     //!
    TBranch *b_fspinconfig; //!
    TBranch *b_ftrigger;    //!
    TBranch *b_fVZ;
    TBranch *b_fvzVpd;
    TBranch *b_fevTime;
    TBranch *b_fverRank;           //!
    TBranch *b_fpT;                //!
    TBranch *b_fp;                 //!
    TBranch *b_feta;               //!
    TBranch *b_fphi;               //!
    TBranch *b_fcharge;            //!
    TBranch *b_fnSigmaPion;        //!
    TBranch *b_fnSigmaKaon;        //!
    TBranch *b_fnSigmaProton;      //!
    TBranch *b_fnSigmaElectron;    //!
    TBranch *b_fnSigmaPionToF;     //!
    TBranch *b_fnSigmaKaonToF;     //!
    TBranch *b_fnSigmaProtonToF;   //!
    TBranch *b_fnSigmaElectronToF; //!
    TBranch *b_fdEdx;              //!
    TBranch *b_fdca;               //!
    TBranch *b_ffitPts;            //!
    TBranch *b_ffitPtsPoss;        //!
    TBranch *b_fhitsdedx;          //!
    TBranch *b_fBetaToF;           //!

    Iff2017(char *ifile);
    virtual ~Iff2017();
    virtual Int_t Cut(Long64_t entry);
    virtual Int_t GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void Init();
    virtual void Loop();
    virtual Bool_t Notify();
    virtual void Finish(char *ofile);
    virtual void Show(Long64_t entry = -1);
    // virtual void     loadPol();
    inline void setHome(TDirectory *pHome)
    {
        cout << "setting home.." << endl;
        home = pHome;
    }

    TH1D *hsigmaPionM[3][2][2][5];
    TH1D *hsigmaKaonM[3][2][2][5];
    TH1D *hsigmaProtonM[3][2][2][5];
    TH1D *hsigmaElectronM[3][2][2][5];
    TH1D *hmsqrM[3][2][2][5];

    TH1D *hsigmaPionPt[3][2][2][5];
    TH1D *hsigmaKaonPt[3][2][2][5];
    TH1D *hsigmaProtonPt[3][2][2][5];
    TH1D *hsigmaElectronPt[3][2][2][5];
    TH1D *hmsqrPt[3][2][2][5];

    TH2D *hPivsKpt[3][2][2][5];
    TH2D *hPivsPpt[3][2][2][5];
    TH2D *hPivsEpt[3][2][2][5];

    TH2D *hPivsKm[3][2][2][5];
    TH2D *hPivsPm[3][2][2][5];
    TH2D *hPivsEm[3][2][2][5];

    TH2D *hPvsKpt[3][2][2][5];
    TH2D *hPvsKm[3][2][2][5];

    TH1D *hsigmaPionEta[3][2][9];
    TH1D *hsigmaKaonEta[3][2][9];
    TH1D *hsigmaElectronEta[3][2][9];
    TH1D *hsigmaProtonEta[3][2][9];
    TH1D *hmsqrEta[3][2][9];

    TH2D *hPivsKeta[3][2][9];
    TH2D *hPivsPeta[3][2][9];
    TH2D *hPivsEeta[3][2][9];
    TH2D *hPvsKeta[3][2][9];

    TH1D *hpM[3][2][2][5];
    TH1D *hptM[3][2][2][5];

    TH1D *hpP[3][2][2][5];
    TH1D *hptP[3][2][2][5];
    TH2D *h_dEdxVsp_pt[2][2][5];
    TH2D *h_dEdxVsp_M[2][2][5];
    TH2D *h_dEdx_p;

    TH2D *hsigmaPion_TpcVsToF_pt[3][2][2][5];
    TH2D *hsigmaKaon_TpcVsToF_pt[3][2][2][5];
    TH2D *hsigmaProton_TpcVsToF_pt[3][2][2][5];
    TH2D *hsigmaPion_TpcVsToF_M[3][2][2][5];
    TH2D *hsigmaKaon_TpcVsToF_M[3][2][2][5];
    TH2D *hsigmaProton_TpcVsToF_M[3][2][2][5];
   
    TH2D *hsigmaPion_TpcVsToF_eta[3][2][9];
    TH2D *hsigmaKaon_TpcVsToF_eta[3][2][9];
    TH2D *hsigmaProton_TpcVsToF_eta[3][2][9];

    TH1D *hnSigmaPionTPC;
    TH1D *hnSigmaPionToF;


protected:
    TNtuple *
        ntuple1;
    TNtuple *ntuple2;
    TNtuple *ntuple3;
    TNtuple *ntuple4;
    TNtuple *ntuple5;
    TNtuple *ntuple6;
    TNtuple *ntuple7;
    TNtuple *ntuple8;
    TNtuple *piplus;
    TNtuple *piminus;
    TNtuple *pions;
    TNtuple *allTrig;
    TDirectory *home;

private:
    ClassDef(Iff2017, 1);
};

#endif
