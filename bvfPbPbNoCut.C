#define bvfPbPbNoCut_cxx
#include "bvfPbPbNoCut.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasUtils.h"
#include "/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasStyle.h"
#include "/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasLabels.h"
#include "/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasStyle.C"

#ifdef __CLING__
// these are not headers - do not treat them as such - needed for ROOT6
#include "/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasLabels.C"
#include "/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasUtils.C"
#endif

#ifdef __CINT__
    gROOT->LoadMacro("/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasLabels.C");
    gROOT->LoadMacro("/Volumes/GoogleDrive/My Drive/ATLAS/atlasstyle-00-04-02/AtlasUtils.C");
#endif

const Float_t FCal_range[] = {0,0.063719,0.14414,0.289595,0.525092,0.87541,1.36875,2.04651,2.98931,5}; // fcal_cuts options
const Float_t Weight[] = {6.7890E+07*6.1692E-06, 6.3996E+05*5.8420E-05, 4.7192E+03*1.1270E-04, 9.2038E-05*2.6602E+01};
const int myColor[] = {kBlue, kViolet, kMagenta, kPink, kOrange, kYellow, kSpring, kTeal, kCyan, kAzure, kGray, kGray+1, kGray+3 };

Float_t Eta_range[] = {-2.1,-1.5,-0.9,-0.3,0.3,0.9,1.5,2.1};
const int Eta_N = sizeof(Eta_range)/sizeof(float)-1;

const int count_max = 1500;
const int pt_min = 50;
const int pt_max = 600;
const int pt_bin = 8;
const int min_dist = 0;
const int max_dist = 20;
const int dist_bin = 50;
const int bHadCut = 5;
const int cHadCut = 5;
const Float_t eta_selection = 2.1;

const int kTruth = kBlue;
const int kReco = kRed;

const bool unique_B = false;
const int which = 1;

const int dist_bin2 = 20;
const float low_log = 1e-4;
const float high_log = 0.6;
const float min_dist2 = 1e-5;
const float max_dist2 = 0.3;
const float mcprob_cut = 0.75;

const int dist_bin3 = 100;

const int highd0 = 6;
const int highz0 = 10;
const int highd0sig = 200;
const int highz0sig = 100;
const int perigeebins = 50;
const int kB = kRed;
const int kF = kBlue;

const int cet[] = {0,2,2,5,5,8};//selected centrality sections
const int cet_N = (sizeof(cet)/sizeof(int))/2;

const char track_selections[3][50] = {"No Selection","IP3D Selection","SV1 Selection"};
const int conesize = 4;

const char* dataType = "NoCut";

//enum TRKORIGIN{ PUFAKE=-1,
//    FROMB,
//    FROMC,
//    FRAG,
//    GEANT };

//truth track from B not from C decayed from B: FROMB (0)
//truth track from C from C decayed from B: FROMC (1)
//truth track from C from C: FROMC(1)
//truth track from fragmentation: FRAG(2)


//Maximum dR between B tracks and jet axis:
//for all tracks within dR < 1.0 of jet axis, those whose truth match originates from a B hadron are measured dR, and the furthest dR is recorded.

//All the b jet classified this way should have a B hadron and thus all FROMC should be a C from B. 
const bool jetTruth = true;
void bvfPbPbNoCut::Loop()
{
//   In a ROOT session, you can do:
//      root> .L bvfPbPbNoCut.C
//      root> bvfPbPbNoCut t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  Float_t ptbins[pt_bin+1];   
  float initial = log(pt_min);
  float incre = log(pt_max/pt_min)/pt_bin;
  for (int i = 0; i < (pt_bin+1);i++)  {
      ptbins[i] = TMath::Power(TMath::E(),initial); 
      initial = initial + incre;
      //cout << ptbins[i] << endl;
  }

if (fChain == 0) return;
   const char* uniqueness = unique_B?"_Unique_B":"";
   const char* TorR = jetTruth?"truth":"reco";
   const char* jetTruthness = jetTruth?"":"_allJet";
   const char* track_selection = track_selections[which];

   Long64_t nentries = fChain->GetEntries();
   std::string fname;
   float weight;

   TFile* out = TFile::Open(Form("PbPb%sNoCalrapidity%.1f%s%s_morehist.root", dataType, eta_selection, uniqueness, track_selection),"RECREATE");

   Long64_t nbytes = 0, nb = 0;

   //TH1F* truth_d0 = hotTH1F("truth_d0","Distribution of d0 for all truth tracks in b jets", 100, -10, 10,"","",kTruth,0.3,21,1,true);
   /*TH1F* reco_d0 = hotTH1F("reco_d0",Form("Distribution of d0 for all reco tracks in b jets %s", dataType), perigeebins, -highd0, highd0,"","",kReco,0.3,21,1,true);
  
   TH1F* reco_d0F = hotTH1F("reco_d0F", Form("Distribution of d0 for Fragmentation reco tracks in b jets track cone %d",conesize), perigeebins, -highd0, highd0,"","",kF,0.3,20,1,true);
   TH1F* reco_d0B = hotTH1F("reco_d0B", Form("Distribution of d0 for B reco tracks in b jets track cone %d",conesize), perigeebins, -highd0, highd0,"","",kB,0.3,20,1,true);
   //TH1F* truth_Z0 = hotTH1F("truth_d0","Distribution of d0 for all truth tracks in b jets", 100, -10, 10,"","",kTruth,0.3,21,1,true);
   TH1F* reco_z0 = hotTH1F("reco_z0",Form("Distribution of z0 for all reco tracks in b jets %s", dataType), perigeebins, -highz0, highz0,"","",kReco,0.3,21,1,true);

   TH1F* reco_z0F = hotTH1F("reco_z0F", Form("Distribution of z0 for Fragmentation reco tracks in b jets track cone %d",conesize), perigeebins, -highz0, highz0,"","",kF,0.3,20,1,true);
   TH1F* reco_z0B = hotTH1F("reco_z0B", Form("Distribution of z0 for B reco tracks in b jets track cone %d",conesize), perigeebins, -highz0, highz0,"","",kB,0.3,20,1,true);*/

   //TH1F* truth_sigd0 = hotTH1F("truth_d0","Distribution of d0 for all truth tracks in b jets", 100, -10, 10,"","",kTruth,0.3,21,1,true);
   /*TH1F* reco_d0sig = hotTH1F("reco_d0sig",Form("Distribution of d0 sigfor all reco tracks in b jets %s", dataType), perigeebins, -highd0sig, highd0sig,"","",kReco,0.3,21,1,true);

   //TH1F* truth_sigZ0 = hotTH1F("truth_d0","Distribution of d0 for all truth tracks in b jets", 100, -10, 10,"","",kTruth,0.3,21,1,true);
   TH1F* reco_z0sig = hotTH1F("reco_z0sig",Form("Distribution of z0 sig for all reco tracks in b jets %s", dataType), perigeebins, -highz0sig, highz0sig,"","",kReco,0.3,21,1,true);

   TH1F* reco_d0sig_f = hotTH1F("reco_d0sig_f",Form("Distribution of d0 sig for reco fragmentation tracks in b jets %s", dataType), perigeebins, -highd0sig, highd0sig,"","",kF,0.3,21,1,true);
   TH1F* reco_d0_f = hotTH1F("reco_d0_f",Form("Distribution of d0 for reco fragmentation tracks in b jets %s", dataType), perigeebins, -highd0, highd0,"","",kF,0.3,21,1,true);

   TH1F* reco_d0sig_b = hotTH1F("reco_d0sig_b",Form("Distribution of d0 sig for reco b tracks in b jets %s", dataType), perigeebins, -highd0sig, highd0sig,"","",kB,0.3,21,1,true);
   TH1F* reco_d0_b = hotTH1F("reco_d0_b",Form("Distribution of d0 for reco b tracks in b jets %s", dataType), perigeebins, -highd0, highd0,"","",kB,0.3,21,1,true);*/

   //TH1F* h;

   //TList* numTot = new TList();
   //TList* numB = new TList();
   //TList* numF = new TList();
   //TList* numBnoC = new TList();
   //TList* farDRB = new TList();
   TH1F* numTot[cet_N][pt_bin];
   TH1F* numB[cet_N][pt_bin];
   TH1F* numF[cet_N][pt_bin];
   TH1F* numBnoC[cet_N][pt_bin];
   TH1F* farDRB[cet_N][pt_bin];

   for (int i = 0; i < cet_N; i++){//loop over centrality
   	for (int j = 0; j < pt_bin; j++){//loop over pt
    numTot[i][j] = hotTH1F(Form("num_tot_cet_%d_pt_%d",i,j),Form("Distribution of total number of tracks in b jets %s", dataType), 11, -0.5, 10.5, "", "", myColor[j],0.3,21,1,true);
    numB[i][j] = hotTH1F(Form("num_b_cet_%d_pt_%d",i,j),Form("Distribution of number of b tracks in b jets %s", dataType), 11, -0.5, 10.5, "", "", myColor[j],0.3,21,1,true);
    numF[i][j] = hotTH1F(Form("num_f_cet_%d_pt_%d",i,j),Form("Distribution of number of fragmentation tracks in b jets %s", dataType), 11, -0.5, 10.5, "", "", myColor[j],0.3,21,1,true);
    numBnoC[i][j] = hotTH1F(Form("num_b_no_c_cet_%d_pt_%d",i,j),Form("Distribution of number of b but not c tracks in b jets %s", dataType), 11, -0.5, 10.5, "", "", myColor[j],0.3,21,1,true);
    farDRB[i][j] = hotTH1F(Form("fardRB_cet_%d_pt_%d",i,j),Form("Distribution of Furthest B track from Jet in b jets track %s", dataType), dist_bin3, 0, 1.0, "", "", myColor[j],0.3,21,1,false);
    //numTot->Add(num_tot);
    //numB->Add(num_b);
    //numF->Add(num_f);
    //numBnoC->Add(num_b_no_c);
    //farDRB->Add(fardRB);
	}
   }

//loop over for weight
   std::vector<float> wgsum;
   for (Long64_t jentry0 = 0; jentry0 < nentries; jentry0++){
    Long64_t ientry0 = LoadTree(jentry0);
    fChain->GetEntry(jentry0);
    if (ientry0 < 0) {
      cout << "invalid index" << endl;
      break;
    }
    //cout << "jentry" << jentry0 << endl;
    //cout << "ientry" << ientry0 << endl;
    if (ientry0 == 0) {
      wgsum.push_back(mcwg); 
      //cout << wgsum.size() << "; " << mcwg << endl;
    } else {
      //cout << wgsum[wgsum.size()-1] << endl;
      wgsum[wgsum.size()-1] = wgsum[wgsum.size()-1]+mcwg;
      //cout << wgsum.size() << ";" << mcwg << endl;
    }
   }

   cout << "number of slices:" << wgsum.size() << endl;
   for (int i = 0; i < wgsum.size(); i++){
    cout << "JZ Slice " << i+1 << ": " << wgsum[i] << endl;
   }

   int multiB = 0;
   int nJets = 0;
   int NJets = 0;

   //TH1F* l3d_truth = new TH1F("l3d_truth","l3d_truth");
   for (Long64_t jentry = 0; jentry < nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //float FCal_et = Fcal;
//loop over each jet to categorize them into different jet type
//rule, if there's no truth jet matched, pass the jet

      float FCal_et = Fcal;
      //cout << "FCal_et: " << FCal_et << endl;
      int central;
      for (int f = 0; f < cet_N; f++){
        if (FCal_et > FCal_range[9-cet[f*2+1]] && FCal_et < FCal_range[9-cet[f*2]]) central = f;
      }
        fname = fChain -> GetCurrentFile() -> GetName();
        for (int i = 0; i < wgsum.size(); i++){
        if (fname.find(Form("JZ%d",i+1))<fname.size()){
          weight = Weight[i]*mcwg/wgsum[i];
       }
      }

      for (int j = 0; j < jet_truthMatch -> size(); j++){

         int jetType = 0;//1 for B jet, 2 for C jet, 3 for T jet, 0 for light jet
         if (jet_truthMatch->at(j)!=1 && jetTruth) continue; 
         float jetPt = jetTruth?(jet_truthPt->at(j)*1e-3):(jet_pt->at(j)*1e-3);
         if (fabs(jet_truthEta->at(j))>2.1) continue;
         if (jetPt<pt_min) continue;

         //if (jet_cH_pdgId->at(j)[0]!=-99) jetType = 2;
         //if (jet_bH_pdgId->at(j)[0]!=-99) jetType = 1;

         if (jet_dRminToT->at(j) < 0.3) jetType = 3;
         //if (jet_nGhostTau) jetType = 3;
         if (jet_cH_pdgId->at(j)[0]!=-99) jetType = 2;
         if (jet_bH_pdgId->at(j)[0]!=-99) jetType = 1;

         //int nbhad_t = 0;//number of truth hadrons matched to each jet
         //int nchad_t = 0;
         //use the closest truth hadron matched to the jet

         if (jetType == 1){
            if (jet_bH_Lxy->at(j)[0] <= -99) {
               cout << "出大问题 (big problem) " << endl;
               continue;
            }
            if (!(jet_bH_Lxy->at(j)[0] >=0)) {
               cout << "Has no decay vertex" << endl;
               continue;
            }
            // 
            int tot_num = 0;
            int b_num = 0;
            int f_num = 0;
            int b_num_no_c = 0;
            //cout << "index: " << (int)(TMath::Log(jet_truthPt->at(j)*1e-3/pt_min)/incre) << endl;
            //h = (TH1F*)hsv_truth->At((int)(TMath::Log(jet_truthPt->at(j)*1e-3/pt_min)/incre));
            //h -> Fill(jet_bH_dRjet->at(j)[0],weight);
            if (jet_trk_orig->at(j).size()!=jet_trk_pt->at(j).size()) {
              cout << "different track selection!!!" << endl;
              break;
            }
            float dRB = -10;
            for (int k = 0; k < jet_trk_pt->at(j).size(); k++){
              //reco_d0->Fill(jet_trk_ip3d_d0->at(j)[k],weight);
              //reco_d0sig ->Fill(jet_trk_ip3d_d0sig->at(j)[k],weight);
              //reco_z0->Fill(jet_trk_ip3d_d0->at(j)[k],weight);
              //reco_z0sig ->Fill(jet_trk_ip3d_d0sig->at(j)[k],weight); 

              //if there's no truth match
              if (jet_trk_pdg_id->at(j)[k] == -999) continue;
              //apply IP3D track selection
              if (which == 1){
              	if (jet_trk_pt->at(j)[k]*1e-3 < 1) continue;
              	if (fabs(jet_trk_ip3d_z0->at(j)[k]*sin(jet_trk_theta->at(j)[k])) > 1.5) continue;
              	if (fabs(jet_trk_ip3d_d0->at(j)[k]) > 1) continue;
              	if ((jet_trk_nPixHits->at(j)[k]+jet_trk_nSCTHits->at(j)[k]) < 7) continue;
              	if ((jet_trk_nPixHoles->at(j)[k]+jet_trk_nSCTHoles->at(j)[k]) > 2) continue;
              	if (jet_trk_nPixHoles->at(j)[k] > 1) continue;
              }
          	  if (which == 2){
          	  	if (jet_trk_pt->at(j)[k]*1e-3 < 0.7) continue;
          	  	if (jet_trk_nSCTHits->at(j)[k] < 4) continue;
          	  	if (jet_trk_nPixHits->at(j)[k] < 1) continue;   	
          	  }

              float etaT = jet_trk_eta->at(j)[k];
              float phiT = jet_trk_phi->at(j)[k];
              float etaJ = jet_eta->at(j);
              float phiJ = jet_phi->at(j);
              float deta = fabs(etaT-etaJ);
              float dphi = fabs(phiT-phiJ) < TMath::Pi() ? fabs(phiT-phiJ) : 2*TMath::Pi() - fabs(phiT-phiJ);
              float dR = sqrt(pow(dphi,2)+pow(deta,2));

              if (jet_trk_orig->at(j)[k] == 0) {//FROMB
              	tot_num = tot_num + 1;
              	b_num = b_num + 1;
              	b_num_no_c = b_num_no_c + 1;
                //reco_d0B->Fill(jet_trk_d0->at(j)[k],weight);
                //reco_z0B->Fill(jet_trk_z0->at(j)[k],weight);
                if (dR > dRB) dRB = dR;
              }
              if (jet_trk_orig->at(j)[k] == 1) {//FROMC
              	tot_num = tot_num + 1;
              	b_num = b_num + 1;
                //reco_d0B->Fill(jet_trk_d0->at(j)[k],weight);
                //reco_z0B->Fill(jet_trk_z0->at(j)[k],weight);
                if (dR > dRB) dRB = dR;
              }
              if (jet_trk_orig->at(j)[k] == 2) {//FRAG
                //reco_d0F->Fill(jet_trk_d0->at(j)[k],weight);
                //reco_z0F->Fill(jet_trk_z0->at(j)[k],weight);
              	tot_num = tot_num + 1;
              	f_num = f_num + 1;
              }
            }

            /*h = (TH1F*)numTot->At((int)(TMath::Log(jetPt/pt_min)/incre));
            h->Fill(tot_num,weight);
            h = (TH1F*)numB->At((int)(TMath::Log(jetPt/pt_min)/incre));
            h->Fill(b_num,weight);
            h = (TH1F*)numF->At((int)(TMath::Log(jetPt/pt_min)/incre));
            h->Fill(f_num,weight);
            h = (TH1F*)numBnoC->At((int)(TMath::Log(jetPt/pt_min)/incre));
            h->Fill(b_num_no_c,weight);
            h = (TH1F*)farDRB->At((int)(TMath::Log(jetPt/pt_min)/incre));
            h->Fill(dRB,weight);*/
            int pT = (int)(TMath::Log(jetPt/pt_min)/incre);
            numTot[central][pT]->Fill(tot_num,weight);
            numB[central][pT]->Fill(b_num,weight);
            numF[central][pT]->Fill(f_num,weight);
            numBnoC[central][pT]->Fill(b_num_no_c,weight);
            farDRB[central][pT]->Fill(dRB,weight);
         }
      }

   }

   //std::vector<float> numBtrks;
   //std::vector<float> numFtrks;

   TGraphErrors* numBtrks[cet_N]; 
   TGraphErrors* numFtrks[cet_N];   
//B & F distribution Graphs

//c0: centrality graphs for each pT interval
//c1: pT graphs for each centrality
   TCanvas* c0 = new TCanvas("c0","c0",500,500);
    TPad* thePad0 = (TPad*)c0->cd();
   TH1F *h0 = thePad0->DrawFrame(-0.5,0,10.5,0.3);
   h0->SetXTitle("Number of B Tracks");
   h0->SetYTitle("Fraction/1");
   h0->SetTitle("Distribution of Number of B Tracks in B jet");

   TCanvas* c1 = new TCanvas("c1","c1",500,500);
   TPad* thePad1 = (TPad*)c1->cd();
   TH1F *h1 = thePad1->DrawFrame(-0.5,0,10.5,0.3);
   h1->SetXTitle("Number of B Tracks");
   h1->SetYTitle("Fraction/1");
   h1->SetTitle("Distribution of Number of B Tracks in B jet");

   TCanvas* c2 = new TCanvas("c2","c2",500,500);
    TPad* thePad2 = (TPad*)c2->cd();
   TH1F *h2 = thePad2->DrawFrame(-0.5,0,10.5,0.3);
   h2->SetXTitle("Number of Fragmentation Tracks");
   h2->SetYTitle("Fraction/1");
   h2->SetTitle("Distribution of Number of Fragmentation Tracks in B jet");

   TCanvas* c3 = new TCanvas("c3","c3",500,500);
   TPad* thePad3 = (TPad*)c3->cd();
   TH1F *h3 = thePad3->DrawFrame(-0.5,0,10.5,0.3);
   h3->SetXTitle("Number of Fragmentation Tracks");
   h3->SetYTitle("Fraction/1");
   h3->SetTitle("Distribution of Number of Fragmentation Tracks in B jet");

   for (int i = 0; i < cet_N; i ++){
   	numBtrks[i] = hotTGraphErrors(Form("numBtrks_cet_%d",i),"numBtrks",kB,0.8,20,1,3003,0.5);
    numFtrks[i]= hotTGraphErrors(Form("numFtrks_cet_%d",i),"numFtrks",kF,0.8,20,1,3003,0.5);
    c1->cd();
    h1->Draw();
    c3->cd();
    h3->Draw();
   myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.02);
   myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.02);
   myText(0.65,0.81,kBlack,"Tracking Selection: ",0.02);
   myText(0.65,0.79,kBlack,Form("%s",track_selection),0.02);
   myText(0.65,0.76,kBlack,Form("%d %% to %d %%", 10*cet[2*i], 10*cet[2*i+1]),0.02);
   myText(0.65,0.73,kBlack,"mcprob > 0.75",0.02);
    for (int j = 0; j < pt_bin; j++){
    	c1->cd();
   		numB[i][j]->Scale(1./numB[i][j]->GetSumOfWeights());
   		numB[i][j]->Draw("SAME");
      numB[i][j]->Write();
   		myBoxText(0.4, 0.87-0.028*j, 0.03, myColor[j],0,Form("%d GeV < p_{T}^{%s} < %d GeV", (int)ptbins[j], TorR, (int)ptbins[j+1]),0.3,myColor[j],21,true,0.02);
   		numBtrks[i]->SetPoint(j,(ptbins[j]+ptbins[j+1])/2.,numB[i][j]->GetMean());
   		numBtrks[i]->SetPointError(j,(ptbins[j+1]-ptbins[j])/2.,numB[i][j]->GetStdDev());
   		c3->cd();
   		numF[i][j]->Scale(1./numF[i][j]->GetSumOfWeights());
   		numF[i][j]->Draw("SAME");
      numF[i][j]->Write();
   		myBoxText(0.4, 0.87-0.028*j, 0.03, myColor[j],0,Form("%d GeV < p_{T}^{%s} < %d GeV", (int)ptbins[j], TorR, (int)ptbins[j+1]),0.3,myColor[j],21,true,0.02);
   		numFtrks[i]->SetPoint(j,(ptbins[j]+ptbins[j+1])/2.,numF[i][j]->GetMean());
   		numFtrks[i]->SetPointError(j,(ptbins[j+1]-ptbins[j])/2.,numF[i][j]->GetStdDev());
    }
    c1->SaveAs(Form("Distribution of Number of B Tracks in B jet PbPb %s rapidity %.1f %s%s%s %d %% %d %%.pdf", dataType, eta_selection, uniqueness, jetTruthness, track_selection, 10*cet[2*i],10*cet[2*i+1]));
    c3->SaveAs(Form("Distribution of Number of Fragmentation Tracks in B jet PbPb %s rapidity %.1f %s%s%s %d %% %d %%.pdf", dataType, eta_selection, uniqueness, jetTruthness, track_selection, 10*cet[2*i],10*cet[2*i+1]));
   }

	for (int j = 0; j < pt_bin; j++){
		c0->cd();
		h0->Draw();
		myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.02);
   		myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.02);
   		myText(0.65,0.81,kBlack,"Tracking Selection: ",0.02);
   		myText(0.65,0.79,kBlack,Form("%s",track_selection),0.02);
   		myText(0.65,0.76,kBlack,Form("%d GeV < p_{T}^{%s} < %d GeV", (int)ptbins[j], TorR, (int)ptbins[j+1]),0.02);
      myText(0.65,0.73,kBlack,"mcprob > 0.75",0.02);
   		c2->cd();
   		h2->Draw();
   		myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.02);
   		myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.02);
   		myText(0.65,0.81,kBlack,"Tracking Selection: ",0.02);
   		myText(0.65,0.79,kBlack,Form("%s",track_selection),0.02);
   		myText(0.65,0.76,kBlack,Form("%d GeV < p_{T}^{%s} < %d GeV", (int)ptbins[j], TorR, (int)ptbins[j+1]),0.02);
      myText(0.65,0.73,kBlack,"mcprob > 0.75",0.02);
		for (int i = 0; i < cet_N; i++){
			c0->cd();
			numB[i][j]->SetLineColor(myColor[j]+2*i);
			numB[i][j]->SetMarkerColor(myColor[j]+2*i);
			numB[i][j]->Draw("SAME");
			myBoxText(0.4, 0.87-0.028*i, 0.03, myColor[j]+2*i,0,Form("%d %% to %d %%", 10*cet[2*i], 10*cet[2*i+1]),0.3,myColor[j]+2*i,21,true,0.02);
			c2->cd();
			numF[i][j]->SetLineColor(myColor[j]+2*i);
			numF[i][j]->SetMarkerColor(myColor[j]+2*i);
			numF[i][j]->Draw("SAME");
			myBoxText(0.4, 0.87-0.028*i, 0.03, myColor[j]+2*i,0,Form("%d %% to %d %%", 10*cet[2*i], 10*cet[2*i+1]),0.3,myColor[j]+2*i,21,true,0.02);
		}
		c0->SaveAs(Form("Distribution of Number of B Tracks in B jet PbPb %s rapidity %.1f %s%s%s %d to %d.pdf", dataType, eta_selection, uniqueness, jetTruthness, track_selection, (int)ptbins[j], (int)ptbins[j+1]));
		c2->SaveAs(Form("Distribution of Number of Fragmentation Tracks in B jet PbPb %s rapidity %.1f %s%s%s %d to %d.pdf", dataType, eta_selection, uniqueness, jetTruthness, track_selection, (int)ptbins[j], (int)ptbins[j+1]));
	}

   TCanvas* c5 = new TCanvas("c5","c5");
   TPad* thePad5 = (TPad*)c5->cd();
   TH1F *h5 = thePad5->DrawFrame(-0.02,0,0.4,0.2);
   h5->SetXTitle("#DeltaR_{max}");
   h5->SetYTitle(Form("Fraction/%.3f",1./dist_bin3));
   h5->SetTitle("Distribution of Maximum #DeltaR between B tracks and jet in B jet");
   

for (int i = 0; i < cet_N; i++){
  h5->Draw();
   myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.02);
   myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.02);
   myText(0.65,0.81,kBlack,"Tracking Selection: ",0.02);
   myText(0.65,0.79,kBlack,Form("%s",track_selection),0.02);
   myText(0.65,0.76,kBlack,Form("%d %% to %d %%", 10*cet[2*i], 10*cet[2*i+1]),0.02);

   for (int j = 0; j < pt_bin; j++){
    c5->cd();
    farDRB[i][j]->Scale(1./farDRB[i][j]->GetSumOfWeights());
    farDRB[i][j]->Draw("SAMELP");
    myBoxText(0.4, 0.87-0.028*j, 0.03, myColor[j],0,Form("%d GeV < Jet p_{T}^{%s} < %d GeV", (int)ptbins[j], TorR, (int)ptbins[j+1]),0.3,myColor[j],21,true,0.02);
   }

   c5->SaveAs(Form("Distribution of Maximum dR of B Tracks in B jet PbPb %s rapidity %.1f %s%s%s centrality %d %% %d %%.pdf", dataType, eta_selection, uniqueness, jetTruthness, track_selection, 10*cet[2*i], 10*cet[2*i+1]));
 }

//number of B vs F
	TCanvas* c6 = new TCanvas("c6","c6",500,500);
  for (int i = 0; i < cet_N; i++){
	TMultiGraph* g = hotTMultiGraph("Average Number of Tracks in b-jet",Form("p_{T}^{%s, jet}",TorR),"Number of Tracks",50,600,0,13);
	g->Add(numBtrks[i]);
	numBtrks[i]->Write();
	g->Add(numFtrks[i]);
	numFtrks[i]->Write();
	gPad->SetTicks(1);
	g->Draw("AP2");
	myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.03);
  myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.03);
  myText(0.65,0.81,kBlack,"Tracking Selection: ",0.03);
  myText(0.65,0.78,kBlack,Form("%s",track_selection),0.03);
  myText(0.65,0.75,kBlack,"mcprob > 0.75",0.03);
  myText(0.65,0.72,kBlack,Form("PbPb: %d %% to %d %%",10*cet[2*i],10*cet[2*i+1]),0.03);
  myText(0.65,0.65,kBlack,Form("%s",dataType),0.03);
  myBoxText(0.2, 0.87, 0.04, kB,0.5,"Number of B Tracks",0.3,kB,20,false,0.03,3003);
	myBoxText(0.2, 0.82, 0.04, kF,0.5,"Number of Fragmentation Tracks",0.3,kF,20,false,0.03,3003);
  c6->SaveAs(Form("Number of Tracks in B jet pp %s rapidity %.1f %s%s%s %d %% to %d %%.pdf", dataType, eta_selection, uniqueness, jetTruthness,track_selection,10*cet[2*i],10*cet[2*i+1]));
}
	


/*TCanvas* c4 = new TCanvas("c4","c4",500,500);
   TPad* thePad4 = (TPad*)c4->cd();
   TH1F *h4 = thePad4->DrawFrame(-highd0,1e-4,highd0,1.0);
   h4->SetXTitle("d0");
   h4->SetYTitle(Form("Fraction/%.2f",2.*highd0/perigeebins));
   h4->SetTitle("Distribution of d0 in Tracks in B jet");
   h4->Draw();
   reco_d0B->Scale(1./reco_z0B->GetSumOfWeights());
   reco_d0F->Scale(1./reco_z0F->GetSumOfWeights());
   reco_d0B->Draw("SAME");
   reco_d0F->Draw("SAME");
    myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.03);
    myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.03);
    myText(0.65,0.81,kBlack,"Tracking Selection: ",0.03);
    myText(0.65,0.78,kBlack,Form("%s",track_selection),0.03);
    myText(0.65,0.75,kBlack,"mcprob > 0.75",0.03);
    myText(0.65,0.70,kBlack,Form("%s",dataType),0.03);
    myBoxText(0.2, 0.87, 0.04, kB, 0, "B Tracks",0,kB,20,true,0.03);
    myBoxText(0.2, 0.82, 0.04, kF, 0, "Fragmentation Tracks",0,kF,20,true,0.03);
    c4->SetLogy();
    c4->SaveAs(Form("d0 of Tracks in B jet pp %s rapidity %.1f %s%s%s.pdf", dataType, eta_selection, uniqueness, jetTruthness,track_selection));

  TCanvas* c5 = new TCanvas("c5","c5",500,500);
   TPad* thePad5 = (TPad*)c5->cd();
   TH1F *h5 = thePad5->DrawFrame(-highz0,1e-4,highz0,1.0);
   h5->SetXTitle("z0");
   h5->SetYTitle(Form("Fraction/%.2f",2.*highz0/perigeebins));
   h5->SetTitle("Distribution of z0 in Tracks in B jet");
   h5->Draw();
   reco_z0B->Scale(1./reco_z0B->GetSumOfWeights());
   reco_z0F->Scale(1./reco_z0F->GetSumOfWeights());
   reco_z0B->Draw("SAME");
   reco_z0F->Draw("SAME");
    myText(0.65,0.87,kBlack,Form("#eta < |%.1f|",eta_selection),0.03);
    myText(0.65,0.84,kBlack,"pp b#bar{b} filtered MC",0.03);
    myText(0.65,0.81,kBlack,"Tracking Selection: ",0.03);
    myText(0.65,0.78,kBlack,Form("%s",track_selection),0.03);
    myText(0.65,0.75,kBlack,"mcprob > 0.75",0.03);
    myText(0.65,0.70,kBlack,Form("%s",dataType),0.03);
    myBoxText(0.2, 0.87, 0.04, kB, 0, "B Tracks",0.3,kB,20,true,0.03);
    myBoxText(0.2, 0.82, 0.04, kF, 0, "Fragmentation Tracks",0.3,kF,20,true,0.03);
    c5->SetLogy();
    c5->SaveAs(Form("z0 of Tracks in B jet pp %s rapidity %.1f %s%s%s.pdf", dataType, eta_selection, uniqueness, jetTruthness,track_selection));

   numB->Write();
   numF->Write();
   numBnoC->Write();
   numTot->Write();

   reco_d0B->Write();
   reco_d0F->Write();
   reco_d0sig->Write();
   reco_z0B->Write();
   reco_z0F->Write();
   reco_z0sig->Write();*/

   out->Close();
}
