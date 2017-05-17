#include "TPaletteAxis.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

void PlotComparaison(){

  TFile* fCheck = new TFile("/home/kswiss/Workspace/workoffline/CompStudy/saved/170509/PreCheckTracesColumn10.3.root");
  TFile* fSimulation = new TFile("/home/kswiss/Workspace/workoffline/CompStudy/saved/170509/tracesSimulation.3.root");
  TFile* fData = new TFile("/home/kswiss/Workspace/workoffline/CompStudy/saved/170509/tracesData.3.root");
  
  vector<TH1F*> hCheck;
  vector<TH1F*> hSimulation;
  vector<TH1F*> hData;

  for(int i = 0; i < 22; i++){
    TString hName("h"); hName+=1;hName+="_m";hName+=4;hName+="_r";hName+=i;hName+="_c";hName+=10;
    if(fCheck->Get(hName)) hCheck.push_back((TH1F*)fCheck->Get(hName));
    if(fData->Get(hName)) hData.push_back((TH1F*)fData->Get(hName));
    if(fSimulation->Get(hName)) hSimulation.push_back((TH1F*)fSimulation->Get(hName));

  }

  cout << "Number of Pixels in Column 10 for Data, Sim, Check: " <<  hData.size() << " " << hSimulation.size() << " " <<  hCheck.size() << endl;

  TCanvas* cTraces = new TCanvas("cTraces","cTraces",800,600);
  THStack* sc = new THStack("sc","");
  TLegend* lc = new TLegend(.15,.35,.45,.85);

  TString sR12C10("h1_m4_r12_c10"),sR10C10("h1_m4_r10_c10"),sR8C10("h1_m4_r8_c10"), sR7C10("h1_m4_r7_c10"), sR6C10("h1_m4_r6_c10");

  for(int i = 0; i < hCheck.size(); i++) {
    hCheck[i]->SetLineWidth(2);
    hCheck[i]->SetLineColor(kBlack);
    //   if(hCheck[i]->GetName() == sR8C10 || hCheck[i]->GetName() == sR10C10 || hCheck[i]->GetName() == sR12C10) sc->Add(hCheck[i]);    
    //     if(hCheck[i]->GetName() == sR10C10 ) sc->Add(hCheck[i]);    
    if(hCheck[i]->GetName() == sR10C10 || hCheck[i]->GetName() == sR8C10) {
      sc->Add(hCheck[i]);
      TString lcName = hCheck[i]->GetName();
      lc->AddEntry(hCheck[i],lcName.Replace(0,5,"check_"),"l");
    }
  }
  for(int i = 0; i < hData.size(); i++) {
    hData[i]->SetLineWidth(2);
    hData[i]->Scale(5);//ADC to Photon Count
    hData[i]->Rebin(5);
    hData[i]->SetLineColor(kRed);
    if(hData[i]->GetName() == sR10C10 || hData[i]->GetName() == sR8C10) {
      sc->Add(hData[i]);
      // lc->AddEntry(hData[i],hData[i]->GetName(),"l");
      TString lcName = hData[i]->GetName();
      lc->AddEntry(hData[i],lcName.Replace(0,5,"data_"),"l");
     }
    //    if(hData[i]->GetName() == sR10C10) sc->Add(hData[i]);    
  }
  for(int i = 0; i < hSimulation.size(); i++) {
    hSimulation[i]->SetLineWidth(2);
    hSimulation[i]->Scale(5);//ADC to Photon Count
    hSimulation[i]->Rebin(5);
    hSimulation[i]->SetLineColor(kBlue);
    // if(hSimulation[i]->GetName() == sR8C10 || hSimulation[i]->GetName() == sR10C10 || hSimulation[i]->GetName() == sR12C10) sc->Add(hSimulation[i]);    
    //    if(hSimulation[i]->GetName() == sR10C10) sc->Add(hSimulation[i]);    
    if(hSimulation[i]->GetName() == sR10C10 ||hSimulation[i]->GetName() == sR8C10) {
      sc->Add(hSimulation[i]);
      //      lc->AddEntry(hSimulation[i],hSimulation[i]->GetName(),"l");
      TString lcName = hSimulation[i]->GetName();
      lc->AddEntry(hSimulation[i],lcName.Replace(0,5,"sim_"),"l");
       
    }
  }
  
  sc->Draw("nostack");
  lc->Draw();
  cTraces->SaveAs("tracesmerged.png");
  
}
