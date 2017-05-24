#include "ComparaisonStudy.h"
#include <evt/Event.h>
#include <fevt/FEvent.h>
#include <fevt/Eye.h>
#include <fevt/TelescopeSimData.h>
#include <fevt/Telescope.h>
#include <fevt/PixelSimData.h>
#include <fevt/Pixel.h>
#include <fevt/FdConstants.h>
#include <det/Detector.h>
#include <fdet/FDetector.h>
#include <fdet/Eye.h>
#include <fdet/Telescope.h>
#include <fdet/Camera.h>
#include <fdet/Pixel.h>
#include <fdet/Channel.h>
#include <fdet/Mirror.h>
#include <fdet/Filter.h>
#include <fdet/Corrector.h>

#include <utl/Point.h>
#include <utl/UTMPoint.h>
#include <utl/Particle.h>
#include <utl/Vector.h>
#include <utl/TimeStamp.h>
#include <utl/Particle.h>
#include <utl/ErrorLogger.h>
#include <utl/ReferenceEllipsoid.h>
#include <utl/RandomEngine.h>
#include <fwk/LocalCoordinateSystem.h>
#include <fwk/CoordinateSystemRegistry.h>
#include <fwk/CentralConfig.h>
#include <fwk/SVNGlobalRevision.h>
#include <fwk/RandomEngineRegistry.h>
#include <utl/Math.h>

#include <utl/MathConstants.h>
#include <utl/CoordinateSystemPtr.h>
#include <utl/AxialVector.h>
#include <utl/Vector.h>
#include <utl/PhysicalConstants.h>
#include <utl/AugerUnits.h>
#include <utl/UTCDateTime.h>
#include <CLHEP/Random/Randomize.h>

#include <TMap.h>
#include <TTree.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TH2D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TF1.h>
#include <TFile.h>
#include <TVector3.h>
#include <vector>
#include <AugerEvent.h>
#include <EyeEvent.hh>
#include <EyePixelList.hh>
#include <FadcData.hh>
#include <TPaletteAxis.h>

using namespace std;
using namespace utl;
using namespace fwk;
using namespace evt;
using namespace fevt;
using namespace fdet;
using namespace det;
using CLHEP::RandFlat;

namespace ComparaisonStudyNS {

// Initialize your static data members.

// Define the methods preferibly in the same order as listed in the header file.

ComparaisonStudy::ComparaisonStudy() 
{
  // Initialize your class.
  // Prefer initialization lists over = to 
  // initialize your members.
}

ComparaisonStudy::~ComparaisonStudy()
{
}

VModule::ResultFlag 
ComparaisonStudy::Init()
{
  INFO("Init()");
  //here we want to start by comparing only the first page of the traces. later we will do a full comparaison when the simulation will be doing followers properly.   
  fEventCounter = 0;
  
  outputPlots = new TFile("FullTraces.root","recreate");

  for(int i = 0; i < fNRows; i++){  
    TString hPixelRowNameSim("hPixelRowSim");hPixelRowNameSim+= (i+1);
    TString hPixelRowNameDat("hPixelRowDat");hPixelRowNameDat+= (i+1);
    TString hPixelRowTitle("Row ");hPixelRowTitle+= (i+1);hPixelRowTitle+= ";Time Bin (100 ns);Column Number;ADC Counts / 20 Time Bins";
    hPixelRow[0][i] = new TH2F(hPixelRowNameSim,hPixelRowTitle,50,0,999,fNColumns,0,fNColumns);
    hPixelRow[1][i] = new TH2F(hPixelRowNameDat,hPixelRowTitle,50,0,999,fNColumns,0,fNColumns);
  }

  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);
  //  detTelGlobal = detTelTMP;
  for(int i = 0; i < fNTels; i++){
    for(int j = 0; j < fNPixels; j++){

      //intialize storing for sim
      TString hnameSim("hSim_tel");hnameSim+=i+1;hnameSim+="_r";hnameSim+=detTelGlobal.GetPixel(j+1).GetRow();hnameSim+="_c";hnameSim+=detTelGlobal.GetPixel(j+1).GetColumn();
      fhRawPixel[0][i][j] = new TH1F(hnameSim,hnameSim,50,0,999);
      //and now for data
      TString hnameDat("hDat_tel");hnameDat+=i+1;hnameDat+="_r";hnameDat+=detTelGlobal.GetPixel(j+1).GetRow();hnameDat+="_c";hnameDat+=detTelGlobal.GetPixel(j+1).GetColumn();
      fhRawPixel[1][i][j] = new TH1F(hnameDat,hnameDat,50,0,999);


    }
  }
  return eSuccess;
}

VModule::ResultFlag 
ComparaisonStudy::Run(evt::Event& event)
{
  const evt::Header& header = event.GetHeader();  
  
  AugerEvent& rawEvent = event.GetRawEvent();
  fEventCounter++;

  //to deal with sim and data in input file
  if(fEventCounter == 1) {
    fData = false;
    iTag = 0;
    INFO("New SIM Event!");
  }else{
    fData = true;
    iTag = 1;
    if (!(header.GetId() == "eye1_run5254_event6336")) return eSuccess;
    INFO("New DATA Event!");
  }
  
  for (AugerEvent::EyeIterator eyeIter = rawEvent.EyesBegin();
       eyeIter != rawEvent.EyesEnd(); ++eyeIter) {
    
    TEyeEvent& eyeEvent = *eyeIter;
    TEyePixelList* eyePixelList = eyeEvent.GetPixelList();
    TEyeFADCData* eyeFADCData = eyeEvent.GetFADCData();
    const unsigned int numpixels = eyePixelList->GetNumPixels();
    TEyeEventHeader *eyeHeader = eyeEvent.GetEventHeader();
    const int eyeId = eyeEvent.GetEventHeader()->GetEyeNo();
    const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
    FEvent& fdEvent = event.GetFEvent();
    fevt::Eye& eye = fdEvent.GetEye(eyeId);
    const fdet::Eye& detEye = detFD.GetEye(eyeId);


    cout << eyeEvent.GetEventHeader()->GetEventNo() << " " << eyePixelList->GetNumPixels() << " " <<  eyePixelList->GetEyeNo() << endl;
    
    TString fName("traces.root");
    TFile *fTMP;
    if( fEventCounter ==1) {
      fTMP = new TFile(fName,"recreate");
    }else{
      fTMP = new TFile(fName,"update");
    }

    for (unsigned int iPixel = 0; iPixel < eyePixelList->GetNumPixels(); ++iPixel) {
      
      FdUtil::Fd::FdPixelNumber pixelNumber = eyePixelList->GetPixel(iPixel);
      const unsigned int mirrorId = FdUtil::Fd::GetEyeMirrorNo(pixelNumber);
      fevt::Telescope& tel = eye.GetTelescope(mirrorId);
      //cout << " detEye " << eyeId<<endl;
      const fdet::Telescope& detTel = detEye.GetTelescope(mirrorId);
      const unsigned int channelId = FdUtil::Fd::GetMirrorPixelNo(pixelNumber);
      //     const unsigned int pixelId = FdUtil::Fd::GetPixelNo(pixelNumber);
      
      const TFADCData* fadcData = eyeFADCData->GetFADCData(pixelNumber);

      
      const fdet::Channel& thisChannel = detTel.GetChannel(channelId);
      const unsigned int pixelId = thisChannel.GetPixelId();


      if (!fadcData || thisChannel.IsVirtual() || pixelId > 440) {
	continue;
      }
      //cout << mirrorId << " " << channelId << " " << pixelId << endl;
      
      
      const unsigned int startBin = fadcData->GetTraceStartBin();
      const unsigned int endBin = fadcData->GetTraceEndBin();
      
      TFADCData::FADCDataWord fadcword = fadcData->GetFADCTrace();
      
      
      baseline=0.;
      baselineRMS=0.;
      
      for (unsigned int pos = startBin; pos <startBin+280; ++pos) baseline += int(FADCDataWordGetData(&fadcword[pos]))/280.;
      
      
      for (unsigned int pos = startBin; pos <startBin+280; ++pos) baselineRMS += pow(int(FADCDataWordGetData(&fadcword[pos])),2)/280.;
      
      baselineRMS -= pow(baseline,2);
      if (baselineRMS>0) baselineRMS=sqrt(baselineRMS);
      else baselineRMS=0.;


      TString hName("h");
      if(fData) hName+="Data";
      if(!fData) hName+="Sim";
      hName+="_m";hName+=mirrorId;hName+="_r";hName+=detTel.GetPixel(pixelId).GetRow();hName+="_c";hName+=detTel.GetPixel(pixelId).GetColumn();
      cout << hName.Data() << endl;
      TH1F* trace = new TH1F(hName,hName,50,startBin,endBin);
      
      for (unsigned int pos = startBin; pos <= endBin; ++pos) {
	charge = FADCDataWordGetData(&fadcword[pos]);
	//	cout << int(pos/20.) << endl;
	trace->SetBinContent(int(pos/20.),charge-baseline);
	  
      }//trace loop

      cout << iTag << " " << trace << " " << mirrorId << " " << pixelId << endl;
      GlueTrace(iTag,trace, mirrorId, pixelId);//here glue means save in memory. 
    }//pixel loop
    fTMP->Write();
    fTMP->Close();

    //    return eSuccess;
      
    }//mirror loop


  //need to fillall arrays before moving on 
  if(fEventCounter == 1) return eSuccess;

  const TimeStamp times = UTCDateTime(2017,05,19,0,0).GetTimeStamp();
  Detector::GetInstance().Update(times);
  const fdet::FDetector& detFD = det::Detector::GetInstance().GetFDetector();
  const fdet::Eye& detEye = detFD.GetEye(4);
  const fdet::Telescope& detTelGlobal = detEye.GetTelescope(4);

  
  for(int k = 0; k < 2; k++){
    for(int i = 0; i < fNTels; i++){
      for(int j = 0; j < fNPixels; j++){
	if(fhRawPixel[k][i][j]->GetEntries() == 0){
	  delete fhRawPixel[k][i][j];
	  continue;
	}
	for(int iBin = 1; iBin <= fhRawPixel[k][i][j]->GetNbinsX(); iBin++){ 
	  hPixelRow[k][detTelGlobal.GetPixel(j+1).GetRow()-1]->SetBinContent(iBin,detTelGlobal.GetPixel(j+1).GetColumn(),fhRawPixel[k][i][j]->GetBinContent(iBin));
	} 
      }
    } 
  }  

  TH2F *hResidual[fNRows];
  for(int iRow = 0; iRow < fNRows; iRow++){
    TString hResidualName("hResidual");hResidualName+=iRow+1;
    TString hResidualTitle("Residual for Row ");hResidualTitle+=iRow+1;
    //get the data histogram

    hResidual[iRow] = (TH2F*) hPixelRow[1][iRow]->Clone();
    hResidual[iRow]->SetNameTitle(hResidualName, hResidualTitle);
    hResidual[iRow]->Add(hPixelRow[0][iRow],-1.);
  }

  //plotting the residual for all selected rows in 2D hist form
  TCanvas* cglues = new TCanvas("cglues","cglues",1000,800);
  gStyle->SetOptStat(0);
  for(int i=4; i <= 14; i++){
    hResidual[i]->Draw("colz");
    hResidual[i]->GetZaxis()->SetTitleOffset(1.02); 
    hResidual[i]->GetZaxis()->SetTitle("Residual ADC/2#mus"); 
    cglues->Update();
    TPaletteAxis *palette =(TPaletteAxis*)hResidual[i]->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.9);
    palette->SetX2NDC(0.92);
    palette->SetY1NDC(0.1);
    palette->SetY2NDC(0.9);
    cglues->Modified();
    
    TString fNameGlues("glues/");fNameGlues+=hResidual[i]->GetName();
    fNameGlues+=".png";
    cglues->SaveAs(fNameGlues.Data());  
  }

  Color_t colorList[11] = {kPink,kRed,kOrange,kSpring+8,kTeal+8,kTeal,kCyan,kAzure-4,kBlue,kViolet+10,kBlack};
  TCanvas* cresidual = new TCanvas("cresidual","cresidual",1000,800);
  TMultiGraph *mg = new TMultiGraph();
  //TH2F * hResidualCanvas = new TH2F("hResidualCanvas",";Event Number;Residual",100,0,12,100,-500,500);
  // hResidualCanvas->Draw();
  for(int i=6; i <= 10; i++){  //look through rows 6-10
    double residualTotal;
    TGraph* gResidual = new TGraph(1);
    for(int iBinx=1;iBinx<=hResidual[i]->GetNbinsX(); iBinx++){
      residualTotal += hResidual[i]->GetBinContent(iBinx, 10);//look at column 10
    }
    residualTotal = residualTotal/hResidual[i]->GetNbinsX();
    cout << fEventCounter << " " << residualTotal << endl;
    //fEventCounter wont work. Need to think on how to implement multiple sim filesx
    gResidual->SetPoint(0,1,residualTotal);
    gResidual->SetMarkerStyle(21);
    gResidual->SetMarkerColor(colorList[i-5]);
    gResidual->SetMarkerSize(2.0);
    //    gResidual.Draw("P same");
    mg->Add(gResidual);
      //    cresidual->Update();
  }
  mg->Draw("AP");
  cresidual->SaveAs("gResidual.png");
  //   hPixelRow[0][i-1]->Draw("colz");
  //   hPixelRow[0][i-1]->GetZaxis()->SetTitleOffset(0.98);
  //   cglues->Update();
  //   TPaletteAxis *palette =
  //     (TPaletteAxis*)hPixelRow[0][i-1]->GetListOfFunctions()->FindObject("palette");
  //   palette->SetX1NDC(0.9);
  //   palette->SetX2NDC(0.92);
  //   palette->SetY1NDC(0.1);
  //   palette->SetY2NDC(0.9);
  //   // palette->SetX1NDC(0.1);
  //   // palette->SetX2NDC(0.12);
  //   // palette->SetY1NDC(0.1);
  //   // palette->SetY2NDC(0.5);
  //   cglues->Modified();
  //   TString fNameGlues("glues/");fNameGlues+=hPixelRow[0][i-1]->GetName();fNameGlues+=".png";
  //   cglues->SaveAs(fNameGlues.Data());
  // }

  
  outputPlots->Write();
  outputPlots->Close();

  
  return eSuccess;
}

VModule::ResultFlag 
ComparaisonStudy::Finish() 
{
  
  return eSuccess;
}

VModule::ResultFlag
ComparaisonStudy::GlueTrace(int iTagtmp, TH1F* theTrace, int telIdtmp, int pixelIdtmp) 
{
  cout << theTrace->GetNbinsX() << " " << fhRawPixel[iTagtmp][telIdtmp-1][pixelIdtmp-1]->GetNbinsX() <<endl;
  for(int iBin = 1; iBin <= theTrace->GetNbinsX(); iBin++){
    double normalizedBinContent = theTrace->GetBinContent(iBin);
    fhRawPixel[iTagtmp][telIdtmp-1][pixelIdtmp-1]->SetBinContent(iBin,normalizedBinContent);
  }
  return eSuccess;
}




}
