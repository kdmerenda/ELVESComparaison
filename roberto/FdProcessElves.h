/**
   \file
   declaration of FdProcessElves

   \author 
   \version $Id: FdProcessElves.h 15247 2017-01-17  $
 \\ modified by rmussa 17/1/2017
   
*/


#ifndef _FdProcessElves_
#define _FdProcessElves_

namespace utl {class AxialVector;}
namespace fevt{
  class Eye;
  class Pixel;
}

//static const char CVSId_FdSDPFinderOG_FdSDPFinder[] = 
//"$Id: FdSDPFinder.h 15247 2009-11-24 19:12:49Z rulrich $";

#include <fwk/VModule.h>
#include <utl/Vector.h>


//namespace FdProcessElves {

  /**
     \class FdProcessElves FdProcessElves.h "FdProcessElves/FdProcessElves.h"

     \brief ELVES Reconstruction code 

     \author Roberto Mussa
     \date 21 Sep 2016
     \ingroup FDRecModules
  */


  class FdProcessElves :public fwk::VModule{

  public:
   
    fwk::VModule::ResultFlag Init();
    fwk::VModule::ResultFlag Run(evt::Event& event);
    fwk::VModule::ResultFlag Finish();

    //    std::string GetSVNId() const 
    //  {return std::string("$Id: FdSDPFinder.h 15247 2009-11-24 19:12:49Z rulrich $");} 

  private:
      TObjArray *  fPixels;
      TObjArray * fPixels2;
      TObjArray * fPixels3;
      TObjArray * fPixels4;
      TObjArray * fPixels5;
      TObjArray * fPixels6;
      TObjArray * fPixels7;
      TObjArray * fPixels8;
      TObjArray * fSites;
      TObjArray * fSites2;
      
    

    std::string    fEventId_for_profile_plot;
    std::string    fEye_number_string;

      
    //utl::AxialVector    FindSDPFirstGuess(fevt::Eye& eye);
    //void RemoveNoise(fevt::Eye& eye);
    //void ReadmitPixel(fevt::Eye& eye);
    //void FindSDP(fevt::Eye& eye);
    //void EstimateChi_i(fevt::Eye& eye);
    //void RemoveOutliers(fevt::Eye& eye);
    bool IsIsolated(const fevt::Pixel&, const fevt::Eye& eye ) ;

    static void  MinuitFitFunc(int& npar, double *gin, double& f,
			       double *par, int iflag);
      
    void ProcessRawData(TEyeEvent& eyeEvent, fevt::Eye& eye);
      void ProcessHistos();
      int FdAllTraces_CameraViewPlot();
      void FdElvesGeometryCorrection();
      void FdElvesAtmosphericCorrection();
      void FdSingleElvesPixel_CameraViewIni();
      void FdSingleElvesPixel_LongLatViewIni(double);

      void FdSingleElvesTime_CameraViewPlot(int);
      void FdSingleElvesPull_CameraViewPlot(int);
      void FdSingleElvesSigmaT_CameraViewPlot(int);
      void FdSingleElvesSigmaR_CameraViewPlot(int);
      void FdSingleElvesQmax_CameraViewPlot(int);
      void FdSingleElvesQCtot_CameraViewPlot(int);
      void FdSingleElvesQCCtot_CameraViewPlot(int);
      
      void FdSingleElvesQCCtot_LongLatViewPlot(TCanvas* );
      void FdSingleElvesQCCtot_LongLatViewPlot2(TCanvas* );
      
      void GlueTraces(TEyeEvent& eyeEvent, fevt::Eye& eye);
      bool WrongEventLabel(int EventLabel);
      bool NewEventFound(fevt::Eye& eye);
      int NbaysEventFoundinPreScan(fevt::Eye& eye);
      void FdElvesPulseFinder(int it, int ip);
      bool FindPulse(int it, int ip);
      int FindSecondPulse(int it, int ip);
      TVector3 FdSingleElvesBoltFinder();
      TVector3 FdDoubleElvesBoltFinder();
      TVector3 CalculateRawBoltLoc(int, int, int);
      int FdSingleElvesBoltFitter(TVector3 , double  );
      int CheckTimeFit5D_Output( TF1* );
      int PlotTimeFit5D_Output( TF1* );
      //TVector3 FdSingleElvesBoltFitter();
      void DrawSkewedDipole(double, double );
          

      void ParametrizePixelPulse(int it,int ip);
      void ParametrizePixelPulse2(int it,int ip, int ic);
      void InitializeRowMarkers(int );
      void InitializeColMarkers(int );   // added 03/04/2017
      double AngDistFromFirstPixel(double);
      double CalculateGeomCorr( int , int , int );
      double CalculateAtmoCorr( int );
      void  LightEmission_vs_AnglePlot();
      void  LightEmission_vs_DistancePlot();  // added mid March 2017
      void  LightEmission_vs_DistanceFit(); // added 14 Apr 2017

      int  LightEmission_vs_EllipticDistancePlot(); // added 11 Apr 2017
      
      void  FdElvesBestEllipseSearch(int );
      double PathFromBolt_to_FD(int , Double_t *);   // Replica of function  t0LinFitter
      double CalculateEllipticDistance(int , double, double);   // Replica of function  t0LinFitter
      void LongLatView_QPLOTS();
      void CameraView_4TPLOTS();
      void CameraView_4QPLOTS();

      double BestFitParameter[9];
      

     // void DrawAmplitudeColor_LongLatMap(TEyeEvent& eyeEvent, fevt::Eye& eye);
     // void DrawAmplitudeColor_CorrectedLongLatMap(TEyeEvent& eyeEvent, fevt::Eye& eye);
      void DrawTotChargeColor_FDview(TEyeEvent& eyeEvent, fevt::Eye& eye);
      void CalculateDelaysBetweenBays(TEyeEvent& eyeEvent,fevt::Eye& eye);
      void LastReset();
      void ResetAllForNewEvent();

      void WriteOutputFile(long unsigned int, long unsigned int);
      void WriteAllHistos();

      void ResetAllCounters();
      void ResetAllVectors();
      void ResetAllHistos();

//utility functions
      //static Double_t ellipse(Double_t *x, Double_t *par);
      static Double_t circumference(Double_t *x, Double_t *par);
      //static Double_t Eq2root(double,double,double,int );
      
      static Double_t asymGaussian(Double_t *x, Double_t *par);
      static Double_t P2exp(Double_t *x, Double_t *par);
      static Double_t t0LinFitter(Double_t *x, Double_t *par);
      //Double_t t0LinFitter(Double_t *x, Double_t *par);

      
      static Double_t MyFunc2D(Double_t *x, Double_t *par);
      void Fit2(const fevt::Eye & eye);
      void Fit3D(const fevt::Eye & eye);
      
      
      static  double RootOf(double a, double b, double c, int i);


      // Output files (rmussa for elves analysis)
      TFile *fOut;
      fstream fOut2;

      
      
    /// minimum number of pixels to perform reconstruction
    unsigned int     fMinPixels;
    /// maximum number of pixels accepted for a good event
    unsigned int     fMaxPixels;

    /// Time limit for pixel isolation 
    double           fIsolationLimit;
    bool             fGoodEvent;
 
      
      
      //fevt::Eye * FdProcessElves::fCurEye =0;
      //double FdProcessElves::fChi2=0;
      //unsigned int FdProcessElves::fNDof=0;
      unsigned int nPages = 0;
      //double pixel_slot[3][NPIXELS][150];
      //double pixPed[3][3][NPIXELS];
      static const int size = 1000;
      static const int nSlotInterval = 20;
      static const int NPIXELS = 440;
      static const int NROWS=22;
      static const int NCOLS =20;
      static const int NMAXPAGES=3;
      static const int NMAXTELS=3;
      static const int NMAXLINES=4000;
      static const int NP2MIN=10;
      static const int NDP1=10;
      static const int NDP2=64;
      static const int NDBbins = 50;
    
      static const int nslots = NMAXPAGES*size/nSlotInterval;
    //  static const int nslots = 450;
      
      static constexpr float DTBIN = 0.1; // bin width in microseconds
      //static const float DTSLOT = DTBIN * nSlotInterval; // bin width in microseconds
      static constexpr float DTSLOT = 2.0; // bin width in microseconds
      
      static constexpr float QMIN = 1500.;  // minimum threshold for a pulse
      static constexpr float DTMIN = 10.;  // minimum time interval between pulses
      
      static constexpr double degree = 3.14159265359/180.;
      static constexpr double km = 1.;
      double m = km/1000.;
      // Display image in canvas and pad.
      //
      double backwall[4]; // = {-30.00,60.03,-171.85,-116.68}; 	// with respect to the Est, in degrees
      double hsite[4];// = {1.45,1.45,1.45,1.7}; // altitude asl
      static constexpr double Rearth = 6371.; // in km
      static constexpr double SpeedOfLight = 0.2997; // in km / usec;
      
      static constexpr double DDELTA_KM = 50.; // in km / usec;

      double Amirror = 11.*m*m;
      double Hd0 =85.*km;
      // FD coordinates (in Easting and Northing)
      
      double longsite[4]; // = {-69.449673,-69.012203,-69.210845,-69.599633};
      double latsite[4]; // = {-35.495759,-35.291974,-34.935916,-35.114138};
      int sitecolor[4];
      
      TGraph *PixelDirs[NMAXTELS];
      TGraph *PixelEdges[NPIXELS*NMAXTELS];
      TGraph *PixelDirsLL[NMAXTELS];
      TGraph *PixelEdgesLL[NPIXELS*NMAXTELS];

      TGraph *fB100km;
      TGraph *fB200km;
      TGraph *fBE1_100km;
      TGraph *fBE1_200km;
      
      //TImage *img;
      
      

      int site1=1; //default in Los Leones
      
      //double Hd=77.*km; // Height of D layer in km
      //double dHd=1.*km; // Height of D layer in km

      int NpsPages;
      int GPSsec[NMAXLINES], GPSns[NMAXLINES], evno[NMAXLINES], nwd[NMAXLINES], bayno[NMAXLINES];
      int dts[NMAXPAGES-1][NMAXLINES];
      
      int ThisBay1;
      int ThisBay2;
      int LeftMostBay=0;
      bool ProcessingEvent = false ;
      int EffNpages=0;
      int PreScanIndex = -1;
      int rawMinGPSns = 1000000000;


      
      int fTelescopeCount, fL7TelescopeCount;
      TVector3 Site1Location = {1.,0.,0.};
      TVector3 East1 = {0.,1.,0.};
      TVector3 PoloNord = {0.,0.,1.};

      double rawPed[3][NPIXELS];
      double rawPed2[3][NPIXELS];
      //double rawTrace[3][NPIXELS][nslots];
      double calTrace[3][NPIXELS][nslots];
      unsigned int EmTime_ns[3][NPIXELS][nslots];
      double LatPixel[3][NPIXELS];
      double LongPixel[3][NPIXELS];
      double QtotPixel[3][NPIXELS];
      double PhiPixel[3][NPIXELS];
      double ElevPixel[3][NPIXELS];
      double Theta_pixel[NPIXELS];
      double LatPixelVert[3*NPIXELS];
      double LongPixelVert[3*NPIXELS];


      bool pixHasData[NMAXPAGES][NPIXELS];
      
      bool TelHit[NMAXPAGES][6];
      bool TelHasElveTrigger[6];
      
      bool IsThisFirstEvent;
      int HitTel[NMAXTELS];
      double Area_Ionosphere[NMAXTELS*NPIXELS]={};
      
      int NpixCount1;
      int NpixCount2;
      
      double BoltLatitude;
      double BoltLongitude;
      double BoltAltitude;
      double BoltGPStime;
      double TLEAltitude;
      

      TGraph *QC_vs_DistBoltRP[NROWS];
      TGraph *QC_vs_DistBoltCP[3*NCOLS];
      TGraph *QC_vs_ElliDistRP[NROWS];
      TGraph *QC_vs_ElliDistCP[3*NCOLS];
      TGraph *QC_vs_SinThBE_RP[NROWS];
      TGraph *QC_vs_SinThBE_CP[3*NCOLS];

      
      TH1F *hRawPixels[3][NPIXELS];
      TH1F *hCalPixels[3][NPIXELS];
      TH1F *hGCorrPixels[3][NPIXELS];
      TH1F *hTCorrPixels[3][NPIXELS];
      TH1F *hCalResiduals[3][NPIXELS];

      TH1D *hPixelPeakTimes[3][NROWS];
      
      TH1F *hRawTmax;
      TH1F *hRawQmax;
      TH1F *hRawTmax2;
      TH1F *hRawQmax2;
      TH1F *hFineTmax;
      
      TH1F *hFineQmax;
      TH1F *hFineQCmax;
      TH1F *hFineQCCmax;

      TH1F *hFineQtot;
      TH1F *hFineQCtot;
      TH1F *hFineQCCtot;

      TH1F *hFineSigmaT;
      TH1F *hFineSkew;
      
      TH1F *hFineTmax2;
      TH1F *hFineQmax2;

      TH1F *hFIT3Tmax;
      TH1F *hFIT3SigmaT;
      TH1F *hFIT3Qmax;
      TH1F *hFIT3Skew;
      TH1F *hFitPullsP;
      TH1F *hFitPullsN;
      
      TCanvas *cFIT3Plot;
      TCanvas *cFDViewPlot;
      TCanvas *cFDViewPlot2;


      TH2F *PeakAmpli = new TH2F("PeakAmpli"," Peak amplitude : row vs col",NCOLS,0,NCOLS,NROWS,0,NROWS);
      TH2F *PeakTime = new TH2F("PeakTime"," Peak Time (us) : row vs col",NCOLS,0,NCOLS,NROWS,0,NROWS);
      TH2F *PulseSigma = new TH2F("PulseSigma"," Pulse Sigma(us) Time : row vs col",NCOLS,0,NCOLS,NROWS,0,NROWS);
      TH2F *PulseSkew= new TH2F("PulseSkew"," Pulse Skewness(Right/Left) : row vs col",NCOLS,0,NCOLS,NROWS,0,NROWS);
      TH2F *PulseUS= new TH2F("PulseUS"," Pulse Undershoot : row vs col",NCOLS,0,NCOLS,NROWS,0,NROWS);
      
      TH1F *htrace = new TH1F("htrace"," Pixel Trace ",NMAXPAGES*50,0.,NMAXPAGES*100.);
      TH1F *hpull = new TH1F("hpull"," Fit Pull ",200,-10.,10.);
      
      
      //stringstream ss3f;
      //ss3f << " Surface Density of light emission ; Long(degrees); Lat(degrees)";
      
      
      //    TH2F *PixelsLL  = new TH2F("PixelsLL",ss3f.str().c_str(),100,-72.122,-59.878,100,-40,-30);
      //TH2F *PixelsLL  = new TH2F("PixelsLL"," Surface Density of light emission ; Long(degrees); Lat(degrees)",100,-78.244,-54.756,100,-45,-25);
      TH2F *PixelsLL  = new TH2F("PixelsLL"," Surface Density of light emission ; Long(degrees); Lat(degrees)",100,-70,-60,100,-38.2,-30);
      
      
      TH1F *hpixels[3*NPIXELS] ;
      TH1F *hresiduals[3*NPIXELS] ;
      TH1F *hpulls[3*NPIXELS] ;
      
      TH1F *hQC_vs_Dist[3];
      TH1F *hQC_vs_SinThBE[3];

    
      
      TH1F *hpixelsUC[3*NPIXELS] ;
      TH1F *TimeMax_vs_Col[NROWS];
      TH1F *NphMax_vs_Col[NROWS];
      TH1F *SigmaT_vs_Col[NROWS];
      TH1F *Skew_vs_Col[NROWS];
      
      TH1F *hChi2BoltFinder = new TH1F("hChi2BoltFinder"," #chi^{2} of the pre-fit for each row ",NROWS,0.5,0.5+NROWS);
      TH1F *hPhiMinBoltFinder = new TH1F("hPhiMinBoltFinder"," Phi of the first pixel in each row ",NROWS,0.5,0.5+NROWS);
      TH1F *hTMinBoltFinder = new TH1F("hTMinBoltFinder"," Time of the first pixel in each row ",NROWS,0.5,0.5+NROWS);
      
      TH2F *hCountElves = new TH2F("hCountElves"," One or two pulses per pixels? ; Npixels with 1 pulse; Npixels with 2 pulses ",NPIXELS,0,NPIXELS,NPIXELS,0,NPIXELS);

      
      TProfile *pQCC_vs_DistBolt =
             new TProfile("pQCC_vs_DistBolt","Corrected Nph/Area vs distance from bolt; ArcR(km); dNph/dArea",NDBbins,0., 500.,0.,3.e15);
      
      TH1F *hQCC_Asym = new TH1F("hQCC_Asym","Elve Azimuthal asymmetry; relative spread ",40,0.,2.);
      
      TH2F *hQCCAveSpread = new TH2F("hQCCAveSpread","Elve angular distribution ; relative spread ",NDP1,25.,525.,NDP2,0.,360.);

      TProfile *pQCC_vs_ElliDist[NDP1*NDP2];
      
      TH1F *hQCC_Spread[NDP1*NDP2];
      
      TF1 *fitFcn;
      
      TGraph *TC_vs_DistBoltRP[NROWS];
      TGraph *DTC_vs_DistBoltRP[NROWS];
      TGraph *DTC_vs_DistBoltCP[2*NCOLS];
      TGraph *DTC_vs_DistB2E_RP[NROWS];
      TGraph *DTC_vs_DistB2E_CP[2*NCOLS];

      TGraphErrors *BestTC_vs_DistBolt;

      
      long unsigned int TelGPSns[6]={};
      long unsigned int ThisGPSsec = 0;
      
      double Qmx[NMAXTELS*NPIXELS]={};
      double Tmx[NMAXTELS*NPIXELS]={};
      double R_EFD[NMAXTELS*NPIXELS]={};
      double R_BE[NMAXTELS*NPIXELS]={};
      double Dist_BE[NMAXTELS*NPIXELS]={};
      double SinTh_BE[NMAXTELS*NPIXELS]={};
      double dT_EFD[NMAXTELS*NPIXELS]={};
      double dT_BFD[NMAXTELS*NPIXELS]={};
      
      double QtotTrace = 0.;
      unsigned int EventNanoTime = 0;
      bool change_phi = false;
      
      double minphi = 360.;
      double maxphi = 0.;
      
      double min_charge = 999999999.;
      double max_charge = 0.;
      
      
      static const Int_t NRGBs = 5;
      static const Int_t MaxColors = 127;//255;
      Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
      Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
      Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
      Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };


      
    REGISTER_MODULE("FdProcessElves",FdProcessElves);
  };





#endif // _FdProcessElves_

// Configure (x)emacs for this file ...
// Local Variables:
// mode:c++
// compile-command: "make -C .. FdSDPFinderOG/FdSDPFinder.o  -k"
// End:
