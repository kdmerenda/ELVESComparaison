#ifndef _ComparaisonStudyNS_ComparaisonStudy_h_
#define _ComparaisonStudyNS_ComparaisonStudy_h_

/**
 * \file
 * \author Kevin-Druis Merenda
 * \date 29 Apr 2017
 */

#include <fwk/VModule.h>
#include <boost/utility.hpp>
#include "TH1D.h"
#include "TFile.h"
/*
 * Avoid using using namespace declarations in your headers,
 * doing so makes all symbols from each namespace visible
 * to the client which includes your header.
 */

namespace ComparaisonStudyNS {

  class ComparaisonStudy : 
      public boost::noncopyable,
      public fwk::VModule {
  public:
    ComparaisonStudy();
    ~ComparaisonStudy();
    VModule::ResultFlag Init();
    VModule::ResultFlag Run(evt::Event& e);
    VModule::ResultFlag Finish();
  private:
    // Declare down here your data following this conventions:
    // fFieldName  : members.
    // fgFieldName : static data members.
    const static int fNPixels = 440;
    const static int fNTels = 6;
    TH1D* fhRawPixelData[fNTels][fNPixels];
    std::vector<std::string> eventChecker;
    int fEventCounter=0;
    TFile* outputPlots;
    bool fData;
    // Declare down here your private functions like this:
    //
    // /// Brief description.
    // /*! Detailed description that can span several 
    //     lines.
    // */
    // type foo();

    // This goes at the end.
    REGISTER_MODULE("ComparaisonStudy",ComparaisonStudy);
  };
}

#endif // _ComparaisonStudyNS_ComparaisonStudy_h_
