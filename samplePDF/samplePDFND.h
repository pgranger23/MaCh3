#ifndef _samplePDFND_h_
#define _samplePDFND_h_

// Do we use TF1 or TSpline3* for spline evaluations
#define USE_TSpline3_red 1
#define USE_Akima_Spline 2
#define USE_Truncated_TSpline3 3
#define USE_Truncated_Akima 4
#define USE_Monotone_Spline 5
#define USE_TSpline3 6
#define USE_TF1 7
#define USE_TF1_red 8

// Can use:
//  TSpline3 (third order spline in ROOT)
//  TSpline3_red (reduced class third order spline in ROOT)
//  TF1 (fifth order poly in ROOT)
//  TF1_red (reduced class fifth order poly)
//  (Experimental) Akima_Spline (crd order spline which is allowed to be discontinuous in 2nd deriv) 
//  (Experimental) Truncated_Spline (third order root spline but if eval value outside of outer knots, coefficients are set to 0)
//  (Experimental) Truncated_Akima_Spline (same as Truncated_Spline but uses akima spline interpolation for internal segments)

#define USE_SPLINE USE_TSpline3_red

// Set the __SPLINE_TYPE__ accordingly
#if USE_SPLINE == USE_TSpline3
#define __SPLINE_TYPE__ TSpline3
#elif USE_SPLINE == USE_TSpline3_red
#define __SPLINE_TYPE__ TSpline3_red
#elif USE_SPLINE == USE_TF1
#define __SPLINE_TYPE__ TF1
#elif USE_SPLINE == USE_TF1_red
#define __SPLINE_TYPE__ TF1_red

#elif USE_SPLINE == USE_Akima_Spline
#define __SPLINE_TYPE__ Akima_Spline
#elif USE_SPLINE == USE_Truncated_TSpline3
#define __SPLINE_TYPE__ Truncated_Spline
#elif USE_SPLINE == USE_Truncated_Akima
#define __SPLINE_TYPE__ Truncated_Akima_Spline
#elif USE_SPLINE == USE_Monotone_Spline
#define __SPLINE_TYPE__ Monotone_Spline
#endif

// Do we want to debug the TF1 treatment? This involves writing TSpline3* and TF1* to file and comparing how well they fit TGraph*
//#define DEBUG_TF1

// C++ includes
#include <string>
#include <cmath>

// ROOT include
#include "TChain.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TMatrixDSym.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"

// MaCh3 samplePDF includes
#include "samplePDFBase.h"
#include "Structs.h"
#include "NDMCStruct.h"

// MaCh3 covariance includes
#include "covariance/covarianceXsec.h"

// Include the manager for an alternate constructor
#include "manager/manager.h"

#ifdef CUDA
#include "splines/SplineMonolith.h"
#endif

class samplePDFND : public samplePDFBase {
  public:
      samplePDFND(manager *FitManager);
      ~samplePDFND();

      //KS: Each experiment uses some specyfic variables, make template fucnion for it
      void InitExperimentSpecific();

      // The different covariance matrices to be associated with the samplePDF
      void setXsecCov(covarianceXsec * const xsec_cov);
      //WARNING FIXME TODO T2K specyfic
      //void setSimpleDetCov(covarianceNDDetPoly * const indet) { NDDetCov = indet; }
      //covarianceNDDetPoly * const GetSimpleDetCov() const { return NDDetCov; }

      // Function to set the Asimov fake-data from within samplePDF
      // Put this back in
      void setAsimovFakeData(bool CustomReWeight = false);
      void setDataFromFile(std::string &FileName);
      void setAsimovFakeData_FromFile(std::string &FileName);
      void setAsimovFakeDataThrow();
      void setAsimovFakeDataFluctuated(bool CustomReWeight = false);

      double getLikelihood();
      double getTestStatLLH(double data, double mc, double w2);
      double getSampleLikelihood(int isample);
      double getEventRate(int index); 

      // The reweighting functions
      void reweight(double *oscpar);
      void reweight(double *oscpar, double *oscpar2) { reweight(oscpar); };
      
      void addData(std::vector<double> &dat) {return;}
      void addData(std::vector< vector <double> > &dat) {return;}
      void addData(TH1D* binneddata){return;}
      void addData(TH2Poly* binneddata){return;}
      
      void addData(std::vector<double> &dat, int index);
      void addData(std::vector< std::vector <double> > &dat, int index);
      void addData(TH1D* binneddata, int index);
      void addData(TH2Poly* binneddata, int index);

      // Randomize the starting position
      void RandomStart();
      
      // Getters for the event histograms
      TH1* getPDF(int Selection);
      TH1* getData(int Selection) { return (TH1*)datapdfs->At(Selection); };
      TH2Poly* getW2(int Selection);
      TH1* getPDFMode(int Selection, int Mode) { return (TH1*)((TObjArray*)samplemodepdfs->At(Selection))->At(Mode); };

      void GetKinVars(int sample, KinematicTypes &TypeX, KinematicTypes &TypeY);

      // Setup the binning for a given sample
      void SetupBinning(int Selection, std::vector<double> &BinningX, std::vector<double> &BinningY);

      void fillDataFromSamples();
      void fillReweightingBins();

      virtual void EnableModeHistograms();

      void printRates(bool dataonly = false);

#ifdef CUDA
      void fillGPUSplines();
      // The monolith
      SMonolith *splineMonolith;
#endif

  protected:
      //KS: Find Detector bin, hist bin and other useful information using multithreading
      inline void FindAdditionalInfo();
      
      // Helper function to find normalisation bins for normalisation parameters that aren't simply a mode scaling
      // e.g. different normalisation parameters for carbon and oxygen, neutrino and anti-neutrino, etc
      inline void FindNormBins();

      //KS: Find pointer for each norm dial to reduce impact of covarianceXsec::calcReweight()
      inline void FindNormPointer();

      // Reserve spline memory 
      void ReserveMemory(int nEve);
      // Prepare weights
      void PrepareWeights();

      inline void LoadSamples();

      // virtual for GPU
      void SetSplines(TGraph** &xsecgraph, const int i);
#if USE_SPLINE < USE_TF1
      void SetSplines_Reduced(TGraph** &xsecgraph, const int i);
#endif
      // Helper function to check if the covariances have been set
      inline void CheckCovariances();

      // Perform the main reweight loop
#ifdef MULTITHREAD
      void ReWeight_MC_MP();
#else
      void ReWeight_MC();
#endif

      virtual void ReconfigureFuncPars(){};
      virtual void CalcFuncPars(int Event){};

      // Helper function to reset histograms
      inline void ResetHistograms();

      //KS: Helper which initlaise PDF
      inline void InitialisePDF();

      //Helper which udpate data arrays from histogram
      inline void UpdateDataPDF();

      // Calculate the cross-section weight for a given event
      double CalcXsecWeight(const int EventNumber);
      // Calculate the spline weight for a given event
      // virtual for GPU
      virtual double CalcXsecWeight_Spline(const int EventNumber);
      // Calculate the norm weight for a given event
      double CalcXsecWeight_Norm(const int EventNumber);

      bool HaveIRandomStart; // Have I random started?

      // Pointer to fit manager
      manager *FitManager;

      // The covariance classes
      covarianceBase* NDDetCov;

      // This is the number of cross-section splines we're loading
      // This one is read from setCovMatrix function which looks at the input covariance and counts the number of cross-section splines
      __int__ nXsecSplines;

      // Number of MC events are there
      unsigned int nEvents;

      // what is the maximum number of bins we have
      __int__* maxBins;



      // Struct containing the ND280 information
      std::vector<ND280EVENT> NDEve;
      std::vector<ND280EVENT_AUXILIARY> NDEve_Aux;

      // Dimensions of the ith psyche selection (2D or 1D)
      int* ndims;
      // The kinematic type we're plotting
      KinematicTypes **kinvars;
      
      // The pdfs we're storing
      // MC pdfs
      TObjArray* samplepdfs;
      // data pdfs
      TObjArray* datapdfs;
      // A vector of TH2Ds to hold the weight^2 for e.g. Barlow Beeston
      std::vector<TH2Poly*> W2Hist;
      // mode of MC pdfs
      TObjArray* samplemodepdfs;
      //PolyBin for each sample and each Pmu/cosTheta Mu bin, we need this for Eb
      TH2PolyBin ***polybins;
    
      //KS:: We use below only for transfering events from private array but one can imagine reaplaing TH2Poly class just with those double**, as all inofrmation we have stored in MaxBin etc, and other stuff.
      // The per-thread array
      double **samplePDF_data_array;
      double **samplePDF_array;
      // array of weight^2
      double **samplePDF_w2_array;

#ifdef MULTITHREAD
      // And the mode
      double ***samplePDF_mode_array;
#endif

      TObjArray** modeobjarray;
      // do we want mode MC pdf to be save
      bool modepdf;

      // Struct containing the cross-section info
      XSecStruct<__SPLINE_TYPE__*>* xsecInfo;

      // number of cross-section normalisation params?
      int nxsec_norm_modes;

      // Number of spline parameters
      int nSplineParams;
      std::vector<int> splineParsIndex;
      std::vector<std::string> splineParsNames;
      std::vector<std::string> splineFileParsNames;

      // Number of spline parameters that aren't repeated
      int nSplineParamsUniq;
      std::vector<int> splineParsUniqIndex;
      std::vector<std::string> splineParsUniqNames;

      // Bit-field comparison
      std::vector<int> xsecBitField;
      std::vector<int> linearsplines;

#ifdef CUDA
      // Clean up the memory temporarily allocated for GPU preparation
      inline void CleanUpMemory();
      #ifdef DEBUG_DUMP
      void CompareCPU_GPU_Splines(const int EventNumber);
      // Holds the CPU and GPU weights
      TFile *DebugFile;
      TH1D** gpu_weights_plot;
      TH1D** cpu_weights_plot;
      TH1D** diff_weights_plot;
      int badWeight;
      int nReconf;
      #endif

      #if USE_SPLINE < USE_TF1
      int *segments;
      #endif
      float *vals;
#endif

      //WARNING T2K Specyfic
      /*
      //KS: Use 2D or 1D Binding Energy
      bool Use2dEb;

      // String with what production we want
      TString *prod;

      bool UseSandMC; //whether to use sand or not
      //set if you want to use sand or not
      void setSandMC(bool useSand){ UseSandMC=useSand;};
      bool UseSandMC; //whether to use sand or not
      // Use covariance matrix for detector
      bool simple;
      */


      // Array of FastSplineInfo structs: keeps information on each xsec spline for fast evaluation
      // Method identical to TSpline3::Eval(double) but faster because less operations
#if USE_SPLINE < USE_TF1
      FastSplineInfo *SplineInfoArray;
      template <class T> 
        double FastSplineEval(T* spline, const int SplineNumber);
      void FindSplineSegment();
#endif

#ifdef DEBUG_TF1
      TFile *dumpfile;
#endif

#if DEBUG > 0
      // Dump spline information
      inline void DumpSplines(double xsecw, double detw, int EventNumber);
      // Print information about the xsec and ndobj structs
      inline void PrintStructs(double xsecw, double detw, int i);

      unsigned int nEventsInTree;
      // Vector of POT weights needed for pot+xsec weights
      std::vector<double> POTWeights;

      // 2D histograms which match BANFF's numbers
      std::vector<TH2Poly*> MCOnly;
      std::vector<TH2Poly*> POTOnly;
      std::vector<TH2Poly*> FluxOnly;
      std::vector<TH2Poly*> XsecOnly;
      std::vector<TH2Poly*> DetOnly;
      std::vector<TH2Poly*> NDCovOnly;
      std::vector<TH2Poly*> AllOnly;

      // TH1D which matches NuMu group's binning
      std::vector<TH1D*> NuMuPmu;

      // Events per normalisation parameter
      std::vector<unsigned int> EventsPerNorm;
      // Events per normalisation mode
      std::vector<unsigned int> EventsPerMode;
#endif

  private:
  int blah;
};

#if USE_SPLINE < USE_TF1
// ***************************************************************************
// Fast spline evaluation
// Inspired by TSpline3::Eval, which we can speed up considerably
// Main reason is that we know that for one parameter (e.g. MAQE) we will have the same number of points, x min, xmax, etc for all MAQE splines, so we can signficantly reduce number of operations
// The curious can find very similar GPU code in splines/gpuSplineUtils.cu and CPU code in spline/SplineMonolith.cpp::Eval
// I've included it here for more transparency: this kind of eval should be possible for SK splines too
template <class T>
double samplePDFND::FastSplineEval(T* spline, const int SplineNumber) {
// ***************************************************************************

  // Check if the spline is NULL
  if (spline == NULL) return 1.0;

  // The segment has already been found in FindSplineSegment()
  int segment = SplineInfoArray[SplineNumber].CurrSegment;

  //KS: Each spline segment is described by 3rd order polynomial which has the following form:
  //f(xvar) = weight = d*dx(xvar)^3 + c*dx(xvar)^2 + b*dx(xvar) + y, 
  //where dx(xvar) = xvar - x 
  //This is the meaning of coefficients that we can extract from spline segment
  double x = -999.99;
  double y = -999.99;
  double b = -999.99;
  double c = -999.99;
  double d = -999.99;

  // Now write the coefficients
  spline->GetCoeff(segment, x, y, b, c, d);

  // Get the variation for this reconfigure for the ith parameter

  double xvar = XsecCov->calcReWeight(splineParsIndex[SplineNumber]);

  // The Delta(x)
  double dx = xvar - x;

  // The spline weight to return
  double weight = y+dx*(b+dx*(c+d*dx));

  // Check that eval on the TSpline3 is the same as our weight
#ifdef DEBUG
  int GlobalIndex = splineParsIndex[SplineNumber];
  // Difference between eval and weight
  double diff = fabs(spline->Eval(XsecCov->calcReWeight(GlobalIndex)) - weight);

  if (diff > 1.E-7) {

    std::cerr << "TSpline3->Eval() != custom eval: Something is wrong with FastSplineEval!" << std::endl;
    std::cerr << "Difference in spline evaluation > 1.E-5!" << std::endl;

    std::cerr << "Eval      = " << spline->Eval(XsecCov->calcReWeight(GlobalIndex)) << std::endl;
    std::cerr << "Cust      = " << weight << std::endl;
    std::cerr << "diff      = " << diff << std::endl;
    std::cerr << "param     = " << splineParsNames[SplineNumber] << std::endl;
    std::cerr << "variation = " << XsecCov->calcReWeight(GlobalIndex) << std::endl;
    std::cerr << "paramVal  = " << XsecCov->getParProp(GlobalIndex) << std::endl;

    // Check we've found the right segment
    int klow = spline->FindX(xvar);
    if (klow >= spline->GetNp()-1 && spline->GetNp() > 1) klow = spline->GetNp()-2;

    std::cerr << "segment   = " << segment << std::endl;
    std::cerr << "spl segm  = " << klow << std::endl;
    std::cerr << "nPoints   = " << SplineInfoArray[SplineNumber].nPts << std::endl;
    std::cerr << "nPoints sp= " << spline->GetNp() << std::endl;


    std::cerr << "Printing x information:" << std::endl;
    for (int i = 0; i < SplineInfoArray[SplineNumber].nPts; ++i) {
      std::cerr << "   " << i << " = " << SplineInfoArray[SplineNumber].xPts[i] << std::endl;
    }
    std::cerr << "From spline: " << std::endl;
    for (int i = 0; i < spline->GetNp(); ++i) {
      double xtmp, ytmp;
      spline->GetKnot(i, xtmp, ytmp);
      std::cerr << "   " << i << " = " << xtmp << std::endl;
    }

    if (klow != segment) {
      std::cerr << "Have found wrong segment in FindSplineSegment!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
#endif
  return weight;
}
#endif
#endif
