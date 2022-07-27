// C++ includes
#include <iostream>
#include <vector>

// ROOT include
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TH2Poly.h"
#include "THStack.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TApplication.h"

// MaCh3 includes
#include "samplePDF/samplePDFND.h"


// *******************
// Class to hold the variations, throws and so on
class SampleSummary {
// *******************
  public:
    SampleSummary(int nSamples, const std::string &Outputfile, int nChainSteps);
    ~SampleSummary();

    void AddData(std::vector<TH2Poly*> &DataHist);
    void AddNominal(std::vector<TH2Poly*> &NominalHist, std::vector<TH2Poly*> &W2Nom);
    void AddThrow(std::vector<TH2Poly*> &MCHist, std::vector<TH2Poly*> &W2Hist, double LLHPenalty = 0.0, int DrawNumber = 0);
    void AddThrowByMode(std::vector<std::vector<TH2Poly*>> &SampleVector_ByMode);
    
    void Write();
    
    void SetLikelihood(int TestStatystic){ likelihood = TestStatystic;};
    void SetSamplePDF(samplePDFND* const sample){ sampleND = sample;};
    
  private:
    TRandom3* rnd;
    bool first_pass;

    //KS: We have two methods for Poissonian fluctuation
    bool StandardFluctuation;
    // Finalise the distributions from the thrown samples
    void inline MakePredictive();

    void inline PrepareTree();

    // Helper functions to calculate likelihoods for TH1D and TH2Ds
    void inline CalcLLH(TH2Poly * const & Data, TH2Poly * const & MC, TH2Poly * const & W2);
    void inline CalcLLH(TH1D * const & Data, TH1D * const & MC, TH1D * const & W2);

    double inline GetLLH(TH2Poly * const & Data, TH2Poly * const & MC, TH2Poly * const & W2);
    double inline GetLLH(TH1D * const & Data, TH1D * const & MC, TH1D * const & W2);

    // Helper functions to change titles etc of finished plots, calculate pvalues etc
    void inline MakeCutLLH();
    void inline MakeCutLLH1D(TH1D *Histogram, double llh_ref = -999);
    void inline MakeCutLLH2D(TH2D *Histogram);
    void inline MakeCutEventRate(TH1D *Histogram, double DataRate);
    void inline MakeChi2Hists();

    // Check the length of psyche samples agrees
    bool inline CheckPsycheSamples(int Length);

    // Helper to make ratio histograms
    template<class HistType> HistType* RatioHists(HistType* NumHist, HistType* DenomHist);
    // Helper to make ratio of TH2Polys
    TH2Poly* RatioPolys(TH2Poly* NumPoly, TH2Poly* DenomPoly);

    // Helper to project TH2D onto axis
    TH1D* ProjectHist(TH2D* Histogram, bool ProjectX);
    // Helper to project TH2Poly onto axis
    TH1D* ProjectPoly(TH2Poly* Histogram, bool ProjectX, int selection, bool MakeErrorHist = false);
  
    // Make Poisson flcutuation of TH1D hist
    TH1D* MakeFluctuatedHistogram(TH1D* PolyHist);
    TH1D* MakeFluctuatedHistogramStandard(TH1D* PolyHist);
    TH1D* MakeFluctuatedHistogramAlternative(TH1D* PolyHist);
    
    // Make Poisson flcutuation of TH2Poly hist
    TH2Poly* MakeFluctuatedHistogram(TH2Poly* PolyHist);
    TH2Poly* MakeFluctuatedHistogramStandard(TH2Poly* PolyHist);
    TH2Poly* MakeFluctuatedHistogramAlternative(TH2Poly* PolyHist); 
    
    //Return 2 random numbers along axis x and y distributed according to the cell-contents
    int GetRandomPoly2(TH2Poly* PolyHist);
    
    // Get the mode error from a TH1D
    double GetModeError(TH1D * const hpost);
    
    // Helper to Normalise histograms
    template<class HistType> HistType* NormaliseHist(HistType* Histogram);

    // Vector of vectors which holds the loaded MC histograms
    std::vector<std::vector<TH2Poly*> > MCVector;
    std::vector<std::vector<TH2Poly*> > W2MCVector;
    // Vector to hold the penalty term
    std::vector<double> LLHPenaltyVector;

    // Number of samples
    __int__ nSamples;
    // Number of Modes
    __int__ nModes;

    // A map to keep track of what psyche sample indices we want. Save some read time
    std::vector<int> PsycheSampleMap;
    std::vector<std::string> SampleNames;

    // The posterior predictive for the whole selection: this gets built after adding in the toys. Now an array of Th1ds, 1 for each poly bin, for each sample, and the same for W2
    TH1D ***PosteriorHist;
    TH1D ***w2Hist;
        
    // The data histogram for the selection
    TH2Poly **DataHist;
    TH1D **DataHist_ProjectX;
    TH1D **DataHist_ProjectY;
    // The nominal histogram for the selection
    TH2Poly **NominalHist;
    // The w2 histograms
    TH2Poly **W2NomHist;
    TH2Poly **W2MeanHist;
    TH2Poly **W2ModeHist;

    // The histogram containing the lnL for each throw
    TH1D *lnLHist;
    // The lnLhist for the draw vs MC fluctuated
    TH1D *lnLHist_drawfluc;
    // The lnLhist for the draw vs draw fluctuated
    TH1D *lnLHist_drawflucdraw;
    // The lnLhist for the draw vs data
    TH1D *lnLHist_drawdata;
    // The 2D lnLhist, showing (draw vs data) and (draw vs fluct), anything above y=x axis is the p-value
    TH2D *lnLDrawHist;
    // The 2D lnLHist, showing (draw vs data) and (draw vs draw fluct), anything above y=x axis is the p-value
    TH2D *lnLFlucHist;
    
    // The 2D lnLHist but for ProjectionX histogmram (pmu), showing (draw vs data) and (draw vs draw fluct), anything above y=x axis is the p-value
    TH2D *lnLFlucHist_ProjectX;

    // The histogram containing the lnL (draw vs data) for each throw for each sample
    TH1D **lnLHist_Sample_DrawData;
    // The histogram containing the lnL (draw vs draw fluct) for each throw for each sample
    TH1D **lnLHist_Sample_DrawflucDraw;
    // The histogram containing the lnL (draw vs pred fluct) for each throw for each sample
    TH1D **lnLHist_Sample_PredflucDraw;

    // The LLH distribution in pmu cosmu for using the mean in each bin
    TH2Poly **lnLHist_Mean;
    // The LLH distribution in pmu cosmu for using the mode in each bin
    TH2Poly **lnLHist_Mode;
    
    // The LLH distribution in pmu using the mean in each bin
    TH1D **lnLHist_Mean_ProjectX;
    
    // The posterior predictive distribution in pmu cosmu using the mean
    TH2Poly **MeanHist;
    // The posterior predictive distribution in pmu cosmu using the mode
    TH2Poly **ModeHist;

    // Holds the event rate for the distribution
    TH1D **SumHist;
    // Holds the total event rate
    TH1D *EventHist;
    // Holds the bin-by-bin LLH for the mean posterior predictive vs the data
    TH1D **lnLHist_Mean1D;
    // Holds the bin-by-bin LLH for the mode posterior predictive vs the data
    TH1D **lnLHist_Mode1D;

    // Holds the history of which entries have been drawn in the MCMC file
    TH1D *RandomHist;

    // Number of throws by user
    int nChainSteps;

    //bool whether we have Prior or Posterior Predictive 
    bool isPriorPredictive;
    
    // Number of throws
    int nThrows;
    
    // Max Number of Bins per each sample
    int* maxBins;

    // Total LLH for the posterior predictive distribution
    double llh_total;

    // Output filename
    std::string OutputName;
    TFile *Outputfile;

    // TTree which we save useful information to
    TTree *OutputTree;
    // Data vs Draw
    double *llh_data_draw;
    // Data vs Fluctuated Draw
    double *llh_data_drawfluc;
    // Data vs Fluctuated Predictive
    double *llh_data_predfluc;

    // Draw vs Predictive
    double *llh_draw_pred;
    // Fluctuated Predicitve vs Draw
    double *llh_predfluc_draw;
    // Fluctuated Draw vs Predictive
    double *llh_drawfluc_pred;
    // Fluctuated Draw vs Draw
    double *llh_drawfluc_draw;
    // Fluctuated Predictive vs Predictive
    double *llh_predfluc_pred;
    // Fluctuated Draw vs Fluctuated Predictive
    double *llh_drawfluc_predfluc;
    // Fluctuated Data vs Draw
    double *llh_datafluc_draw;

    // Projection X (most likely muon momentum) of LLH
    double *llh_data_draw_ProjectX;
    double *llh_drawfluc_draw_ProjectX;
    
    // LLH penalty for each throw
    double llh_penalty;
    // Event rate for each throw
    double event_rate;

    // Data vs Draw
    double total_llh_data_draw;
    // Fluctuated Draw vs Draw
    double total_llh_drawfluc_draw;
    // Fluctuated Predicitve vs Draw
    double total_llh_predfluc_draw;
    
    // Data vs Fluctuated Predictive
    double total_llh_data_predfluc;
    // Data vs Fluctuated Draw
    double total_llh_data_drawfluc;
    // Draw vs Predictive
    double total_llh_draw_pred;
    // Fluctuated Draw vs Predictive
    double total_llh_drawfluc_pred;
    // Fluctuated Draw vs Fluctuated Predictive
    double total_llh_drawfluc_predfluc;
    // Fluctuated Data vs Draw
    double total_llh_datafluc_draw;
    // Fluctuated Predictive vs Predictive
    double total_llh_predfluc_pred;
    
    // Data vs Draw for projection X (most likely muon momentum)
    double total_llh_data_draw_ProjectX;
    // Fluctuated Draw vs Draw for projection X (most likely muon momentum)
    double total_llh_drawfluc_draw_ProjectX;
    
    //By mode variables
    bool DoByModePlots; 
    // The posterior predictive distribution in pmu cosmu using the mean
    TH2Poly ***MeanHist_ByMode;
    TH1D ****PosteriorHist_ByMode;

    samplePDFND* sampleND;
    
    //Type of likelihood for example Poisson, Barlow-Beestion or Ice Cube
    int likelihood;
};
