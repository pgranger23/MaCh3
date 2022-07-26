#include "SampleSummary.h"

// *******************
// The constructor
SampleSummary::SampleSummary(const int PsycheSamples, const std::string &Filename, int nSteps) {
// *******************

  std::cout << "Making sample summary class..." << std::endl;
  #ifdef MULTITHREAD
  std::cout << "with OpenMP and " << omp_get_max_threads() << " threads" << std::endl;
  #endif
  
  StandardFluctuation = true;
  
  if(StandardFluctuation) std::cout << "Using standard method of statisticall fluctuation" << std::endl;
  else std::cout << "Using alternative method of statisticall fluctuation, which is much slower" << std::endl;
  
  nChainSteps = nSteps;
  //KS: nChainSteps == 0 means we run PriorPredcitive
  
  if(nChainSteps == 0) isPriorPredictive = true;
  else isPriorPredictive = false;
  
  OutputName = Filename;
  nSamples = PsycheSamples;
  nThrows = 0;
  first_pass = true;
  Outputfile = NULL;
  OutputTree = NULL;

  rnd = new TRandom3();

  nModes = GetMaCh3Modes();

  DataHist = new TH2Poly*[nSamples];
  DataHist_ProjectX = new TH1D*[nSamples];
  DataHist_ProjectY = new TH1D*[nSamples];
  NominalHist = new TH2Poly*[nSamples];
  PosteriorHist = new TH1D**[nSamples];
  W2NomHist = new TH2Poly*[nSamples];
  w2Hist = new TH1D**[nSamples];

  maxBins = new int[nSamples];
  
  lnLHist_Mean = new TH2Poly*[nSamples];
  lnLHist_Mode = new TH2Poly*[nSamples];
  lnLHist_Mean_ProjectX = new TH1D*[nSamples];
  MeanHist = new TH2Poly*[nSamples];
  ModeHist = new TH2Poly*[nSamples];
  W2MeanHist = new TH2Poly*[nSamples];
  W2ModeHist = new TH2Poly*[nSamples];
  SumHist = new TH1D*[nSamples];
  lnLHist_Mean1D = new TH1D*[nSamples];
  lnLHist_Mode1D = new TH1D*[nSamples];
  lnLHist_Sample_DrawData = new TH1D*[nSamples];
  lnLHist_Sample_DrawflucDraw = new TH1D*[nSamples];
  lnLHist_Sample_PredflucDraw = new TH1D*[nSamples];
    
  //KS: When a histogram is created with an axis lower limit greater or equal to its upper limit ROOT will automitaclly adjust histogram range
  // https://root.cern.ch/doc/master/classTH1.html#auto-bin
  lnLHist = new TH1D("lnLHist_predpredfluc", "lnLHist_predpredfluc", 100, 1, -1);
  lnLHist->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Pred)");
  lnLHist->GetYaxis()->SetTitle("Counts");

  lnLHist_drawdata = new TH1D("lnLHist_drawdata", "lnLHist_drawdata", 100, 1, -1);
  lnLHist_drawdata->GetXaxis()->SetTitle("-2LLH (Data, Draw)");
  lnLHist_drawdata->GetYaxis()->SetTitle("Counts");

  lnLHist_drawfluc = new TH1D("lnLHist_drawpredfluc", "lnLHist_drawpredfluc", 100, 1, -1);
  lnLHist_drawfluc->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Draw)");
  lnLHist_drawfluc->GetYaxis()->SetTitle("Counts");

  lnLHist_drawflucdraw = new TH1D("lnLHist_drawflucdraw", "lnLHist_drawflucdraw", 100, 1, -1);
  lnLHist_drawflucdraw->GetXaxis()->SetTitle("-2LLH (Draw Fluc, Draw)");
  lnLHist_drawflucdraw->GetYaxis()->SetTitle("Counts");

  lnLDrawHist = new TH2D("lnLDrawHist", "lnLDrawHist", 50, 1, -1, 50, 1, -1);
  lnLDrawHist->GetXaxis()->SetTitle("-2LLH_{Pred Fluc, Draw}");
  lnLDrawHist->GetYaxis()->SetTitle("-2LLH_{Data, Draw}");

  lnLFlucHist = new TH2D("lnLFlucHist", "lnLFlucHist", 50, 1, -1, 50, 1, -1);
  lnLFlucHist->GetXaxis()->SetTitle("-2LLH_{Draw Fluc, Draw}");
  lnLFlucHist->GetYaxis()->SetTitle("-2LLH_{Data, Draw}");

  EventHist = new TH1D("EventHist", "Total Event Rate", 100, 1, -1);
  EventHist->GetXaxis()->SetTitle("Total event rate");
  EventHist->GetYaxis()->SetTitle("Counts");
  EventHist->SetLineWidth(2);

  //KS: at this point we haven't passed samplePDF ND yet, could be easyly changed but I am lazy right now
  KinematicTypes Type1 = kLeptonMomentum;
  lnLFlucHist_ProjectX = new TH2D("lnLFlucHist_ProjectX", "lnLFlucHist_ProjectX", 50, 1, -1, 50, 1, -1);
  lnLFlucHist_ProjectX->GetXaxis()->SetTitle(("-2LLH_{Draw Fluc, Draw} for " + Kinematic_ToString(Type1)).c_str());
  lnLFlucHist_ProjectX->GetYaxis()->SetTitle(("-2LLH_{Data, Draw} for " + Kinematic_ToString(Type1)).c_str());
  
  // Holds the hist of random number draws, only works for posterior predictive
  if(!isPriorPredictive)
  {
    RandomHist = new TH1D("RandomHist", "RandomHist", 100, 0, nChainSteps);
    RandomHist->GetXaxis()->SetTitle("Step");
    double binwidth = nChainSteps/RandomHist->GetNbinsX();
    std::stringstream ss;
    ss << "Draws/" << binwidth;
    RandomHist->GetYaxis()->SetTitle(ss.str().c_str());
    RandomHist->SetLineWidth(2);
  }
  for (int i = 0; i < nSamples; ++i) 
  {
    PosteriorHist[i] = NULL;
    w2Hist[i] = NULL;
      
    DataHist[i] = NULL;
    DataHist_ProjectX[i] = NULL;
    DataHist_ProjectY[i] = NULL;
    NominalHist[i] = NULL;

    SumHist[i] = NULL;
    MeanHist[i] = NULL;
    W2MeanHist[i] = NULL;
    W2ModeHist[i] = NULL;
    lnLHist_Mean[i] = NULL;
    lnLHist_Mode[i] = NULL;
    lnLHist_Mean_ProjectX[i] = NULL;
    lnLHist_Mean1D[i] = NULL;
    lnLHist_Mode1D[i] = NULL;
    lnLHist_Sample_DrawData[i] = NULL;
    lnLHist_Sample_DrawflucDraw[i] = NULL;
    lnLHist_Sample_PredflucDraw[i] = NULL;
  }//end loop over samples
  
  DoByModePlots = false;
  MeanHist_ByMode = NULL;
  PosteriorHist_ByMode = NULL;
}

// *******************
//  Destructor
SampleSummary::~SampleSummary() {
// *******************
    delete rnd;
    
    delete lnLHist;
    delete lnLHist_drawdata;
    delete lnLHist_drawfluc;
    delete lnLHist_drawflucdraw;
    delete lnLDrawHist;
    delete lnLFlucHist;
    delete EventHist;

    delete lnLFlucHist_ProjectX;
    for (int i = 0; i < nSamples; ++i) 
    {
        if(DataHist[i] == NULL) continue;
        delete DataHist[i];
        delete NominalHist[i];
        delete SumHist[i];
        delete MeanHist[i];
        delete W2MeanHist[i];
        delete W2ModeHist[i];
        
        delete lnLHist_Mean[i];
        delete lnLHist_Mode[i];
        delete lnLHist_Mean_ProjectX[i];
        delete lnLHist_Mean1D[i];
        delete lnLHist_Mode1D[i];
        delete lnLHist_Sample_DrawData[i];
        delete lnLHist_Sample_DrawflucDraw[i];
        delete lnLHist_Sample_PredflucDraw[i];
        
        for (int j = 0; j < maxBins[i]; ++j) 
        {   
            delete PosteriorHist[i][j];
            delete w2Hist[i][j];
        }
        delete[] PosteriorHist[i];
        delete[] w2Hist[i];
    }
    delete[] DataHist;
    delete[] NominalHist;
    delete[] SumHist;
    delete[] MeanHist;
    delete[] W2MeanHist;
    delete[] W2ModeHist;
    
    delete[] lnLHist_Mean;
    delete[] lnLHist_Mode;
    delete[] lnLHist_Mean_ProjectX;
    delete[] lnLHist_Mean1D;
    delete[] lnLHist_Mode1D;
    delete[] lnLHist_Sample_DrawData;
    delete[] lnLHist_Sample_DrawflucDraw;
    delete[] lnLHist_Sample_PredflucDraw;
    
    delete[] maxBins;
        
    delete[] PosteriorHist;
    delete[] w2Hist;
    
    if(!isPriorPredictive)  delete RandomHist;

    delete[] llh_data_draw;
    delete[] llh_data_drawfluc;
    delete[] llh_data_predfluc;
    delete[] llh_draw_pred;
    delete[] llh_drawfluc_pred;
    delete[] llh_drawfluc_predfluc;
    delete[] llh_drawfluc_draw;
    delete[] llh_predfluc_pred;
    delete[] llh_predfluc_draw;
    delete[] llh_datafluc_draw;
    
    delete[] llh_data_draw_ProjectX;
    delete[] llh_drawfluc_draw_ProjectX;
    
    if(DoByModePlots)
    {
        for (int i = 0; i < nSamples; ++i) 
        {
            for (int j = 0; j < nModes+1; j++)
            {
                for (int k = 0; k < maxBins[i]; ++k) 
                {   
                    delete PosteriorHist_ByMode[i][j][k];
                }
                delete[] PosteriorHist_ByMode[i][j];
                delete MeanHist_ByMode[i][j];
            }
            delete[] PosteriorHist_ByMode[i];
            delete[] MeanHist_ByMode[i];
        }
        delete[] PosteriorHist_ByMode;
        delete[] MeanHist_ByMode;
    }
}

// *******************
// Check size of psyche sample against size of vectors
bool SampleSummary::CheckPsycheSamples(int Length) {

// *******************
  bool ok = (sampleND->GetNsamples() == Length);
  if (!ok) {
    std::cerr << "Size of SampleVector input != number of psyche samples" << std::endl;
    std::cout << "Size of SampleVector: " << Length << std::endl;
    std::cout << "Size of SamplePDF samples: " << sampleND->GetNsamples() << std::endl;
    std::cerr << "Something has gone wrong with making the Samples" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }
  return ok;
}

// *******************
// Add a data histogram to the list (will have N_psyche_samples of these)
// Since the data doesn't change with varying the MC
void SampleSummary::AddData(std::vector<TH2Poly*> &Data) {
// *******************
  int Length = Data.size();
  // Check length of psyche samples are OK
  if (!CheckPsycheSamples(Length)) throw;
  for (int i = 0; i < Length; ++i) {
    if (Data[i] == NULL) {
      DataHist[i] = NULL;
      DataHist_ProjectX[i] = NULL;
      DataHist_ProjectY[i] = NULL;
      maxBins[i] = 0;
    } else {
      PsycheSampleMap.push_back(i);
      DataHist[i] = (TH2Poly*)(Data[i]->Clone());
      DataHist_ProjectX[i] = ProjectPoly(DataHist[i], true, i);
      DataHist_ProjectY[i] = ProjectPoly(DataHist[i], false, i);
      maxBins[i] = DataHist[i]->GetNumberOfBins();
    }
  }
}

// *******************
// Add the nominal histograms to the list (will have N_psyche_samples of these)
void SampleSummary::AddNominal(std::vector<TH2Poly*> &Nominal, std::vector<TH2Poly*> &NomW2) {
// *******************

  int Length = Nominal.size();
  if (!CheckPsycheSamples(Length)) throw;
  // Loop over the lenght of nominal and set the initial distributions up
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < Length; ++i) {
  // If NULL it indicates the selection was turned off, so initialise all the hists to NULL
    if (Nominal[i] == NULL) {
      PosteriorHist[i] = NULL;
      NominalHist[i] = NULL;
      W2NomHist[i] = NULL;
      lnLHist_Mean[i] = NULL;
      lnLHist_Mode[i] = NULL;
      lnLHist_Mean_ProjectX[i] = NULL;
      MeanHist[i] = NULL;
      ModeHist[i] = NULL;
      W2MeanHist[i] = NULL;
      W2ModeHist[i] = NULL;
      SumHist[i] = NULL;
      lnLHist_Sample_DrawData[i] = NULL;
      lnLHist_Sample_DrawflucDraw[i] = NULL;
      lnLHist_Sample_PredflucDraw[i] = NULL;
      // If not NULL it indicates the selection was turned on, so initialise the privates
    } else {
      NominalHist[i] = (TH2Poly*)(Nominal[i]->Clone());
      std::string name = std::string(NominalHist[i]->GetName());
      name = name.substr(0, name.find("_nom"));
      W2NomHist[i] = (TH2Poly*)(NomW2[i]->Clone());

      PosteriorHist[i] = new TH1D*[maxBins[i]];
      w2Hist[i] = new TH1D*[maxBins[i]];

      lnLHist_Mean[i] = (TH2Poly*)(NominalHist[i]->Clone());
      lnLHist_Mean[i]->SetNameTitle((name+"_MeanlnL").c_str(), (name+"_MeanlnL").c_str());
      lnLHist_Mean[i]->Reset("");
      lnLHist_Mean[i]->GetZaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");
      
      lnLHist_Mode[i] = (TH2Poly*)(NominalHist[i]->Clone());
      lnLHist_Mode[i]->SetNameTitle((name+"_ModelnL").c_str(), (name+"_ModelnL").c_str());
      lnLHist_Mode[i]->Reset("");
      lnLHist_Mode[i]->GetZaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");

      lnLHist_Mean_ProjectX[i] = (TH1D*)(DataHist_ProjectX[i]->Clone());
      lnLHist_Mean_ProjectX[i]->SetNameTitle((name+"_MeanlnL_ProjectX").c_str(), (name+"_MeanlnL_ProjectX").c_str()); 
      lnLHist_Mean_ProjectX[i]->Reset("");
      lnLHist_Mean_ProjectX[i]->GetYaxis()->SetTitle("-2lnL_{sample} #times sign(MC-data)");
      
      MeanHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      MeanHist[i]->SetNameTitle((name+"_mean").c_str(), (name+"_mean").c_str());
      MeanHist[i]->Reset("");
      MeanHist[i]->GetZaxis()->SetTitle("Mean");

      ModeHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      ModeHist[i]->SetNameTitle((name+"_mode").c_str(), (name+"_mode").c_str());
      ModeHist[i]->Reset("");
      ModeHist[i]->GetZaxis()->SetTitle("Mode");

      W2MeanHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      W2MeanHist[i]->SetNameTitle((name+"_w2mean").c_str(), (name+"_w2mean").c_str());
      W2MeanHist[i]->Reset("");
      W2MeanHist[i]->GetZaxis()->SetTitle("W2Mean");

      W2ModeHist[i] = (TH2Poly*)(NominalHist[i]->Clone());
      W2ModeHist[i]->SetNameTitle((name+"_w2mode").c_str(), (name+"_w2mode").c_str());
      W2ModeHist[i]->Reset("");
      W2ModeHist[i]->GetZaxis()->SetTitle("W2Mode");

      SumHist[i] = new TH1D((name+"_sum").c_str(),(name+"_sum").c_str(), 100, 1, -1);
      SumHist[i]->GetXaxis()->SetTitle("N_{events}");
      SumHist[i]->GetYaxis()->SetTitle("Counts");
      double Integral = NoOverflowIntegral(DataHist[i]);
      std::stringstream ss;
      ss << Integral;
      SumHist[i]->SetTitle((std::string(SumHist[i]->GetTitle())+"_"+ss.str()).c_str());

      // Declare the lnL histograms
      lnLHist_Mean1D[i] = new TH1D((name+"_MeanlnL1D").c_str(),(name+"_MeanlnL1D").c_str(), 50, 1, -1);
      lnLHist_Mean1D[i]->GetXaxis()->SetTitle("-2LLH (Data, Pred)");
      lnLHist_Mean1D[i]->GetYaxis()->SetTitle("Counts");

      lnLHist_Mode1D[i] = new TH1D((name+"_ModelnL1D").c_str(),(name+"_ModelnL1D").c_str(), 50, 1, -1);
      lnLHist_Mode1D[i]->GetXaxis()->SetTitle("-2LLH (Data, Pred)");
      lnLHist_Mode1D[i]->GetYaxis()->SetTitle("Counts");

      lnLHist_Sample_DrawData[i] = new TH1D((name+"_lnLdrawdata").c_str(),(name+"_lnL").c_str(), 100, 1, -1);
      lnLHist_Sample_DrawData[i]->GetXaxis()->SetTitle("-2LLH (Data, Draw)");
      lnLHist_Sample_DrawData[i]->GetYaxis()->SetTitle("Counts");
      
      lnLHist_Sample_DrawflucDraw[i] = new TH1D((name+"_lnLdrawfluc").c_str(),(name+"_lnL").c_str(), 100, 1, -1);
      lnLHist_Sample_DrawflucDraw[i]->GetXaxis()->SetTitle("-2LLH (Draw Fluc, Draw)");
      lnLHist_Sample_DrawflucDraw[i]->GetYaxis()->SetTitle("Counts");
      
      lnLHist_Sample_PredflucDraw[i] = new TH1D((name+"_lnLpredfluc").c_str(),(name+"_lnL").c_str(), 100, 1, -1);
      lnLHist_Sample_PredflucDraw[i]->GetXaxis()->SetTitle("-2LLH (Pred Fluc, Draw)");
      lnLHist_Sample_PredflucDraw[i]->GetYaxis()->SetTitle("Counts");
      
      //KS: We copy histograms so delete original
      delete Nominal[i];
      delete NomW2[i];
    }
  }
}

// *******************
// Add a throw from the MCMC to the posterior predictive
// The input here is nSamples long
void SampleSummary::AddThrow(std::vector<TH2Poly*> &SampleVector, std::vector<TH2Poly*> &W2Vec, double LLHPenalty, int DrawNumber) {
// *******************

  nThrows++;
  //KS: Only make sense for PosteriorPredictive
  if( !isPriorPredictive )RandomHist->Fill(DrawNumber);

  int size = SampleVector.size();
  if (!CheckPsycheSamples(size)) throw;

  // Push back the throw
  MCVector.push_back(SampleVector);
  LLHPenaltyVector.push_back(LLHPenalty);
  W2MCVector.push_back(W2Vec);

  // Loop over the sameples
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int SampleNum = 0;  SampleNum < PsycheSampleMap.size(); SampleNum++) 
  {
    int PsycheSampleEnum = PsycheSampleMap[SampleNum];
    if (SampleVector[PsycheSampleEnum] == NULL) continue;

    // Initalise the posterior hist
    if (first_pass) 
    {
       const int nXBins = 5000;
      //Initialise TH1D which corresponds to each bin in the sample's th2poly
      std::string name = std::string(SampleVector[PsycheSampleEnum]->GetName());
      for (int i = 0; i <= maxBins[PsycheSampleEnum]; i++)
      {
        std::stringstream ss;
        ss << "_" << i;
        PosteriorHist[PsycheSampleEnum][i] = new TH1D((name+ss.str()).c_str(), (name+ss.str()).c_str(),nXBins, 1, -1);
        w2Hist[PsycheSampleEnum][i] = new TH1D((name+ss.str()+"_w2").c_str(), (name+ss.str()+"_w2").c_str(),nXBins, 1, -1);
      } 
    }
    // Fill the sum histogram with the integral of the sampled distribution
    SumHist[PsycheSampleEnum]->Fill(NoOverflowIntegral(SampleVector[PsycheSampleEnum]));
    // Loop over the distribution and fill the prior/posterior predictive
    for (int i = 1; i <= maxBins[PsycheSampleEnum]; ++i) {
      double Content = SampleVector[PsycheSampleEnum]->GetBinContent(i);
      PosteriorHist[PsycheSampleEnum][i]->Fill(Content);
      double w2 = W2Vec[PsycheSampleEnum]->GetBinContent(i);
      w2Hist[PsycheSampleEnum][i]->Fill(w2);
    } // end bin loop
  } // end nd280 samples loop
  first_pass = false;
} // end AddThrow

// *******************
// Add a throw from the MCMC to the posterior predictive
// The input here is has dimension [nsample][nMaCh3Modes]
void SampleSummary::AddThrowByMode(std::vector<std::vector<TH2Poly*>> &SampleVector_ByMode) {
// *******************

    //KS: This means this is first time
    if(!DoByModePlots)
    {
        std::cout<< "Turning reaction breadkwon mode, brum brum"<<std::endl;
        PosteriorHist_ByMode = new TH1D***[nSamples];
        MeanHist_ByMode = new TH2Poly**[nSamples];
    }
    
    // Loop over the sameples
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (unsigned int SampleNum = 0;  SampleNum < PsycheSampleMap.size(); SampleNum++) 
    {
        int PsycheSampleEnum = PsycheSampleMap[SampleNum];
        if (DataHist[PsycheSampleEnum] == NULL) continue;
        
        if(!DoByModePlots)//Do only first time
        {
            PosteriorHist_ByMode[PsycheSampleEnum] = new TH1D**[nModes+1];
            MeanHist_ByMode[PsycheSampleEnum] = new TH2Poly*[nModes+1];
        }
        for (int j = 0; j < nModes+1; j++)
        {
            if(!DoByModePlots) //Do only first time
            { 
                PosteriorHist_ByMode[PsycheSampleEnum][j] = new TH1D*[maxBins[PsycheSampleEnum]];
                const int nXBins = 5000;

                std::string name = std::string(NominalHist[PsycheSampleEnum]->GetName());
                name = name.substr(0, name.find("_nom"));
                name = name + "_"+MaCh3mode_ToString(MaCh3_Mode(j));
                
                for (int i = 0; i <= maxBins[PsycheSampleEnum]; i++)
                {
                    std::stringstream ss;
                    ss << "_" << i;
                    //Initialise TH1D which corresponds to each bin in the sample's th2poly
                    PosteriorHist_ByMode[PsycheSampleEnum][j][i] = new TH1D((name+ss.str()).c_str(),(name+ss.str()).c_str(),nXBins, 1, -1);
                } 
                MeanHist_ByMode[PsycheSampleEnum][j] = (TH2Poly*)(NominalHist[PsycheSampleEnum]->Clone());
                MeanHist_ByMode[PsycheSampleEnum][j]->SetNameTitle((name+"_mean").c_str(), (name+"_mean").c_str());
                MeanHist_ByMode[PsycheSampleEnum][j]->Reset("");
                MeanHist_ByMode[PsycheSampleEnum][j]->GetZaxis()->SetTitle("Mean");
            }
            // Loop over the distribution and fill the prior/posterior predictive
            for (int i = 1; i <= maxBins[PsycheSampleEnum]; ++i) 
            {
                double Content = SampleVector_ByMode[PsycheSampleEnum][j]->GetBinContent(i);
                PosteriorHist_ByMode[PsycheSampleEnum][j][i]->Fill(Content);
            }
        }
    }
    DoByModePlots = true;
} // end AddThrowByMode

// **********************
void SampleSummary::PrepareTree() {
// **********************

  // Number of good samples
  int nPsycheSamples = PsycheSampleMap.size();

  // The array of doubles we write to the TTree
  // Data vs Draw
  llh_data_draw = new double[nPsycheSamples];
  // Data vs Fluctuated Draw
  llh_data_drawfluc = new double[nPsycheSamples];
  // Data vs Fluctuated Predictive
  llh_data_predfluc = new double[nPsycheSamples];

  // Draw vs Predictive
  llh_draw_pred = new double[nPsycheSamples];

  // Fluctuated Draw vs Predictive
  llh_drawfluc_pred = new double[nPsycheSamples];
  // Fluctuated Draw vs Fluctuated Predictive
  llh_drawfluc_predfluc = new double[nPsycheSamples];
  // Fluctuated Draw vs Draw
  llh_drawfluc_draw = new double[nPsycheSamples];
  // Fluctuated Predictive vs Predictive
  llh_predfluc_pred = new double[nPsycheSamples];
  // Fluctuated Predicitve vs Draw
  llh_predfluc_draw = new double[nPsycheSamples];
   
  // Fluctuated Data vs Draw
  llh_datafluc_draw = new double[nPsycheSamples];
  
  // Data vs Draw for 1D projection
  llh_data_draw_ProjectX = new double[nPsycheSamples];
  llh_drawfluc_draw_ProjectX = new double[nPsycheSamples];
    
  // The output tree we're going to write to
  OutputTree = new TTree("LLH_draws", "LLH_draws");
  SampleNames.resize(nPsycheSamples);
  // Loop over the samples and set the addresses of the variables to write to file
  for (int i = 0; i < nPsycheSamples; ++i) {
    // Get the psyche sample number
    int PsycheSample = PsycheSampleMap[i];
    // Get the name
    std::string SampleName = sampleND->GetSampleName(PsycheSample);
    // Strip out spaces
    while (SampleName.find(" ") != std::string::npos) {
      SampleName.replace(SampleName.find(" "), 1, std::string("_"));
    }
    // Also strip out - signs because it messes up TBranches
    while (SampleName.find("-") != std::string::npos) {
      SampleName.replace(SampleName.find("-"), 1, std::string("_"));
    }

    SampleNames[i] = SampleName;

    OutputTree->Branch((SampleName+"_data_draw").c_str(),     &llh_data_draw[i]);
    OutputTree->Branch((SampleName+"_drawfluc_draw").c_str(), &llh_drawfluc_draw[i]);
    OutputTree->Branch((SampleName+"_predfluc_draw").c_str(), &llh_predfluc_draw[i]);
//    All LLH below are for validation reason but not used for final P-Value
 
    OutputTree->Branch((SampleName+"_data_drawfluc").c_str(), &llh_data_drawfluc[i]);
    OutputTree->Branch((SampleName+"_data_predfluc").c_str(), &llh_data_predfluc[i]);
    OutputTree->Branch((SampleName+"_draw_pred").c_str(),     &llh_draw_pred[i]);
    OutputTree->Branch((SampleName+"_drawfluc_pred").c_str(), &llh_drawfluc_pred[i]);
    OutputTree->Branch((SampleName+"_drawfluc_predfluc").c_str(), &llh_drawfluc_predfluc[i]);
    OutputTree->Branch((SampleName+"_predfluc_pred").c_str(), &llh_predfluc_pred[i]);
    OutputTree->Branch((SampleName+"_datafluc_draw").c_str(), &llh_datafluc_draw[i]);
 
//    All LLH below are used for calcauting P-Value but suign 1D porjections  
    OutputTree->Branch((SampleName+"_data_draw_ProjectX").c_str(), &llh_data_draw_ProjectX[i]);
    OutputTree->Branch((SampleName+"_drawfluc_draw_ProjectX").c_str(), &llh_drawfluc_draw_ProjectX[i]);
  }

  OutputTree->Branch("LLH_Penalty",         &llh_penalty);
  OutputTree->Branch("Total_LLH_Data_Draw", &total_llh_data_draw);
  OutputTree->Branch("Total_LLH_DrawFluc_Draw", &total_llh_drawfluc_draw);
  OutputTree->Branch("Total_LLH_PredFluc_Draw", total_llh_predfluc_draw);
  //    All LLH below are for validation reason but not used for final P-Value

  OutputTree->Branch("Total_LLH_Data_DrawFluc", &total_llh_data_drawfluc);
  OutputTree->Branch("Total_LLH_Data_PredFluc", &total_llh_data_predfluc);
  OutputTree->Branch("Total_LLH_Draw_Pred",     &total_llh_draw_pred);
  OutputTree->Branch("Total_LLH_DrawFluc_Pred", &total_llh_drawfluc_pred);
  OutputTree->Branch("Total_LLH_DrawFluc_PredFluc", &total_llh_drawfluc_predfluc);
  OutputTree->Branch("Total_LLH_PredFluc_Pred", &total_llh_predfluc_pred);
  OutputTree->Branch("Total_LLH_DataFluc_Draw", &total_llh_datafluc_draw);

  OutputTree->Branch("Event_Rate", &event_rate);
//    All LLH below are used for calcauting P-Value but suign 1D porjections 
  OutputTree->Branch("total_llh_data_draw_ProjectX", &total_llh_data_draw_ProjectX);
  OutputTree->Branch("total_llh_drawfluc_draw_ProjectX", &total_llh_drawfluc_draw_ProjectX);
}

// *******************
// Write the contents to the file
void SampleSummary::Write() {
// *******************

  // Make the output file (MakePosterioPredictive call writes to this)
  std::string TempString = OutputName;
  TempString.replace(TempString.find(".root"), 5, std::string("_procsW2.root"));
  Outputfile = new TFile(TempString.c_str(), "RECREATE");
  // Prepare the output tree
  PrepareTree();

  std::cout << "Summarising " << nThrows << " throws..." << std::endl;
  // After all the throws are added finalise the sample
  TStopwatch timer;
  timer.Start();
  MakePredictive();
  timer.Stop();
  std::cout << "Made prior/posterior predictive, it took " << timer.RealTime() << "s" <<", now writing..." << std::endl;

  Outputfile->cd();
  TDirectory **Dir = NULL; Dir = new TDirectory*[nSamples];

  //make fancy event rate histogram
  double DataRate = 0.0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:DataRate)
  #endif
  for (int i = 0; i < nSamples; ++i)
  {
    if (DataHist[i] == NULL) continue;
    DataRate += NoOverflowIntegral(DataHist[i]);
  }
  MakeCutEventRate(EventHist, DataRate);
  
  OutputTree->Write();

  // Make the various distributions
  lnLHist->Write();
  lnLHist_drawfluc->Write();
  lnLHist_drawflucdraw->Write();
  lnLHist_drawdata->Write();
  lnLDrawHist->Write();
  lnLFlucHist->Write();
  //KS: Only avaible for Posterior Predictive
  if(!isPriorPredictive) RandomHist->Write();

  lnLFlucHist_ProjectX->Write();
  
  // Loop over each sample and write to file
  //KS: Multithreading is tempting here but we also write to ROOT file, seprating all LLH and poly projections from write could work well
  int nPsycheSamples = PsycheSampleMap.size();
  for (int j = 0; j < nPsycheSamples; ++j) 
  {
    int i = PsycheSampleMap[j];
    
     // Skip the null histograms
    if (DataHist[i] == NULL || NoOverflowIntegral(DataHist[i]) == 0) continue;

    // Make a new direcotry
    Dir[i] = Outputfile->mkdir(sampleND->GetSampleName(i).c_str());
    Dir[i]->cd();

    // Make the data/MC ratio histogram
    TH2Poly *RatioHistMean = RatioPolys(DataHist[i], MeanHist[i]);
    RatioHistMean->GetZaxis()->SetTitle("Data/Mean");
    TH2Poly *RatioHistMode = RatioPolys(DataHist[i], ModeHist[i]);
    RatioHistMode->GetZaxis()->SetTitle("Data/Mode");
    TH2Poly *RatioHistNom = RatioPolys(DataHist[i], NominalHist[i]);
    RatioHistNom->GetZaxis()->SetTitle("Data/Nom");

    // And the normalised data histogram
    TH2Poly *DataNormHist = NormalisePoly(DataHist[i]);
    // Last true refers to if project along x or y
    TH2Poly *MeanNormHist = NormalisePoly(MeanHist[i]);
    TH2Poly *ModeNormHist = NormalisePoly(ModeHist[i]);
    TH1D *MeanProjectX = ProjectPoly(MeanHist[i], true, i, true);
    TH1D *MeanProjectY = ProjectPoly(MeanHist[i], false, i, true);
    TH1D *ModeProjectX = ProjectPoly(ModeHist[i], true, i, true);
    TH1D *ModeProjectY = ProjectPoly(ModeHist[i], false, i, true);

    TH1D *W2MeanProjectX = ProjectPoly(W2MeanHist[i], true, i);
    TH1D *W2MeanProjectY = ProjectPoly(W2MeanHist[i], false, i);
    TH1D *W2ModeProjectX = ProjectPoly(W2ModeHist[i], true, i);
    TH1D *W2ModeProjectY = ProjectPoly(W2ModeHist[i], false, i);

    TH2Poly *NomNormHist = NormalisePoly(NominalHist[i]);
    TH1D *NomProjectX = ProjectPoly(NominalHist[i], true, i);
    TH1D *NomProjectY = ProjectPoly(NominalHist[i], false, i);

    TH1D *W2NomProjectX = ProjectPoly(W2NomHist[i], true, i);
    TH1D *W2NomProjectY = ProjectPoly(W2NomHist[i], false, i);

    // Same for the TH2Ds
    CalcLLH(DataHist[i], NominalHist[i], W2NomHist[i]);
    CalcLLH(DataHist[i], MeanHist[i], W2MeanHist[i]);
    CalcLLH(DataHist[i], ModeHist[i], W2ModeHist[i]);

    // Calculate the log likelihood for the 1D dists
    // Sets the title of the second TH1D to the -2LLH
    CalcLLH(DataHist_ProjectX[i], NomProjectX, W2NomProjectX);
    CalcLLH(DataHist_ProjectX[i], MeanProjectX, W2MeanProjectX);
    CalcLLH(DataHist_ProjectX[i], ModeProjectX, W2ModeProjectX);
    CalcLLH(DataHist_ProjectY[i], NomProjectY, W2NomProjectY);
    CalcLLH(DataHist_ProjectY[i], MeanProjectY, W2MeanProjectY);
    CalcLLH(DataHist_ProjectY[i], ModeProjectY, W2ModeProjectY);
    
    //Make fancy event rate histogram
    MakeCutEventRate(SumHist[i], NoOverflowIntegral(DataHist[i]));
        
    OutputTree->Draw((SampleNames[j]+"_data_draw:"+SampleNames[j]+"_drawfluc_draw>>htemp").c_str());
    TH2D *TempHistogram = (TH2D*)((gDirectory->Get("htemp"))->Clone());
    TempHistogram->GetXaxis()->SetTitle("-2LLH(Draw Fluc, Draw)");
    TempHistogram->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram->SetNameTitle((SampleNames[j]+"_drawfluc_draw").c_str(), (SampleNames[j]+"_drawfluc_draw").c_str());
    MakeCutLLH2D(TempHistogram);
    TempHistogram->Write();
    delete TempHistogram;

    // Also write the 2D histograms for the p-value
    OutputTree->Draw((SampleNames[j]+"_data_draw:"+SampleNames[j]+"_predfluc_draw>>htemp2").c_str());
    TH2D *TempHistogram2 = (TH2D*)((gDirectory->Get("htemp2"))->Clone());
    TempHistogram2->GetXaxis()->SetTitle("-2LLH(Pred Fluc, Draw)");
    TempHistogram2->GetYaxis()->SetTitle("-2LLH(Data, Draw)");
    TempHistogram2->SetNameTitle((SampleNames[j]+"_predfluc_draw").c_str(), (SampleNames[j]+"_predfluc_draw").c_str());
    MakeCutLLH2D(TempHistogram2);
    TempHistogram2->Write();
    delete TempHistogram2;
   
    // finally p-value for 1D projection
    OutputTree->Draw((SampleNames[j]+"_data_draw_ProjectX:"+SampleNames[j]+"_drawfluc_draw_ProjectX>>htemp3").c_str());
    TH2D *TempHistogram3 = (TH2D*)((gDirectory->Get("htemp3"))->Clone());
    KinematicTypes Type1, Type2;
    sampleND->GetKinVars(i, Type1, Type2);
    TempHistogram3->GetXaxis()->SetTitle(("-2LLH_{Draw Fluc, Draw} for " + Kinematic_ToString(Type1)).c_str());
    TempHistogram3->GetYaxis()->SetTitle(("-2LLH_{Data, Draw} for " + Kinematic_ToString(Type1)).c_str());
    TempHistogram3->SetNameTitle((SampleNames[j]+"_drawfluc_draw_ProjectX").c_str(), (SampleNames[j]+"_drawfluc_draw_ProjectX").c_str());
    MakeCutLLH2D(TempHistogram3);
    TempHistogram3->Write();
    delete TempHistogram3;
    
    // Write the Histgorams to each folder
    DataHist[i]->Write();
    NominalHist[i]->Write();
    MeanHist[i]->Write();
    ModeHist[i]->Write();
    RatioHistMean->Write();
    RatioHistMode->Write();
    RatioHistNom->Write();

    W2NomHist[i]->Write();
    W2MeanHist[i]->Write();
    W2ModeHist[i]->Write();

    DataNormHist->Write();
    NomNormHist->Write();
    MeanNormHist->Write();
    ModeNormHist->Write();

    DataHist_ProjectX[i]->Write();
    NomProjectX->Write();
    MeanProjectX->Write();
    ModeProjectX->Write();

    DataHist_ProjectY[i]->Write();
    NomProjectY->Write();
    MeanProjectY->Write();
    ModeProjectY->Write();

    W2NomProjectX->Write();
    W2MeanProjectX->Write();
    W2ModeProjectX->Write();

    W2NomProjectY->Write();
    W2MeanProjectY->Write();
    W2ModeProjectY->Write();
    
    //KS: This will dump lots of hists, use it only for debuging 
    /*NOTE
    for (int b = 1; b <= maxBins[i]; ++b) 
    {
      PosteriorHist[i][b]->Write();
      //w2Hist[i][b]->Write(); //This isn't useful check only in desperation
    } 
    */
    lnLHist_Mean[i]->Write();
    lnLHist_Mode[i]->Write();
    
    lnLHist_Mean_ProjectX[i]->Write();
    
    lnLHist_Mean1D[i]->Write();
    lnLHist_Mode1D[i]->Write();
    
    MakeCutLLH1D(lnLHist_Sample_DrawData[i], GetLLH(DataHist[i], MeanHist[i], W2MeanHist[i]));
    lnLHist_Sample_DrawData[i]->Write();
    MakeCutLLH1D(lnLHist_Sample_DrawflucDraw[i], GetLLH(DataHist[i], MeanHist[i], W2MeanHist[i]));
    lnLHist_Sample_DrawflucDraw[i]->Write();
    MakeCutLLH1D(lnLHist_Sample_PredflucDraw[i], GetLLH(DataHist[i], MeanHist[i], W2MeanHist[i]));
    lnLHist_Sample_PredflucDraw[i]->Write();
    
    if(DoByModePlots)
    {
        for (int j = 0; j < nModes+1; j++)
        {
            MeanHist_ByMode[i][j]->Write();
            TH1D *MeanProjectX_ByMode = ProjectPoly(MeanHist_ByMode[i][j], true, i, true);
            TH1D *MeanProjectY_ByMode = ProjectPoly(MeanHist_ByMode[i][j], false, i, true);
            MeanProjectX_ByMode->Write();
            MeanProjectY_ByMode->Write();
            //KS: This will dump lots of hists, use it only for debuging 
            /*NOTE
            for (int b = 1; b <= maxBins[i]; ++b) 
            {
                PosteriorHist_ByMode[i][j][b]->Write();
            } 
            */
            
            delete MeanProjectX_ByMode;
            delete MeanProjectY_ByMode;
        } // End loop over bins
    }
    
    // Delete temporary objects
    delete RatioHistMean;
    delete RatioHistMode;
    delete RatioHistNom;

    delete DataNormHist;
    delete MeanNormHist;
    delete ModeNormHist;
    delete NomNormHist;

    delete DataHist_ProjectX[i];
    delete MeanProjectX;
    delete ModeProjectX;
    delete NomProjectX;

    delete DataHist_ProjectY[i];
    delete MeanProjectY;
    delete ModeProjectY;
    delete NomProjectY;

    delete W2NomProjectX;
    delete W2MeanProjectX;
    delete W2ModeProjectX;
  
    delete W2NomProjectY;
    delete W2MeanProjectY;
    delete W2ModeProjectY;
    
    std::cout << std::endl;
  }
  delete[] DataHist_ProjectX;
  delete[] DataHist_ProjectY;

  std::cout << "Wrote to " << Outputfile->GetName() << std::endl;
  Outputfile->Close();
  delete Outputfile;
}

// *******************
// Make the posterior predictive distributions: fit Poissons etc
void SampleSummary::MakePredictive() {
// *******************

  // First make the projection on the z axis of the TH3D* for every pmu cosmu bin
  llh_total = 0.0;

  // Loop over the psyche samples
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:llh_total)
  #endif
  for (unsigned int SampleNum = 0;  SampleNum < PsycheSampleMap.size(); SampleNum++) 
  {
    int PsycheSampleEnum = PsycheSampleMap[SampleNum];
    // Skip disabled samples
    if (DataHist[PsycheSampleEnum] == NULL || NoOverflowIntegral(DataHist[PsycheSampleEnum]) == 0) continue;

    // Count the -2LLH for each histogram
    double negLogL_Mean = 0.0;
    double negLogL_Mode = 0.0;

    // Loop over each pmu cosmu bin
    for (int j = 1; j < maxBins[PsycheSampleEnum]+1; ++j) 
    {
      //Get PolyBin
      TH2PolyBin* bin = (TH2PolyBin*)DataHist[PsycheSampleEnum]->GetBins()->At(j-1)->Clone();
      
      // Just make a little fancy name
      std::string name = PosteriorHist[PsycheSampleEnum][j]->GetName();
      std::stringstream ss2;
      ss2 << name << "_";
      ss2 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
      ss2 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";
      PosteriorHist[PsycheSampleEnum][j]->SetNameTitle(ss2.str().c_str(), ss2.str().c_str());
      w2Hist[PsycheSampleEnum][j]->SetNameTitle(("w2_"+ss2.str()).c_str(), ("w2_"+ss2.str()).c_str());
      
      TH1D *Projection = (TH1D*)PosteriorHist[PsycheSampleEnum][j]->Clone();
      TH1D *W2Projection = (TH1D*)w2Hist[PsycheSampleEnum][j]->Clone();

      // Data content for the j,kth bin
      double nData = DataHist[PsycheSampleEnum]->GetBinContent(j);
      
      // Get the mean for this projection for all the samples
      // This is the mean prediction for this given j,k bin
      double nMean = Projection->GetMean();
      double nMeanError = Projection->GetRMS();
      double nMode = Projection->GetBinCenter(Projection->GetMaximumBin());
      double nModeError = GetModeError(Projection);
      
      double nW2Mean = W2Projection->GetMean();
      double nW2Mode = W2Projection->GetBinCenter(W2Projection->GetMaximumBin());
      
      double TempLLH_Mean = 0.0;
      double TempLLH_Mode = 0.0;
      
      //KS:Get LLH contibution getTestStatLLH can calcualte Barlow Beeston/IceCube or Poisson
      TempLLH_Mean = sampleND->getTestStatLLH(nData, nMean, nW2Mean);
      TempLLH_Mode = sampleND->getTestStatLLH(nData, nMode, nW2Mode);

      // Increment -2LLH
      //KS: do times 2 because banff reports chi2
      negLogL_Mean += 2*TempLLH_Mean;
      negLogL_Mode += 2*TempLLH_Mode;
      
      // Set the content and error to the mean in the bin
      MeanHist[PsycheSampleEnum]->SetBinContent(j, MeanHist[PsycheSampleEnum]->GetBinContent(j)+nMean);
      MeanHist[PsycheSampleEnum]->SetBinError(j, nMeanError);
      // Set the content to the mode in the bin
      ModeHist[PsycheSampleEnum]->SetBinContent(j, ModeHist[PsycheSampleEnum]->GetBinContent(j)+nMode);
      ModeHist[PsycheSampleEnum]->SetBinError(j, nModeError);
      // Set the content to the mean in the bin
      W2MeanHist[PsycheSampleEnum]->SetBinContent(j, W2MeanHist[PsycheSampleEnum]->GetBinContent(j)+nW2Mean);
      // Set the content to the mode in the bin
      W2ModeHist[PsycheSampleEnum]->SetBinContent(j, W2ModeHist[PsycheSampleEnum]->GetBinContent(j)+nW2Mode);
      
      // Set the mean and average LLH for this given bin
      // Can use these hists to see where the largest -2LLH hists come from
      lnLHist_Mean[PsycheSampleEnum]->SetBinContent(j, 2.0*TempLLH_Mean);
      lnLHist_Mode[PsycheSampleEnum]->SetBinContent(j, 2.0*TempLLH_Mode);
      
      lnLHist_Mean1D[PsycheSampleEnum]->Fill(2.0*TempLLH_Mean);
      lnLHist_Mode1D[PsycheSampleEnum]->Fill(2.0*TempLLH_Mode);
      
      delete Projection;
      delete W2Projection;
      delete bin;
    } // End loop over bins
    
    //KS:: Might consider caching it as we use it once agian much later
    TH1D *MeanProjectX = ProjectPoly(MeanHist[PsycheSampleEnum], true, PsycheSampleEnum, true);
    TH1D *W2MeanProjectX = ProjectPoly(W2MeanHist[PsycheSampleEnum], true, PsycheSampleEnum);
    // Loop over each pmu bin for 1D projection
    for (int j = 1; j <= lnLHist_Mean_ProjectX[PsycheSampleEnum]->GetXaxis()->GetNbins(); ++j) 
    {
        // Data content for the j,kth bin
        double nData = DataHist_ProjectX[PsycheSampleEnum]->GetBinContent(j);
        double nMean = MeanProjectX->GetBinContent(j);
        double nW2Mean = W2MeanProjectX->GetBinContent(j);

        double TempLLH_Mean = 0.0;
        TempLLH_Mean = sampleND->getTestStatLLH(nData, nMean, nW2Mean);
        
        //KS: do times 2 because banff reports chi2
        lnLHist_Mean_ProjectX[PsycheSampleEnum]->SetBinContent(j, 2.0*TempLLH_Mean);
    }// End loop over  bins
    
    
    if(DoByModePlots)
    {
        for (int j = 0; j < nModes+1; j++)
        {
            // Loop over each pmu cosmu bin
            for (int i = 1; i < maxBins[PsycheSampleEnum]+1; ++i) 
            {
                //Get PolyBin
                TH2PolyBin* bin = (TH2PolyBin*)DataHist[PsycheSampleEnum]->GetBins()->At(i-1)->Clone();
                
                // Just make a little fancy name
                std::string name = PosteriorHist_ByMode[PsycheSampleEnum][j][i]->GetName();
                std::stringstream ss2;
                ss2 << name << "_";
                ss2 << "p_{#mu} (" << bin->GetXMin() << "-" << bin->GetXMax() << ")";
                ss2 << " cos#theta_{#mu} (" << bin->GetYMin() << "-" << bin->GetYMax() << ")";
                PosteriorHist_ByMode[PsycheSampleEnum][j][i]->SetNameTitle(ss2.str().c_str(), ss2.str().c_str());
                
                // Make the posterior/prior predictive projection on z
                // The z axis of Predictive is the bin content
                // Essentailly zooming in on one bin and looking at the mean and mode of that bin
                TH1D *Projection = (TH1D*)PosteriorHist_ByMode[PsycheSampleEnum][j][i]->Clone();

                // Get the mean for this projection for all the samples
                double nMean = Projection->GetMean();
                double nMeanError = Projection->GetRMS();

                // Set the content and error to the mean in the bin
                MeanHist_ByMode[PsycheSampleEnum][j]->SetBinContent(i, MeanHist_ByMode[PsycheSampleEnum][j]->GetBinContent(i)+nMean);
                MeanHist_ByMode[PsycheSampleEnum][j]->SetBinError(i, nMeanError);

                delete Projection;
                delete bin;
            } // End loop over bins
        }
    }
    llh_total += negLogL_Mean;
    
    delete MeanProjectX;
    delete W2MeanProjectX;
  } // End loop over samples

  // Now we have our posterior predictive histogram and it's LLH
  std::cout << "Prior/Posterior predictive LLH mean (sample only) = " << llh_total << std::endl;
  std::stringstream ss;
  ss << llh_total;
  lnLHist->SetTitle((std::string(lnLHist->GetTitle())+"_"+ss.str()).c_str());

  // Now make the fluctuated hists of the MeanHist and ModeHist
  MakeChi2Hists();

  // Get the 1D LLH dists
  MakeCutLLH();

} // End MakePredictive() function

// *******************
// Make the fluctuated histograms (2D and 1D) for the chi2s
// Essentially taking the MCMC draws and calculating their LLH to the Posterior predicitve distribution
// And additionally taking the data histogram and calculating the LLH to the predictive distribution
// Additionally we calculate the chi2 of the draws (fluctuated) of  the MC with the prior/posterior predictive and plot it vs the chi2 from the draws of MCMC and the data
void SampleSummary::MakeChi2Hists() {
// *******************

  std::cout << "Making the chi2 histograms..." << std::endl;
  // Have this to signify if we're doing the first pass
  first_pass = true;

  double AveragePenalty = 0;

  // Vectors to hold exact LLH
  std::vector<double> LLH_PredFluc_V;
  std::vector<double> LLH_DataDraw_V;
  std::vector<double> LLH_DrawFlucDraw_V;
    
  LLH_PredFluc_V.resize(nThrows);
  LLH_DataDraw_V.resize(nThrows);
  LLH_DrawFlucDraw_V.resize(nThrows);
    
  // Loop over the draws
  // Should look into multi-threading this. Would require temporary THxx structures from arrays
  //KS: Update above would be ideal but currently we loop over samples (see loop below) which isn't as effcient as loop over throws but much much easier to implement
  for (int i = 0; i < nThrows; ++i) 
  {
    if (i % (nThrows/10) == 0) {
      std::cout << "   On throw " << i << "/" << nThrows << " (" << int(double(i)*100.0/double(nThrows)) << "%)"<< std::endl;
    }

    // Set the total LLH to zero to initialise
    total_llh_data_draw = 0.0;
    total_llh_drawfluc_draw = 0.0;
    total_llh_predfluc_draw = 0.0;
    
    total_llh_data_drawfluc = 0.0;
    total_llh_data_predfluc = 0.0;
    total_llh_draw_pred = 0.0;
    total_llh_drawfluc_pred = 0.0;
    total_llh_drawfluc_predfluc = 0.0;
    total_llh_predfluc_pred = 0.0;
    total_llh_datafluc_draw = 0.0;

    event_rate = 0.0;

    total_llh_data_draw_ProjectX = 0.0;
    total_llh_drawfluc_draw_ProjectX = 0.0;
    // Save the double that gets written to file
    llh_penalty = LLHPenaltyVector[i];
    AveragePenalty += llh_penalty;

    #ifdef MULTITHREAD 
    //KS: might be most obscure OMP reduction I have ever made...
    #pragma omp parallel for reduction(+:total_llh_data_draw, total_llh_drawfluc_draw, total_llh_predfluc_draw, total_llh_data_drawfluc, total_llh_data_predfluc, total_llh_draw_pred, total_llh_drawfluc_pred, total_llh_drawfluc_predfluc, total_llh_predfluc_pred, total_llh_datafluc_draw, event_rate, total_llh_data_draw_ProjectX, total_llh_drawfluc_draw_ProjectX)
    #endif
    // Loop over the samples
    for (unsigned int SampleNum = 0;  SampleNum < PsycheSampleMap.size(); SampleNum++) 
    {
      int PsycheSampleEnum = PsycheSampleMap[SampleNum];
      // Get the ith draw for the jth psyche sample
      TH2Poly *DrawHist = (TH2Poly*)(MCVector[i][PsycheSampleEnum])->Clone();
      TH2Poly *DrawW2Hist = (TH2Poly*)(W2MCVector[i][PsycheSampleEnum])->Clone();
      // Skip empty samples
      if (DrawHist == NULL) continue;

      // Add LLH penalties from the systematics to the LLH that use the drawn histogram
      // Data vs Draw
      llh_data_draw[SampleNum] = llh_penalty;
      // Fluctuated Draw vs Draw
      llh_drawfluc_draw[SampleNum] = llh_penalty;
      // Fluctuated Predicitve vs Draw
      llh_predfluc_draw[SampleNum] = llh_penalty;
      
       // Data vs Fluctuated Draw
      llh_data_drawfluc[SampleNum] = llh_penalty;
      // Draw vs Predictive
      llh_draw_pred[SampleNum] = llh_penalty;
      // Fluctuated Draw vs Predictive
      llh_drawfluc_pred[SampleNum] = llh_penalty;
      // Fluctuated Draw vs Fluctuated Predictive
      llh_drawfluc_predfluc[SampleNum] = llh_penalty;
      // Fluctuated Data vs Draw
      llh_datafluc_draw[SampleNum] = llh_penalty;
      
      //Some LLH for 1D projections
      llh_data_draw_ProjectX[SampleNum] = llh_penalty;
      llh_drawfluc_draw_ProjectX[SampleNum] = llh_penalty;
      
      //Other get 0 penalty term
      // Fluctuated Predictive vs Predictive
      llh_predfluc_pred[SampleNum] = 0.0;
      // Data vs Fluctuated Predictive
      llh_data_predfluc[SampleNum] = 0.0;

      // Make the Poisson fluctuated hist
      TH2Poly *FluctHist = MakeFluctuatedHistogram(MeanHist[PsycheSampleEnum]);
      // Also Poisson fluctuate the drawn MCMC hist
      TH2Poly *FluctDrawHist = MakeFluctuatedHistogram(DrawHist);
      // Finally Poisson fluctuate the data hsitogram
      TH2Poly *DataFlucHist = MakeFluctuatedHistogram(DataHist[PsycheSampleEnum]);

      // Likelihood between the drawn histogram and the data
      double DataDrawLLH = GetLLH(DataHist[PsycheSampleEnum], DrawHist, DrawW2Hist);
      llh_data_draw[SampleNum] += DataDrawLLH;
      total_llh_data_draw += DataDrawLLH;
      
      // Likelihood between drawn histogram and fluctuated drawn histogram
      double DrawFlucDrawLLH = GetLLH(FluctDrawHist, DrawHist, DrawW2Hist);
      llh_drawfluc_draw[SampleNum] += DrawFlucDrawLLH;
      total_llh_drawfluc_draw += DrawFlucDrawLLH;
      
      // Likelihood between drawn histogram and fluctuated posterior predictive distribution
      double PredFlucDrawLLH = GetLLH(FluctHist, DrawHist, DrawW2Hist);
      llh_predfluc_draw[SampleNum] += PredFlucDrawLLH;
      total_llh_predfluc_draw += PredFlucDrawLLH;

//    All LLH below are for validation reason but not used for final P-Value

      // Likelihood between the fluctuated drawn histogram and the data
      double DataDrawFlucLLH = GetLLH(DataHist[PsycheSampleEnum], FluctDrawHist, DrawW2Hist);
      llh_data_drawfluc[SampleNum] += DataDrawFlucLLH;
      total_llh_data_drawfluc += DataDrawFlucLLH;

      // Likelihood between the drawn histogram and the data
      double DataPredFlucLLH = GetLLH(DataHist[PsycheSampleEnum], FluctHist, W2MeanHist[PsycheSampleEnum]);
      llh_data_predfluc[SampleNum] += DataPredFlucLLH;
      total_llh_data_predfluc += DataPredFlucLLH;

      // Likelihood between the drawn hist and the Posterior Predictive
      double DrawPredLLH = GetLLH(DrawHist, MeanHist[PsycheSampleEnum], W2MeanHist[PsycheSampleEnum]);
      llh_draw_pred[SampleNum] += DrawPredLLH;
      total_llh_draw_pred += DrawPredLLH;

      // Likelihood between fluctuated drawn and predictive
      double DrawFlucPredLLH = GetLLH(FluctDrawHist, MeanHist[PsycheSampleEnum], W2MeanHist[PsycheSampleEnum]);
      llh_drawfluc_pred[SampleNum]  += DrawFlucPredLLH;
      total_llh_drawfluc_pred  += DrawFlucPredLLH;

      // Likelihood between drawn histogram and fluctuated drawn histogram
      double DrawFlucPredFlucLLH = GetLLH(FluctDrawHist, FluctHist, W2MeanHist[PsycheSampleEnum]);
      llh_drawfluc_predfluc[SampleNum] += DrawFlucPredFlucLLH;
      total_llh_drawfluc_predfluc += DrawFlucPredFlucLLH;

      // Likelihood between the fluctuated drawn histogram and the posterior predictive
      double PredFlucPredLLH = GetLLH(FluctHist, MeanHist[PsycheSampleEnum], W2MeanHist[PsycheSampleEnum]);
      llh_predfluc_pred[SampleNum] += PredFlucPredLLH;
      total_llh_predfluc_pred += PredFlucPredLLH;

      // Likelihood between fluctuated data histogram and drawn histogram 
      double DataFlucDrawLLH = GetLLH(DataFlucHist, DrawHist, DrawW2Hist);
      llh_datafluc_draw[SampleNum] += DataFlucDrawLLH;
      total_llh_datafluc_draw += DataFlucDrawLLH;
      
      event_rate += NoOverflowIntegral(DrawHist);
      
      lnLHist_Sample_DrawData[PsycheSampleEnum]->Fill(DataDrawLLH);
      lnLHist_Sample_DrawflucDraw[PsycheSampleEnum]->Fill(DrawFlucDrawLLH);
      lnLHist_Sample_PredflucDraw[PsycheSampleEnum]->Fill(PredFlucDrawLLH);
      
      
//    At the end we leave LLH for 1D projections  
      TH1D *DrawHistProjectX = ProjectPoly(DrawHist, true, PsycheSampleEnum);
      TH1D *DrawW2HistProjectX = ProjectPoly(DrawW2Hist, true, PsycheSampleEnum);

      TH1D *FluctDrawHistProjectX = MakeFluctuatedHistogram(DrawHistProjectX);
      
      // Likelihood between the drawn histogram and the data for muon momentum
      double DataDrawLLH_ProjectX = GetLLH(DataHist_ProjectX[PsycheSampleEnum], DrawHistProjectX, DrawW2HistProjectX);
      llh_data_draw_ProjectX[SampleNum] += DataDrawLLH_ProjectX;
      total_llh_data_draw_ProjectX += DataDrawLLH_ProjectX;
      
      double DrawFlucDrawLLH_ProjectX = GetLLH(FluctDrawHistProjectX, DrawHistProjectX, DrawW2HistProjectX);
      llh_drawfluc_draw_ProjectX[SampleNum] += DrawFlucDrawLLH_ProjectX;
      total_llh_drawfluc_draw_ProjectX += DrawFlucDrawLLH_ProjectX;
      
      
      // Delete the fluctuated histograms
      delete FluctHist;
      delete FluctDrawHist;
      delete DataFlucHist;
      delete DrawHist;
      delete DrawW2Hist;
      
      delete DrawHistProjectX;
      delete DrawW2HistProjectX;
      delete FluctDrawHistProjectX;
    } // End loop over samples (still looping throws)

    // Add LLH penalties from the systematics to the LLH that use the drawn histogram
    total_llh_data_draw     += llh_penalty;
    total_llh_drawfluc_draw += llh_penalty;
    total_llh_predfluc_draw += llh_penalty;
    
    total_llh_data_drawfluc += llh_penalty;
    total_llh_draw_pred     += llh_penalty;
    total_llh_drawfluc_pred += llh_penalty;
    total_llh_drawfluc_predfluc += llh_penalty;

    total_llh_data_draw_ProjectX += llh_penalty;
    total_llh_drawfluc_draw_ProjectX += llh_penalty;

    lnLHist->Fill(total_llh_predfluc_pred);
    lnLHist_drawdata->Fill(total_llh_data_draw);
    lnLHist_drawfluc->Fill(total_llh_predfluc_draw);
    lnLHist_drawflucdraw->Fill(total_llh_drawfluc_draw);

    lnLDrawHist->Fill(total_llh_predfluc_draw, total_llh_data_draw);
    lnLFlucHist->Fill(total_llh_drawfluc_draw, total_llh_data_draw);

    EventHist->Fill(event_rate);

    lnLFlucHist_ProjectX->Fill(total_llh_drawfluc_draw_ProjectX, total_llh_data_draw_ProjectX);
    
    // Also save to arrays to make sure we have the utmost super accuracy
    LLH_PredFluc_V[i] = total_llh_predfluc_draw;
    LLH_DataDraw_V[i] = total_llh_data_draw;
    LLH_DrawFlucDraw_V[i] = total_llh_drawfluc_draw;

    // Write to the output tree
    OutputTree->Fill();
  } // End loop over throws

  AveragePenalty = AveragePenalty/double(nThrows);
  std::cout << "Average LLH penalty over toys is " << AveragePenalty << std::endl;

  // Calculate exact p-value instead of binned
  unsigned int Accept_PredFluc = 0;
  unsigned int Accept_DrawFluc = 0;
  for (int i = 0; i < nThrows; ++i) 
  {
    if (LLH_DataDraw_V[i] > LLH_DrawFlucDraw_V[i]) Accept_DrawFluc++;
    if (LLH_DataDraw_V[i] > LLH_PredFluc_V[i]) Accept_PredFluc++;
  }
  double pvalue_DrawFluc = double(Accept_DrawFluc)/double(nThrows);
  double pvalue_PredFluc = double(Accept_PredFluc)/double(nThrows);

  std::cout << "Calculated exact p-value using Fluctuation of Draw: "       << pvalue_DrawFluc << std::endl;
  std::cout << "Calculated exact p-value using Fluctuation of Prediction: " << pvalue_PredFluc << std::endl;
}

// *******************
// Make the cut LLH histogram
void SampleSummary::MakeCutLLH() {
// *******************
  Outputfile->cd();
  MakeCutLLH1D(lnLHist);
  MakeCutLLH1D(lnLHist_drawfluc);
  MakeCutLLH1D(lnLHist_drawdata);
  MakeCutLLH1D(lnLHist_drawflucdraw);
  
  MakeCutLLH2D(lnLDrawHist);
  MakeCutLLH2D(lnLFlucHist);
  MakeCutLLH2D(lnLFlucHist_ProjectX);
}

// ****************
// Make the 1D cut distribution and give the 1D p-value
void SampleSummary::MakeCutLLH1D(TH1D *Histogram, double llh_ref) {
// ****************
  double TotalIntegral = Histogram->Integral();
  double Above = 0.0;
  // Get the LLH reference from total llh or some reference histogram
  double llh_reference = 0.0;
  if (llh_ref >= 0) {
    llh_reference = llh_ref;
  } else {
    llh_reference = llh_total;
  }
  for (int i = 0; i < Histogram->GetXaxis()->GetNbins(); ++i) {
    double xvalue = Histogram->GetBinCenter(i+1);
    if (xvalue >= llh_reference) {
      Above += Histogram->GetBinContent(i+1);
    }
  }
  double pvalue = Above/TotalIntegral;
  std::stringstream ss;
  ss << int(Above) << "/" << int(TotalIntegral) << "=" << pvalue;
  Histogram->SetTitle((std::string(Histogram->GetTitle())+"_"+ss.str()).c_str());

  // Write a TCanvas and make a line and a filled histogram
  TLine *TempLine = new TLine(llh_reference , Histogram->GetMinimum(), llh_reference, Histogram->GetMaximum());
  TempLine->SetLineColor(kBlack);
  TempLine->SetLineWidth(2);

  // Make the fill histogram
  TH1D *TempHistogram = (TH1D*)(Histogram->Clone());
  TempHistogram->SetFillStyle(1001);
  TempHistogram->SetFillColor(kRed);
  for (int i = 0; i < TempHistogram->GetNbinsX(); ++i) 
  {
    if (TempHistogram->GetBinCenter(i+1) < llh_reference) 
    {
      TempHistogram->SetBinContent(i+1, 0.0);
    }
  }

  TLegend *Legend = new TLegend(0.6, 0.6, 0.9, 0.9);
  Legend->SetFillColor(0);
  Legend->SetFillStyle(0);
  Legend->SetLineWidth(0);
  Legend->SetLineColor(0);
  Legend->AddEntry(TempLine, Form("Reference LLH, %.0f, p-value=%.2f", llh_reference, pvalue), "l");
  Legend->AddEntry(Histogram, Form("LLH, #mu=%.1f#pm%.1f", Histogram->GetMean(), Histogram->GetRMS()), "l");
  std::string Title = Histogram->GetName();
  Title += "_canv";
  TCanvas *TempCanvas = new TCanvas(Title.c_str(), Title.c_str(), 1024, 1024);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();
  Histogram->Draw();
  TempHistogram->Draw("same");
  TempLine->Draw("same");
  Legend->Draw("same");

  TempCanvas->Write();

  delete TempLine;
  delete TempHistogram;
  delete TempCanvas;
  delete Legend;
}


// ****************
// Make the 2D cut distribution and give the 2D p-value
void SampleSummary::MakeCutLLH2D(TH2D *Histogram) {
// ****************

  double TotalIntegral = Histogram->Integral();
  // Count how many fills are above y=x axis
  // This is the 2D p-value
  double Above = 0.0;
  for (int i = 0; i < Histogram->GetXaxis()->GetNbins(); ++i) {
    double xvalue = Histogram->GetXaxis()->GetBinCenter(i+1);
    for (int j = 0; j < Histogram->GetYaxis()->GetNbins(); ++j) {
      double yvalue = Histogram->GetYaxis()->GetBinCenter(j+1);
      // We're only interested in being _ABOVE_ the y=x axis
      if (xvalue >= yvalue) {
        Above += Histogram->GetBinContent(i+1, j+1);
      }
    }
  }
  double pvalue = Above/TotalIntegral;
  std::stringstream ss;
  ss << int(Above) << "/" << int(TotalIntegral) << "=" << pvalue;
  Histogram->SetTitle((std::string(Histogram->GetTitle())+"_"+ss.str()).c_str());

  // Now add the TLine going diagonally
  double minimum = Histogram->GetXaxis()->GetBinLowEdge(1);
  if (Histogram->GetYaxis()->GetBinLowEdge(1) > minimum) {
    minimum = Histogram->GetYaxis()->GetBinLowEdge(1);
  }
  double maximum = Histogram->GetXaxis()->GetBinLowEdge(Histogram->GetXaxis()->GetNbins());
  if (Histogram->GetYaxis()->GetBinLowEdge(Histogram->GetYaxis()->GetNbins()) < maximum) {
    maximum = Histogram->GetYaxis()->GetBinLowEdge(Histogram->GetYaxis()->GetNbins());
    //KS: Extend by bin width to perfectly fit canvas
    maximum += Histogram->GetYaxis()->GetBinWidth(Histogram->GetYaxis()->GetNbins());
  }
  else maximum += Histogram->GetXaxis()->GetBinWidth(Histogram->GetXaxis()->GetNbins());
  TLine *TempLine = new TLine(minimum, minimum, maximum, maximum);
  TempLine->SetLineColor(kRed);
  TempLine->SetLineWidth(2);

  std::string Title = Histogram->GetName();
  Title += "_canv";
  TCanvas *TempCanvas = new TCanvas(Title.c_str(), Title.c_str(), 1024, 1024);
  TempCanvas->SetGridx();
  TempCanvas->SetGridy();
  TempCanvas->cd();
  gStyle->SetPalette(81);
  Histogram->SetMinimum(-0.01);
  Histogram->Draw("colz");
  TempLine->Draw("same");

  TempCanvas->Write();
  delete TempLine;
  delete TempCanvas;
}

// ****************
// Make the 1D Event Rate Hist
void SampleSummary::MakeCutEventRate(TH1D *Histogram, double DataRate) {
// ****************
    // For the event rate histogram add a TLine to the data rate
    TLine *TempLine = new TLine(DataRate, Histogram->GetMinimum(), DataRate, Histogram->GetMaximum());
    TempLine->SetLineColor(kRed);
    TempLine->SetLineWidth(2);
    // Also fit a Gaussian because why not?
    TF1 *Fitter = new TF1("Fit", "gaus", Histogram->GetBinLowEdge(1), Histogram->GetBinLowEdge(Histogram->GetNbinsX()+1));
    Histogram->Fit(Fitter, "RQ");
    Fitter->SetLineColor(kRed-5);
    // Calculate a p-value
    double Above = 0.0;
    for (int z = 0; z < Histogram->GetNbinsX(); ++z) {
      double xvalue = Histogram->GetBinCenter(z+1);
      if (xvalue >= DataRate) {
        Above += Histogram->GetBinContent(z+1);
      }
    }
    double pvalue = Above/Histogram->Integral();
    TLegend *Legend = new TLegend(0.4, 0.75, 0.98, 0.90);
    Legend->SetFillColor(0);
    Legend->SetFillStyle(0);
    Legend->SetLineWidth(0);
    Legend->SetLineColor(0);
    Legend->AddEntry(TempLine, Form("Data, %.0f, p-value=%.2f", DataRate, pvalue), "l");
    Legend->AddEntry(Histogram, Form("MC, #mu=%.1f#pm%.1f", Histogram->GetMean(), Histogram->GetRMS()), "l");
    Legend->AddEntry(Fitter, Form("Gauss, #mu=%.1f#pm%.1f", Fitter->GetParameter(1), Fitter->GetParameter(2)), "l");
    std::string TempTitle = std::string(Histogram->GetName());
    TempTitle += "_canv";
    TCanvas *TempCanvas = new TCanvas(TempTitle.c_str(), TempTitle.c_str(), 1024, 1024);
    TempCanvas->SetGridx();
    TempCanvas->SetGridy();
    TempCanvas->SetRightMargin(0.03);
    TempCanvas->SetBottomMargin(0.08);
    TempCanvas->SetLeftMargin(0.10);
    TempCanvas->SetTopMargin(0.06);
    TempCanvas->cd();
    Histogram->Draw();
    TempLine->Draw("same");
    Fitter->Draw("same");
    Legend->Draw("same");
    TempCanvas->Write();
    Histogram->Write();
        
    delete TempLine;
    delete TempCanvas;
    delete Fitter;
    delete Legend;
}

// ****************
// Make a ratio histogram
template<class HistType>
HistType* SampleSummary::RatioHists(HistType *NumHist, HistType *DenomHist) {
// ****************

  HistType *NumCopy = (HistType*)(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());
  NumCopy->Divide(DenomHist);

  return NumCopy;
}

// ****************
// Make a ratio th2poly
TH2Poly* SampleSummary::RatioPolys(TH2Poly *NumHist, TH2Poly *DenomHist) {
// ****************

  TH2Poly *NumCopy = (TH2Poly*)(NumHist->Clone());
  std::string title = std::string(DenomHist->GetName()) + "_ratio";
  NumCopy->SetNameTitle(title.c_str(), title.c_str());

  for(int i=1; i < NumCopy->GetNumberOfBins()+1; i++)
  {
    NumCopy->SetBinContent(i,NumHist->GetBinContent(i)/DenomHist->GetBinContent(i));
  }

  return NumCopy;
}

// ****************
// Normalise a histogram
template<class HistType>
HistType* SampleSummary::NormaliseHist(HistType *Histogram) {
// ****************

  HistType* HistCopy = (HistType*)(Histogram->Clone());
  HistCopy->Scale(1./HistCopy->Integral("width"), "width");
  std::string title = std::string(HistCopy->GetName())+"_norm";
  HistCopy->SetNameTitle(title.c_str(), title.c_str());

  return HistCopy;
}

// ****************
// Calculate the LLH for TH1D, set the LLH to title of MCHist
void SampleSummary::CalcLLH(TH1D * const & DataHist, TH1D * const & MCHist, TH1D * const & W2Hist) {
// ****************
  double llh = GetLLH(DataHist, MCHist, W2Hist);
  std::stringstream ss;
  ss << "_2LLH=" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  std::cout << std::setw(55) << std::left << MCHist->GetName() << std::setw(10) << DataHist->Integral() << std::setw(10) << MCHist->Integral() << std::setw(10) << llh << std::endl;
}

// ****************
// Calculate the LLH for TH1D, set the LLH to title of MCHist
void SampleSummary::CalcLLH(TH2Poly * const & DataHist, TH2Poly * const & MCHist, TH2Poly * const & W2Hist) {
// ****************
  double llh = GetLLH(DataHist, MCHist, W2Hist);
  std::stringstream ss;
  ss << "_2LLH=" << llh;
  MCHist->SetTitle((std::string(MCHist->GetTitle())+ss.str()).c_str());
  std::cout << std::setw(55) << std::left << MCHist->GetName() << std::setw(10) << NoOverflowIntegral(DataHist) << std::setw(10) << NoOverflowIntegral(MCHist) << std::setw(10) << llh << std::endl;
}

// ****************
double SampleSummary::GetLLH(TH2Poly * const & DataHist, TH2Poly * const & MCHist, TH2Poly * const & W2Hist) {
  // ****************
  double llh = 0.0;
  for (int i = 1; i < DataHist->GetNumberOfBins()+1; ++i) 
  {
    double data = DataHist->GetBinContent(i);
    double mc = MCHist->GetBinContent(i);
    double w2 = W2Hist->GetBinContent(i);
    llh += sampleND->getTestStatLLH(data, mc, w2);
  }
  //KS: do times 2 because banff reports chi2
  return 2*llh;
}

// ****************
double SampleSummary::GetLLH(TH1D * const & DataHist, TH1D * const & MCHist, TH1D * const & W2Hist) {
// ****************
  double llh = 0.0;
  for (int i = 1; i <= DataHist->GetXaxis()->GetNbins(); ++i) 
  {
    double data = DataHist->GetBinContent(i);
    double mc = MCHist->GetBinContent(i);
    double w2 = W2Hist->GetBinContent(i);
    llh += sampleND->getTestStatLLH(data, mc, w2);
  }
  //KS: do times 2 because banff reports chi2
  return 2*llh;
}

// ****************
// Make a projection
TH1D* SampleSummary::ProjectHist(TH2D* Histogram, bool ProjectX) {
// ****************

  TH1D* Projection = NULL;
  std::string name;
  if (ProjectX) {
    name = std::string(Histogram->GetName()) + "_x";
    Projection = Histogram->ProjectionX(name.c_str(), 1, Histogram->GetYaxis()->GetNbins(), "e");
  } else {
    name = std::string(Histogram->GetName()) + "_y";
    Projection = Histogram->ProjectionY(name.c_str(), 1, Histogram->GetXaxis()->GetNbins(), "e");
  }

  return Projection;
}

// ****************
// Make a projection
TH1D* SampleSummary::ProjectPoly(TH2Poly* Histogram, bool ProjectX, int selection, bool MakeErrorHist) {
// ****************

  std::vector<double> xbins;
  std::vector<double> ybins;

  sampleND->SetupBinning(selection, xbins, ybins);
  TH1D* Projection = NULL;
  std::string name;
  if (ProjectX) {
    name = std::string(Histogram->GetName()) + "_x";
    Projection = PolyProjectionX(Histogram,name.c_str(),xbins, MakeErrorHist);
  } 
  else {
    name = std::string(Histogram->GetName()) + "_y";
    Projection = PolyProjectionY(Histogram,name.c_str(),ybins, MakeErrorHist);
  }

  return Projection;
}


// ****************
//KS: We have two methdos how to apply statystical flcutuation standard is faster hence is default
TH1D* SampleSummary::MakeFluctuatedHistogram(TH1D* PolyHist){
// ****************
    TH1D *FluctHist = NULL;
    if(StandardFluctuation) FluctHist = MakeFluctuatedHistogramStandard(PolyHist);
    else FluctHist = MakeFluctuatedHistogramAlternative(PolyHist);
    
    return FluctHist;
}

// ****************
//KS: We have two methdos how to apply statystical flcutuation standard is faster hence is default
TH2Poly* SampleSummary::MakeFluctuatedHistogram(TH2Poly* PolyHist){
// ****************
    //KS: For TH2Poly alternative method is 2100 times slower!!! However resuts are similiar hence we normally don't use it
    TH2Poly *FluctHist = NULL;
    if(StandardFluctuation) FluctHist = MakeFluctuatedHistogramStandard(PolyHist);
    else FluctHist = MakeFluctuatedHistogramAlternative(PolyHist);
       
    return FluctHist;
}

// ****************
// Make Poisson Fluctuation of TH1D hist
TH1D* SampleSummary::MakeFluctuatedHistogramStandard(TH1D* PolyHist){
// ****************
      TH1D *FluctHist = (TH1D*)(PolyHist->Clone());
      // Make the Poisson fluctuated hist
      FluctHist->Reset("");

      for (int i = 1; i <= PolyHist->GetXaxis()->GetNbins(); ++i)  
      {
          // Get the posterior predictive bin content
          double MeanContent = PolyHist->GetBinContent(i);
          // Get a Poisson fluctuation of the content
          double Random = rnd->PoissonD(MeanContent);
          // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
          FluctHist->SetBinContent(i,Random);
    }
    return FluctHist;
}

// ****************
// Make Poisson Fluctuation of TH2Poly hist
TH2Poly* SampleSummary::MakeFluctuatedHistogramStandard(TH2Poly* PolyHist){
// ****************
      TH2Poly *FluctHist = (TH2Poly*)(PolyHist->Clone());
      // Make the Poisson fluctuated hist
      FluctHist->Reset("");

      for (int i = 1; i < FluctHist->GetNumberOfBins()+1; ++i) 
      {
          // Get the posterior predictive bin content
          double MeanContent = PolyHist->GetBinContent(i);
          // Get a Poisson fluctuation of the content
          double Random = rnd->PoissonD(MeanContent);
          // Set the fluctuated histogram content to the Poisson variation of the posterior predictive histogram
          FluctHist->SetBinContent(i,Random);
    }
    return FluctHist;
}

// ****************
// Make Poisson Fluctuation of TH1D hist
TH1D* SampleSummary::MakeFluctuatedHistogramAlternative(TH1D* PolyHist){
// ****************
    TH1D *FluctHist = (TH1D*)(PolyHist->Clone());
    // Make the Poisson fluctuated hist
    FluctHist->Reset("");

    double evrate = PolyHist->Integral();
    int num = rnd->PoissonD(evrate);
    int count = 0;
    while(count < num)
    {
        double candidate = PolyHist->GetRandom();
        FluctHist->Fill(candidate);                 
        count++;
    }
    return FluctHist;
}

// ****************
// Make Poisson fluctuation of TH2Poly hist
TH2Poly* SampleSummary::MakeFluctuatedHistogramAlternative(TH2Poly* PolyHist){
// ****************
    TH2Poly *FluctHist = (TH2Poly*)(PolyHist->Clone());
    // Make the Poisson fluctuated hist
    FluctHist->Reset("");
    FluctHist->Fill(0.0, 0.0, 0.0);

    double evrate = NoOverflowIntegral(PolyHist);
    int num = rnd->PoissonD(evrate);
    int count = 0;
    while(count < num)
    {
        int iBin = GetRandomPoly2(PolyHist);                
        FluctHist->SetBinContent(iBin, FluctHist->GetBinContent(iBin) + 1);                
        count++;
    }
    return FluctHist;
}

// ****************
//KS: ROOT developers were too lazy do develop getRanom2 for TH2Poly, this implementation is based on:
// https://root.cern.ch/doc/master/classTH2.html#a883f419e1f6899f9c4255b458d2afe2e
int SampleSummary::GetRandomPoly2(TH2Poly* PolyHist){
// ****************
   int nbins = PolyHist->GetNumberOfBins();
   double r1 = rnd->Rndm();
   
    double fIntegral[nbins+2];
    fIntegral[0] = 0.0;
 
    //KS: This is custom version of ComputeIntegral, once again ROOT was lazy :(
    for (int i = 1; i < nbins+1; ++i) 
    {
        fIntegral[i] = 0.0;
        double content = PolyHist->GetBinContent(i);
        fIntegral[i] += fIntegral[i - 1] + content;
    }
    for (Int_t bin = 1; bin < nbins+1; ++bin)  fIntegral[bin] /= fIntegral[nbins];
    fIntegral[nbins+1] = PolyHist->GetEntries();
    
   //KS: We just return one rather then X and Y, this way we can use SetBinContent rather than Fill, which is faster 
   int iBin = TMath::BinarySearch(nbins, fIntegral,(Double_t) r1);
   //KS: Have to increment because TH2Poly has stupid offset arghh
   iBin += 1;
   return iBin;
}

// ****************
// Get the mode error from a TH1D
double SampleSummary::GetModeError(TH1D * const hpost){
// ****************  
  // Get the bin which has the largest posterior density
  int MaxBin = hpost->GetMaximumBin();

  // The total integral of the posterior
  double integral = hpost->Integral();

  double sum = 0.0;
  int LowBin = MaxBin-1;
  int HighBin = MaxBin+1;
  double LowCon = 0.0;
  double HighCon = 0.0;
  while (sum/integral < 0.6827 && (LowBin >= 0 && HighBin < hpost->GetNbinsX()+1) ) {

    // Get the slice
    //KS:: If each slice reached histogram end then set value to 0, then other slice will be able to move further
    if(LowBin >= 0)LowCon = hpost->GetBinContent(LowBin);
    else LowCon = 0.0;
        
    if(HighBin < hpost->GetNbinsX()+1){HighCon = hpost->GetBinContent(HighBin);}
    else HighCon = 0.0;

    // If we're on the last slice and the lower contour is larger than the upper
    if ((sum+LowCon+HighCon)/integral > 0.6827 && LowCon > HighCon) {
      sum += LowCon;
      break;
      // If we're on the last slice and the upper contour is larger than the lower
    } else if ((sum+LowCon+HighCon)/integral > 0.6827 && HighCon >= LowCon) {
      sum += HighCon;
      break;
    } else {
      sum += LowCon + HighCon;
    }

    //KS:: Move further only if you haven't reached histogram end
    if(LowBin >= 0) LowBin--;
    if(HighBin < hpost->GetNbinsX()+1) HighBin++;
  }

  double Mode_Error = 0.0;
  if (LowCon > HighCon) {
    Mode_Error = fabs(hpost->GetBinCenter(LowBin)-hpost->GetBinCenter(MaxBin));
  } else {
    Mode_Error = fabs(hpost->GetBinCenter(HighBin)-hpost->GetBinCenter(MaxBin));
  }

return Mode_Error;
}
