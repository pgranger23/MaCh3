#include "samplePDFND.h"

// Debug this by defining DEBUG_DUMP if you find differences between CPU and GPU weights
// Warning: will write lots of ROOT files to your directory!
// Also turns on when DEBUG_ND280_DUMP environment variables is defined
//#define DEBUG_DUMP
//#define DEBUG

// ***************************************************************************
// Alternative constructor to make the sample aware of the manager
samplePDFND::samplePDFND(manager *Manager) : samplePDFBase(Manager->GetPOT()) {
// ***************************************************************************
 
  std::cout << "Creating ND instance" << std::endl;
#if USE_SPLINE == USE_TSpline3
  std::cout << "   Using TSpline3" << std::endl;
#elif USE_SPLINE == USE_TSpline3_red   
  std::cout << "   Using TSpline3 reduced" << std::endl;
#elif USE_SPLINE == USE_TF1
  std::cout << "   Using TF1" << std::endl;  
#elif USE_SPLINE == USE_TF1_red
  std::cout << "   Using TF1 reduced" << std::endl;
#elif USE_SPLINE == USE_Akima_Spline
  std::cout << "   Using Akima Spline" << std::endl;
#elif USE_SPLINE == USE_Truncated_TSpline
  std::cout << "   Using Truncated TSpline3" << std::endl;
#elif USE_SPLINE == USE_Truncated_Akima
  std::cout << "   Using Truncated Akima Spline" << std::endl;
#elif USE_SPLINE == USE_Monotone_Spline
  std::cout << "   Using Monotone Spline" << std::endl;
#endif
  // Save a pointer to the manager
  FitManager = Manager;

  // Random start
  HaveIRandomStart = false;

  nSamples = 0;

  samplePDF_data_array = NULL;
  samplePDF_array = NULL;
  samplePDF_w2_array = NULL;
#ifdef MULTITHREAD
  samplePDF_mode_array = NULL;
#endif

  // Pointer to covarianceNDDetPoly
  NDDetCov = NULL;
  // Pointer to covarianceXsec
  XsecCov    = NULL;
  
  // A nice random number to use at our leisure (seeded on TIME)
  rnd = new TRandom3(0);
  // Do we care about the MC pdfs by mode
  samplepdfs  = new TObjArray(0);

  // Holds the data pdfs
  datapdfs    = new TObjArray(0);
  
  // Include making distributions by mode?
  modepdf = false;
  // Holds the MC pdfs by mode
  samplemodepdfs = NULL;
  modeobjarray = NULL;

  nModes = 0;
  nModes = GetNModes();

  ndims = NULL;
  kinvars = NULL;
   
  // Variables related to MC stat
  firsttime = true;
  UpdateW2 = true;
  UpdateW2 = FitManager->GetUpdateW2();
  
  // Set the ND test-statistic
  SetTestStatistic(static_cast<TestStatistic>(FitManager->GetMCStatLLH()));

  // Get which splines should use a linear function
  linearsplines = FitManager->GetXsecLinearSpline();

  // Then load the samples
  LoadSamples();
  
  // Enable the mode histograms AFTER addSelection is called
  if (FitManager->getPlotByMode()) EnableModeHistograms();

  InitExperimentSpecific();

#ifdef CUDA
std::cout << "- Using ND GPU version " << std::endl;
  #ifdef DEBUG_DUMP
  badWeight = 0;
  nReconf = 0;
  #endif
#endif
// Check the preprocessors that have been defined
// Can currently only run in verbose DEBUG mode in single thread
#if DEBUG > 0
    #ifdef MULTITHREAD
    #pragma omp parallel
   {
    #pragma omp critical
     if (omp_get_num_threads() != 1) {
       std::cerr << "I'm running in DEBUG on multiple threads, this won't work" << std::endl;
       std::cerr << "export OMP_NUM_THREADS=1 and try again please!" << std::endl;
       throw;
     }
    #endif
   }
#endif
}

// ***************************************************************************
// Destructor
samplePDFND::~samplePDFND() {
// ***************************************************************************
  std::cout << "Destroying samplePDFND object" << std::endl;
  for(__int__ i = 0; i < nSamples; i++)
  {
    if (samplepdfs->At(i) == NULL) continue;
    for(__int__ j = 0; j < maxBins[i]; j++)
    {
        delete polybins[i][j];
    }
    delete[] polybins[i];
    delete W2Hist[i];
    delete[] kinvars[i];;
  }

   // Delete the master array
   for (__int__ i = 0; i < nSamples; i++)
   {
     delete[] samplePDF_data_array[i];
     delete[] samplePDF_array[i];
     delete[] samplePDF_w2_array[i];
   }
   delete[] samplePDF_data_array;
   delete[] samplePDF_array;
   delete[] samplePDF_w2_array;
#ifdef MULTITHREAD
   // Delete the master mode array
   if (modepdf) {
     for (__int__ i = 0; i < nSamples; i++) {
       for (__int__ j = 0; j < nModes+1; j++) {
         delete[] samplePDF_mode_array[i][j];
       }
       delete[] samplePDF_mode_array[i];
     }
     delete[] samplePDF_mode_array;
   }
#endif
   
  delete rnd;
  delete samplepdfs;
  delete datapdfs;

  if(modepdf)
  {
      delete samplemodepdfs;
      delete[] modeobjarray;
  }

  delete[] xsecInfo;
  delete[] maxBins;
  delete[] ndims;
  delete[] kinvars;
  delete[] polybins;
  delete FitManager;

  if (XsecCov)  delete XsecCov;
  if (NDDetCov) delete NDDetCov;
#ifdef CUDA
  delete splineMonolith;
#endif
}

// ***************************************************************************
//KS: Each experiment uses some specyfic variables, make template function for it
void samplePDFND::InitExperimentSpecific() {
// ***************************************************************************

  //WARNING FIXME TODO This is for the time being
  /*
  // What production we're using
  TString production = FitManager->getNDRuns();
  std::cout << "production: " << production << std::endl;
  // Currently not set in the manager. Dictates if we run with simple ND280 detector treatment, i.e. treat ND280 systematics as normalisations in pmu cosmu bins for each selection
  simple = true;


  Use2dEb = false;

  UseSandMC = true;
  setSandMC(FitManager->getUseSand());
  // Do some checks on the requested production
  if (production == "") {
    std::cerr << "You gave me an empty production!\nExiting!" << std::endl;
    throw;
  }
  if (!production.Contains("P6")) {
    std::cout << "Production: " << production.Data() << std::endl;
    std::cerr << "I only support P6 sorry!" << std::endl;
    throw;
  }

  TString prod_copy(production);
  prod_copy.Remove(prod_copy.First(' '), prod_copy.Length() - prod_copy.First(' '));
  if (prod_copy != "P6" && prod_copy != "P6_S17" && prod_copy != "P6T") {
    std::cerr << "Did not find good good production, exiting" << std::endl;
    std::cerr << "You gave " << prod_copy << ", I need P6 or P6_S17 or P6T" << std::endl;
    throw;
  }
  // Do we really need a new here? Surely setting the string should be fine
  // Warning, something does indeed break, investigate this in the future
  prod = new TString(production);

  if (std::getenv("NIWG_ROOT") == NULL) {
    std::cerr << "Need NIWG environment variable to run with Eb" << std::endl;
    std::cerr << "EXPORT NIWG to your NIWGReWeight environment" << std::endl;
    throw;
  }
  */
}

// ***************************************************************************
// Load up the samples and binning specified by the user
void samplePDFND::LoadSamples() {
// ***************************************************************************
  std::cerr<<"Function LoadSamples is experiment specific however core code uses it"<<std::endl;
  std::cerr<<"Since you haven't implemented it I have to stop it"<<std::endl;
  throw;
} // end LoadSamples

// ***************************************************************************
// Set the datapdfs to the MC-pdfs and treat them as Asimov data
// This might require throwing the starting parameters of the Markov Chain to avoid a stationary chain at the first n steps (n might be as large as 10k!)
void samplePDFND::setAsimovFakeData(bool CustomReWeight) {
// ***************************************************************************

   std::cout << "Setting ND sample to use Asimov as data..." << std::endl;
   if (HaveIRandomStart) {
     std::cerr << "samplePDFND has been set to random start and been reweighted to it" << std::endl;
     std::cerr << "And now being asked to set Asimov data to the random start... I think this is wrong!" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   // Check if the relevant covariances have been set
   CheckCovariances();

   // Set Asimov parameters
   // Could use config here instead but will hard-code for now
   if (CustomReWeight) {
     XsecCov->setParProp(6, 1.0);
     XsecCov->setParProp(7, 1.0);
     double *fake = NULL;
     reweight(fake);
     // Reset the fParCurr and nominal parameters to nominal after the reweight
     std::cout << "Set samplePDFND asimov with cross-section parameters: " << std::endl;
     XsecCov->printNominalCurrProp();
     std::cout << "Now resetting..." << std::endl;
     XsecCov->setParameters();
   }

   // Loop over all the samples
   for (int i = 0; i < nSamples; ++i) {

     // Move to next sample if we haven't enabled this one
     if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

     // Can't see why we wouldn't run 2D
     if (ndims[i] != 2) {
       std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleName[i] << std::endl;
       std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
       throw;
     }

     // Strip out the whitespaces and replace them with underscores in the sample names
     // Will be used in Clone
     std::string temp = SampleName[i];
     while (temp.find(" ") != std::string::npos) {
       temp.replace(temp.find(" "), 1, std::string("_"));
     }

     // Just clone the MC for the simple setAsimovFakeData()
     TH2Poly* MCpdf = (TH2Poly*)(getPDF(i))->Clone(temp.c_str());
     // Pass the pointer to addData
     addData(MCpdf, i);
   } // end for loop

   //KS: Update data pdf from new inserted histograms
   UpdateDataPDF();

   std::cout << "------------------------------------------" << std::endl;
   std::cout << "Have set Asimov fake-data to reweighted MC" << std::endl;
   std::cout << "With cross-section settings: " << std::endl;
   XsecCov->printNominalCurrProp();
   //NDDetCov->printNominalCurrProp();
   std::cout << "------------------------------------------" << std::endl;

   // Print the post-Asimov rates
   printRates();
} // end setAsimovFakeData()

// ***************************************************************************
// Setup the binning for a given sample
// mosly used by SigmaVar, could be called in addSelection...
void samplePDFND::SetupBinning(int Selection, std::vector<double> & BinningX, std::vector<double> & BinningY) {
// ***************************************************************************
  std::cerr<<"Function SetupBinning is experiment specific however core code uses it"<<std::endl;
  std::cerr<<"Since you haven't implemented it I have to stop it"<<std::endl;
  throw;
} // end SetupBinning

// ***************************************************************************
// Set the datapdfs to Poisson fluctuated MC-pdfs and treat them as Asimov data
// This might require throwing the starting parameters of the Markov Chain to avoid a stationary chain at the first n steps (n might be as large as 10k!)
void samplePDFND::setAsimovFakeDataFluctuated(bool CustomReWeight) {
// ***************************************************************************

   std::cout << "Setting ND sample to use statistically fluctuated Asimov as data..." << std::endl;
   if (HaveIRandomStart) {
     std::cerr << "samplePDFND has been set to random start and been reweighted to it" << std::endl;
     std::cerr << "And now being asked to set Asimov data to the random start... I think this is wrong!" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   // Check if the relevant covariances have been set
   CheckCovariances();

   // Set Asimov parameters
   // Could use config here instead but will hard-code for now
   if (CustomReWeight) {
     XsecCov->setParProp(6, 1.0);
     XsecCov->setParProp(7, 1.0);
     double *fake = NULL;
     reweight(fake);
     // Reset the fParCurr and nominal parameters to nominal after the reweight
     std::cout << "Set samplePDFND asimov with cross-section parameters: " << std::endl;
     XsecCov->printNominalCurrProp();
     std::cout << "Now resetting..." << std::endl;
     XsecCov->setParameters();
   }

   // Loop over all the samples
   for (int i = 0; i < nSamples; ++i) {

     // Move to next sample if we haven't enabled this one
     if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

     // Can't see why we wouldn't run 2D
     if (ndims[i] != 2) {
       std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleName[i] << std::endl;
       std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
       throw;
     }

     // Strip out the whitespaces and replace them with underscores in the sample names
     // Will be used in Clone
     std::string temp = SampleName[i];
     while (temp.find(" ") != std::string::npos) {
       temp.replace(temp.find(" "), 1, std::string("_"));
     }

     // Clone the MC, then Poisson fluctuate each bin content
     TH2Poly* MCpdf = (TH2Poly*)(getPDF(i))->Clone(temp.c_str());
     for (int ibin = 1; ibin < MCpdf->GetNumberOfBins()+1; ibin++){
       MCpdf->SetBinContent(ibin, rnd->PoissonD((MCpdf->GetBinContent(ibin))));
     }
    
     // Pass the pointer to addData
     addData(MCpdf, i);
   } // end for loop

   //KS: Update data pdf from new inserted histograms
   if(samplePDF_data_array != NULL) UpdateDataPDF();

   std::cout << "------------------------------------------" << std::endl;
   std::cout << "Have set Asimov fake-data to fluctuated MC" << std::endl;
   std::cout << "With cross-section settings: " << std::endl;
   XsecCov->printNominalCurrProp();
   //NDDetCov->printNominalCurrProp();
   std::cout << "------------------------------------------" << std::endl;

   // Print the post-Asimov rates
   printRates();
} // end setAsimovFakeData()

// ***************************************************************************
void samplePDFND::setAsimovFakeData_FromFile(std::string &FileName) {
// ***************************************************************************
   std::cout << "Setting ND sample to use fakedata from file " << FileName << " as data..." << std::endl;

   if (HaveIRandomStart) {
     std::cerr << "samplePDFND has been set to random start and been reweighted to it" << std::endl;
     std::cerr << "And now being asked to set Asimov data to the random start... I think this is wrong!" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   // Check if the relevant covariances have been set
   CheckCovariances();

   // Open the TFile
   TFile *AsimovFile = new TFile(FileName.c_str(), "OPEN");
   if (!AsimovFile->IsOpen()) {
     std::cerr << "************" << std::endl;
     std::cerr << "Provided Asimov file " << FileName << " does not exist, exiting" << std::endl;
     std::cerr << "************" << std::endl;
     throw;
   }
   AsimovFile->ls();

   // Loop over all the samples
   for (int i = 0; i < nSamples; ++i) {

     // Move to next sample if we haven't enabled this one
     if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

     // Can't see why we wouldn't run 2D
     if (ndims[i] != 2) {
       std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleName[i] << std::endl;
       std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
       throw;
     }

     // Strip out the whitespaces and replace them with underscores in the sample names
     // Will be used in Clone
     std::string temp = SampleName[i];
     while (temp.find(" ") != std::string::npos) {
       temp.replace(temp.find(" "), 1, std::string("_"));
     }

     // Clone the MC from the TFile
     TH2Poly* MCpdf = (TH2Poly*)(AsimovFile->Get((std::string("MC_")+temp).c_str())->Clone());

     if (MCpdf == NULL) {
       std::cerr << "Did not file input Asimov in " << AsimovFile << " for " << SampleName[i] << std::endl;
       throw;
     }
     MCpdf->SetDirectory(0);

     // Pass the pointer to addData
     addData(MCpdf, i);
   } // end for loop

   //KS: Update data pdf from new inserted histograms
   if(samplePDF_data_array != NULL) UpdateDataPDF();

   std::cout << "------------------------------------------" << std::endl;
   std::cout << "Have set Asimov fake-data to reweighted MC" << std::endl;
   std::cout << "With cross-section settings: " << std::endl;
   XsecCov->printNominalCurrProp();
   //NDDetCov->printNominalCurrProp();
   std::cout << "------------------------------------------" << std::endl;

   // Print the post-Asimov rates
   printRates();

   delete AsimovFile;
}

// ***************************************************************************
// This function is used for setting data from a file; it is slightly different
// from the function used to set fake data for asimovs, as it doesn't check the
// MC rates at the same time.
void samplePDFND::setDataFromFile(std::string &FileName) {
// ***************************************************************************

  if (HaveIRandomStart) {
    std::cerr << "samplePDFND has been set to random start and been reweighted to it" << std::endl;
    std::cerr << "And now being asked to set Asimov data to the random start... I think this is wrong!" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  // Check if the relevant covariances have been set
  CheckCovariances();

  // Open the TFile
  TFile *file = new TFile(FileName.c_str(), "OPEN");
  if (!file->IsOpen()) {
    std::cerr << "************" << std::endl;
    std::cerr << "Provided data file " << FileName << " does not exist, exiting" << std::endl;
    std::cerr << "************" << std::endl;
    throw;
  }
  file->ls();

  // Loop over all the samples
  for (int i = 0; i < nSamples; ++i) {

    // Move to next sample if we haven't enabled this one
    if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

    // Can't see why we wouldn't run 2D
    if (ndims[i] != 2) {
      std::cerr << "Can not set data for other than 2D histograms for " << SampleName[i] << std::endl;
      std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }

    // Strip out the whitespaces and replace them with underscores in the sample names
    // Will be used in Clone
    std::string temp = SampleName[i];
    while (temp.find(" ") != std::string::npos) {
      temp.replace(temp.find(" "), 1, std::string("_"));
    }

    // Clone the MC from the TFile
    TH2Poly* MCpdf = (TH2Poly*)(file->Get((std::string("DATA_")+temp).c_str())->Clone());

    if (MCpdf == NULL) {
      std::cerr << "Did not find input data in " << file << " for " << SampleName[i] << std::endl;
      throw;
    }
    MCpdf->SetDirectory(0);

    // Pass the pointer to addData
    addData(MCpdf, i);
  } // end for loop

   //KS: Update data pdf from new inserted histograms
   if(samplePDF_data_array != NULL) UpdateDataPDF();

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Have set data from file. Here are the rates:" << std::endl;

  // Print the post-Asimov rates
  printRates(true);

  delete file;
}


// ***************************************************************************
// Similar to samplePDFND::setAsimovFakeData()
// Throws the parameters that the fake-data is generated at first so that fake-data is not the nominal parameter set
// For this to work we need the associated xsec and ND covariances
void samplePDFND::setAsimovFakeDataThrow() {
// ***************************************************************************

   std::cout << "Setting ND sample to use systematically fluctuated Asimov as data..." << std::endl;

   // Check if the relevant covariances have been set
   CheckCovariances();

   // Don't want nominal values so set bool to false
   // Throw xsec
   XsecCov->throwNominal(false);
   XsecCov->setParameters();
   XsecCov->printNominalCurrProp();
   // Throw detector
   NDDetCov->throwNominal(false);
   NDDetCov->setParameters();
   // Don't necessarily wan't to print this bad-boy (580 parameters)
   //NDDetCov->printNominalCurrProp();

   // Now reweight the samplePDF with the new thrown parameters
   double *fake = 0;
   reweight(fake);

   // Loop over all the samples
   for (int i = 0; i < nSamples; ++i) {

     // Move to next sample if we haven't enabled this one
     if (samplepdfs->At(i) == NULL || datapdfs->At(i) == NULL) continue;

     if (ndims[i] != 2) {
       std::cerr << "Can not set Asimov data for other than 2D histograms for " << SampleName[i] << std::endl;
       std::cerr << "Simply a matter of implementation; to edit this see " << __FILE__ << ":" << __LINE__ << std::endl;
       throw;
     }

     // Strip out the whitespaces and replace them with underscores in the sample names
     // Will be used in Clone
     std::string temp = SampleName[i];
     while (temp.find(" ") != std::string::npos) {
       temp.replace(temp.find(" "), 1, std::string("_"));
     }

     // Just clone the MC for the simple setAsimovFakeData()
     TH2Poly* MCpdf = (TH2Poly*)(getPDF(i))->Clone(temp.c_str());
     // Pass the pointer to addData
     addData(MCpdf, i);
   } 
   // end for loop

   //KS: Update data pdf from new inserted histograms
   if(samplePDF_data_array != NULL) UpdateDataPDF();

   // Let's throw the parameters so that MC != Asimov data
   // Don't want nominal values so set bool to false
   // Throw xsec
   XsecCov->throwNominal(false);
   XsecCov->setParameters();
   XsecCov->printNominalCurrProp();
   // Throw detector
   NDDetCov->throwNominal(false);
   NDDetCov->setParameters();
   // Don't necessarily wan't to print this bad-boy
   //NDDetCov->printNominalCurrProp();
   reweight(fake);

   std::cout << "------------------------------------------" << std::endl;
   std::cout << "Have set Asimov fake-data to varied reweighted MC" << std::endl;
   std::cout << "------------------------------------------" << std::endl;

   // Print the post-Asimov rates
   printRates();
} // end setAsimovFakeData()

// ***************************************************************************
// Function which enables the samplePDFs to be plotted in accordance to interaction mode
void samplePDFND::EnableModeHistograms() {
// ***************************************************************************

   if (modepdf == true) {
     std::cout << "Already enabled modepdf but you're trying to do so again" << std::endl;
     std::cout << "Just returning to caller..." << std::endl;
     return;
   }

   modepdf = true;

   // Expand to make all samples
   samplemodepdfs = new TObjArray(0);
   samplemodepdfs->Expand(nSamples);
   samplemodepdfs->SetOwner(true);
   modeobjarray = new TObjArray*[nSamples]();
   bool valid = false;

   for (int i = 0; i < nSamples; i++)
   {
     if (samplepdfs->At(i) == NULL) continue;
     // Make sure at least one sample is added (unfortunately the calling of EnableModeHistograms() needs to happen _AFTER_ addSelection, or else we have no added selections at all
     valid = true;

     modeobjarray[i] = new TObjArray(0);
     modeobjarray[i]->Expand(nModes+1);
     modeobjarray[i]->SetOwner(true);

     for (int j = 0; j < nModes+1; j++) {
       TString name = SampleName[i].c_str();
       name += j;
       if (ndims[i] == 1) {
         TH1D* tmp = ((TH1D*)samplepdfs->At(i));
         TH1D* clone = (TH1D*)tmp->Clone(name);
         modeobjarray[i]->AddAt(clone, j);
       } else if (ndims[i] == 2) {
         TH2Poly* tmp = ((TH2Poly*)samplepdfs->At(i));
         TH2Poly* clone = (TH2Poly*)tmp->Clone(name);
         modeobjarray[i]->AddAt(clone, j);
         // Complain if ndims is some weird thing
       } else {
         std::cerr << "ndims[" << i << "] != 1 or 2 in EnableModeHistograms()" << std::endl;
         std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
         throw;
       }
     } // End loop over interaction modes
     samplemodepdfs->AddAt(modeobjarray[i], i);
   } // End loop over samples

   if (valid == false) {
     std::cerr << "Tried enabling mode pdfs without having any selections enabled" << std::endl;
     std::cerr << "Did you call EnableModeHistograms() before you did addSelection()?" << std::endl;
     std::cerr << "Call addSelection() first, then EnableModeHistograms()" << std::endl;
     throw;
   }
}


// ***************************************************************************
// Helper function to print rates for the samples with LLH
void samplePDFND::printRates(bool dataonly) {
// ***************************************************************************
    std::cout << std::setw(40) << std::left << "Sample" << std::setw(10) << "Data" << std::setw(10);
    if(!dataonly)
       std::cout<< "MC" << std::setw(10) << "-LLH" << std::setw(10) << "|" << std::endl;
    else
       std::cout << std::endl;

   double sumData  = 0.0;
   double sumMC    = 0.0;
   double likelihood = 0.0;
   for (int i = 0; i < nSamples; i++) {
     std::string name = SampleName[i].c_str();

     if (datapdfs->At(i) != NULL) {
       sumData += NoOverflowIntegral(((TH2Poly*)datapdfs->At(i)));
       if(!dataonly)
       {
        sumMC   += NoOverflowIntegral((TH2Poly*)(getPDF(i)));
        /*for (int k = 0; k < nModes+1; k++) {
          std::string ModeName = std::string("MC")+name+"_"+Mode_ToString(k);
          std::cout << std::setw(40) << std::left << ModeName <<  NoOverflowIntegral(((TH2Poly*)getPDFMode(i,k))) << std::setw(10) << "|"<<std::endl;
        }*/
        likelihood = getSampleLikelihood(i);
       }
 
       std::cout << std::setw(40) << std::left << SampleName[i] << std::setw(10) << NoOverflowIntegral(((TH2Poly*)datapdfs->At(i)));
       if(!dataonly) std::cout << std::setw(10) << NoOverflowIntegral((TH2Poly*)(getPDF(i))) << std::setw(10) << likelihood << std::setw(10) << "|" << std::endl;
       else std::cout << std::endl;
     }
     //KS: Check make sure data and MC binning is the same
     CheckBinningMatch();
   }
   if(!dataonly)
      likelihood = getLikelihood();

   std::cout << std::setw(40) << std::left << "Total" << std::setw(10) << sumData << std::setw(10) ;
   if(!dataonly)
      std::cout<< sumMC << std::setw(10) << likelihood << std::setw(10) << "|" << std::endl;
   else
      std::cout << std::endl;
   //NDDetCov->printNominalCurrProp();
}


// ***************************************************************************
//KS: Helper function check if data and MC binning matches
void samplePDFND::CheckBinningMatch() {
// ***************************************************************************

  //WARNING remove it later it later
  /*
  //KS: Since in psyche independent we separetly give data PDF there is danger that data and MC will not match
  #ifndef PSYCHESETUP
  if(!dataonly && datapdfs->At(i) != NULL)
  {
    for(int j = 1; j < ((TH2Poly*)datapdfs->At(i))->GetNumberOfBins()+1; j++)
    {
      //KS: There is weird offset between bin content and GetBins so this is correct, inspite of looking funny
      TH2PolyBin* polybinMC = (TH2PolyBin*)(((TH2Poly*)samplepdfs->At(i))->GetBins()->At(j-1)->Clone());
      TH2PolyBin* polybinData = (TH2PolyBin*)(((TH2Poly*)datapdfs->At(i))->GetBins()->At(j-1)->Clone());

      if( std::fabs(polybinData->GetXMin() - polybinMC->GetXMin()) > 0.001 ||
          std::fabs(polybinData->GetXMax() - polybinMC->GetXMax()) > 0.001 ||
          std::fabs(polybinData->GetYMin() - polybinMC->GetYMin()) > 0.001 ||
          std::fabs(polybinData->GetYMax() - polybinMC->GetYMax()) > 0.001  )
      {
            std::cout<<"Sample "<<name<<" has different bin edges for data and MC "<<std::endl;
            std::cout<<"data  "<<" x min "<<polybinData->GetXMin()<<" x max "<<polybinData->GetXMax()<<" y min "<<polybinData->GetYMin()<<" y max "<<polybinData->GetYMax()<<std::endl;
            std::cout<<"mc    "<<" x min "<<polybinMC->GetXMin()  <<" x max "<<polybinMC->GetXMax()  <<" y min "<<polybinMC->GetYMin()  <<" y max "<<polybinMC->GetYMax()<<std::endl;
            std::cout<<" Most likely wrong psyche independent file was provided, contact local ND expert "<<std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
            throw;
      }
      delete polybinMC;
      delete polybinData;
    }
  }
#endif
*/
}
// ***************************************************************************
// Change the starting position of the chain by throwing the systematics
void samplePDFND::RandomStart() {
// ***************************************************************************

   HaveIRandomStart = true;

   CheckCovariances();
   
   std::vector<covarianceBase*> Systematics;
   if (XsecCov) Systematics.push_back(XsecCov);

   if (NDDetCov) Systematics.push_back(NDDetCov);

   // Reconfigure the systematics randomly
   for (std::vector<covarianceBase*>::iterator it = Systematics.begin(); it != Systematics.end(); it++) {
     (*it)->RandomConfiguration();
   }

   // Then reweight the simulation to the random configuration
   std::cout << "Reweighting the sample for a random start..." << std::endl;
   double *fake = NULL;
   reweight(fake);

   // Print the rates
   std::cout << "Printing rates after the random start:" << std::endl;
   printRates();
}

// ***************************************************************************
// Reweight the Monte-Carlo at ND for one step of the Markov Chain
void samplePDFND::reweight(double* oscpar) {
// ***************************************************************************

    // Reset the histograms before reweight
    ResetHistograms();
        
    // Prepare the weights
    // e.g. find the relevant spline segments, pass stuff to GPU
    PrepareWeights();

    // Call entirely different routine if we're running with openMP
    #ifdef MULTITHREAD
    ReWeight_MC_MP();
    #else 
    ReWeight_MC();
    #endif
   
    //KS: If you want to not update W2 wights then uncomment this line
    if(!UpdateW2) firsttime = false;
}

#ifndef MULTITHREAD
// ***************************************************************************
void samplePDFND::ReWeight_MC() {
// ***************************************************************************

   ReconfigureFuncPars();
   
   // Loop over all the events
   // Fill their weights for this reconfiguration
   for (unsigned int i = 0; i < nEvents; ++i) 
   {
      // If sample isn't enabled
      if (samplepdfs->At(NDEve[i].event_sample) == NULL) continue;

      double detw   = 1.;
      double xsecw  = 1.;

      // Get the cross-section weights
      xsecw = CalcXsecWeight(i);

      //KS:If weight is 0 or less skip event as it will have no effect on fit and we can decrease step time roughly by 8%
      if(xsecw <= 0) continue;

      // Get the detector weight
      //a bit hacky: if the event is in the overflow, we assign the detector bin as the last one in that row. BANFF are just not applying a weight though, but it would be messy to not assign a bin at all. So now we're still assigning that last bin in the row, but now applying the weight if we're above 30000. This only affects the overflow, and is a very small effect anyway. See covariance/covarianceNDDetPoly:getBin() for how the bin is assignment is done
      if(NDEve[i].mom < 30000 && NDEve[i].mom > 0) detw = NDDetCov->calcReWeight(NDEve[i].det_bin);

      CalcFuncPars(i);

      // Fill the MC
      double weight = NDEve[i].flux_w*detw*xsecw;
      int SampleId = NDEve[i].event_sample;
      int HistBin = NDEve[i].hist_bin;

      samplePDF_array[EventSample][NDEve[i].hist_bin] += weight;

      // This one is needed for Barlow-Beeston calculation
      // Could use TH2D::GetSumw2() instead but that wouldn't work for the multi-threaded case, so use this one for validation, and check against samplepdfs
      if (firsttime) samplePDF_w2_array[EventSample][NDEve[i].mode][maxBins[EventSample]+NDEve[i].hist_bin] += weight*weight;

      if (modepdf)
      {
        int Mode = NDEve[i].mode;
        ((TH2Poly*)((TObjArray*)(samplemodepdfs->At(SampleId)))->At(Mode))->SetBinContent(HistBin, weight+((TH2Poly*)((TObjArray*)(samplemodepdfs->At(SampleId)))->At(Mode))->GetBinContent(HistBin));
      }
   } // end for loop
}
#else
// ***************************************************************************
// Routine for ReWeight_MC_MP on openMP
void samplePDFND::ReWeight_MC_MP() {
// ***************************************************************************

   ReconfigureFuncPars();

   // The per-thread array
   double **samplePDF_array_private = NULL;
   // And the mode
   double ***samplePDF_mode_array_private = NULL;
   // The per-thread array of weight^2
   double **samplePDF_w2_array_private = NULL;

   // Declare the omp parallel region
   // The parallel region needs to stretch beyond the for loop!
  #pragma omp parallel private(samplePDF_array_private, samplePDF_mode_array_private, samplePDF_w2_array_private)
  {
     // private to each thread
     samplePDF_array_private = new double*[nSamples]();
     samplePDF_w2_array_private = new double*[nSamples]();
     for (__int__ i = 0; i < nSamples; ++i) {
       samplePDF_array_private[i] = new double[maxBins[i]]();
       samplePDF_w2_array_private[i] = new double[maxBins[i]]();
       for(__int__ j = 0; j < maxBins[i]; ++j) {
         samplePDF_array_private[i][j] = 0.0;
         samplePDF_w2_array_private[i][j] = 0.0;
       }
     }

     // Do the same for the mode arrays
     if (modepdf) {
       samplePDF_mode_array_private = new double**[nSamples]();
       for (__int__ i = 0; i < nSamples; ++i) {
         samplePDF_mode_array_private[i] = new double*[nModes+1]();
         for (__int__ j = 0; j < nModes+1; ++j) {
           samplePDF_mode_array_private[i][j] = new double[maxBins[i]]();
           for (__int__ k = 0; k < maxBins[i]; ++k) {
             samplePDF_mode_array_private[i][j][k] = 0.0;
           }
         }
       }
       // Set the samplePDF_mode_array_private to zero, even if we're aren't running modepdf
       // This is purely to avoid compiler warning
     } else {
       samplePDF_mode_array_private = new double**[1]();
       samplePDF_mode_array_private[0] = new double*[1]();
       samplePDF_mode_array_private[0][0] = new double[1]();
       samplePDF_mode_array_private[0][0][0] = 0.0;
     }

     // The master for loop: this is where we gain all performace by multi-threading
 #pragma omp for
     for (unsigned int i = 0; i < nEvents; ++i) 
     {
       int EventSample = NDEve[i].event_sample;
       if (samplepdfs->At(EventSample) == NULL) continue;

       double detw    = 1.;
       double xsecw   = 1.;
       
       //clock.Start();
       // Get the cross-section weights
       xsecw = CalcXsecWeight(i);
       
//KS: We can't skip in DEBUG mode otherwise ND COV only Event Rates will be wrong        
#ifndef DEBUG
       //KS:If weight is 0 or less skip event as it will have no effect on fit and we can decrease step time roughly by 8%
       if(xsecw <= 0) continue;
#endif      
       // Get the detector weight
       //a bit hacky: if the event is in the overflow, we assign the detector bin as the last one in that row. BANFF are just not applying a weight though, but it would be messy to not assign a bin at all. So now we're still assigning that last bin in the row, but now applying the weight if we're above 30000. This only affects the overflow, and is a very small effect anyway. See covariance/covarianceNDDetPoly:getBin() for how the bin is assignment is done
       if(NDEve[i].mom < 30000 && NDEve[i].mom > 0) detw = NDDetCov->calcReWeight(NDEve[i].det_bin);
       
       CalcFuncPars(i);
       
       // Is simply the total weight to apply to one event
       double WeightTot = NDEve[i].flux_w*detw*xsecw;
  
       // Fill each threads' private array up
       if(NDEve[i].hist_bin >= 0){
         samplePDF_array_private[EventSample][NDEve[i].hist_bin] += WeightTot;
         // And the w2
         samplePDF_w2_array_private[EventSample][NDEve[i].hist_bin] += WeightTot*WeightTot;
       } else {
         samplePDF_array_private[EventSample][maxBins[EventSample]+NDEve[i].hist_bin] += WeightTot;
         // And the w2
         samplePDF_w2_array_private[EventSample][maxBins[EventSample]+NDEve[i].hist_bin] += WeightTot*WeightTot;
       }
       // Fill each threads' private array up
       if (modepdf) 
       {
         if(NDEve[i].hist_bin >= 0){
           samplePDF_mode_array_private[EventSample][NDEve[i].mode][NDEve[i].hist_bin] += WeightTot;
         } else {
           samplePDF_mode_array_private[EventSample][NDEve[i].mode][maxBins[EventSample]+NDEve[i].hist_bin] += WeightTot;
         }
       }

 #if DEBUG > 0
       // Get limits of plot
       //This is a bit mad, SetBinContent does different things for the overflow than to the rest of the bins for a TH2Poly
       //The gift that keeps on giving...
       if(NDEve[i].hist_bin < 0)
       {
         XsecOnly[EventSample]->SetBinContent(NDEve[i].hist_bin, POTWeights[i]*xsecw);
         NDCovOnly[EventSample]->SetBinContent(NDEve[i].hist_bin, POTWeights[i]*detw);
         AllOnly[EventSample]->SetBinContent(NDEve[i].hist_bin, NDEve[i].flux_w*detw*xsecw);
       } else {
         XsecOnly[EventSample]->SetBinContent(NDEve[i].hist_bin, POTWeights[i]*xsecw+XsecOnly[EventSample]->GetBinContent(NDEve[i].hist_bin));
         NDCovOnly[EventSample]->SetBinContent(NDEve[i].hist_bin, POTWeights[i]*detw+NDCovOnly[EventSample]->GetBinContent(NDEve[i].hist_bin));
         AllOnly[EventSample]->SetBinContent(NDEve[i].hist_bin, NDEve[i].flux_w*detw*xsecw+AllOnly[EventSample]->GetBinContent(NDEve[i].hist_bin));
       }
       // Check for negative or large weights!
       bool foundweird = false;
       if (NDEve[i].flux_w < 0 || NDEve[i].flux_w > __LARGE_WEIGHT__) {
         std::cerr << "NEGATIVE WEIGHT flux weight = " << NDEve[i].flux_w << std::endl;
         std::cerr << "Event # " << i << std::endl;
         foundweird = true;
       }
       if (xsecw <= 0 || xsecw > __LARGE_WEIGHT__) {
         std::cerr << "NEGATIVE WEIGHT xsecw = " << xsecw << std::endl;
         std::cerr << "Event # " << i << std::endl;
         foundweird = true;
       }
       if (detw < 0 || detw > __LARGE_WEIGHT__) {
         std::cerr << "NEGATIVE WEIGHT det weight = " << detw << std::endl;
         std::cerr << "Event # " << i << std::endl;
         foundweird = true;
       }
       if (std::isnan(WeightTot)) {
         std::cerr << "Found a nan total weight!" << std::endl;
         foundweird = true;
       }
 #endif

       // Print the structs if we're interested in event by event breakdowns
 #if DEBUG == 2
       PrintStructs(xsecw, detw, i);
 #endif

       // If we found something weird print everything we've got!
 #if DEBUG > 0
       if (foundweird) DumpSplines(xsecw, detw, i);
 #endif
     } // end for loop and openMP for, but still keep the omp region!
     // Now we can write the individual arrays from each thread to the main array
       for (__int__ i = 0; i < nSamples; ++i) {
         for (__int__ j = 0; j < maxBins[i]; ++j)
         {
           #pragma omp atomic
           samplePDF_array[i][j] += samplePDF_array_private[i][j];
           if(firsttime)
           {
            #pragma omp atomic
            samplePDF_w2_array[i][j] += samplePDF_w2_array_private[i][j];
           }
         }
       }
       if (modepdf) 
       {
         for (__int__ i = 0; i < nSamples; ++i) {
           for (__int__ j = 0; j < nModes+1; ++j) {
             for (__int__ k = 0; k < maxBins[i]; ++k) {
               // Do the same summation for the mode arrays
 #pragma omp atomic
               samplePDF_mode_array[i][j][k] += samplePDF_mode_array_private[i][j][k];
             }
           }
         }
      }

     // Delete each thread's private array (still in OMP block)
     for (__int__ i = 0; i < nSamples; ++i) {
       delete[] samplePDF_array_private[i];
       delete[] samplePDF_w2_array_private[i];
     }
     delete[] samplePDF_array_private;
     delete[] samplePDF_w2_array_private;

     // Delete each thread's private mode array (we're still in a OMP block)
     // This have been assigned even if we aren't running with modepdf
     if (modepdf) {
       for (__int__ i = 0; i < nSamples; ++i) 
       {
         for (__int__ j = 0; j < nModes+1; ++j) {
           delete[] samplePDF_mode_array_private[i][j];
         }
         delete[] samplePDF_mode_array_private[i];
       }
       delete[] samplePDF_mode_array_private;
     } else {
       delete[] samplePDF_mode_array_private[0][0];
       delete[] samplePDF_mode_array_private[0];
       delete[] samplePDF_mode_array_private;
     }
   } // end the #pragma omp parallel region

   if (modepdf) 
   {
     #pragma omp parallel for
     for (__int__ i = 0; i < nSamples; ++i) 
     {
       if (samplepdfs->At(i) == NULL) continue;
       for (__int__ j = 0; j < nModes+1; ++j)
       {
         for (__int__ k = 0; k < maxBins[i]-__TH2PolyOverflowBins__; ++k) 
         {
           ((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->SetBinContent(k, samplePDF_mode_array[i][j][k]+((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->GetBinContent(k));
           //reset arrays
           samplePDF_mode_array[i][j][k] = 0.0;
         }
         for (__int__ k = -__TH2PolyOverflowBins__; k < 0; ++k) 
         {
           ((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->SetBinContent(k, samplePDF_mode_array[i][j][maxBins[i]+k]+((TH2Poly*)((TObjArray*)(samplemodepdfs->At(i)))->At(j))->GetBinContent(k));
           //reset arrays
           samplePDF_mode_array[i][j][maxBins[i]+k] = 0.0;
         }
       }
     }
   }
 #if DEBUG > 0
   std::cout << "**************PRINTING DETAILED SAMPLE INFORMATION AFTER LOOP************" << std::endl;
   std::string FileName;
   if (FitManager != NULL) {
     FileName = FitManager->getOutputFilename();
     // Strip out the ending .root extension
     FileName = FileName.substr(0, FileName.find(".root"))+"_diagnostics.root";
   } else {
     FileName = "Diagnostics.root";
   }
   // Get the TFile
   TFile *DiagFile = new TFile(FileName.c_str(), "RECREATE");
   double MCOnlyBin = 0.0;
   double MCOnlyUn = 0.0;
   double POTOnlyBin = 0.0;
   double POTOnlyUn = 0.0;
   double FluxOnlyBin = 0.0;
   double FluxOnlyUn = 0.0;
   double DetOnlyBin = 0.0;
   double DetOnlyUn = 0.0;
   double XsecOnlyBin = 0.0;
   double XsecOnlyUn = 0.0;
   double NDCovOnlyBin = 0.0;
   double NDCovOnlyUn = 0.0;
   double AllOnlyBin = 0.0;
   double AllOnlyUn = 0.0;
   // Loop over the samples
   for (int i = 0; i < nSamples; ++i)
   {
     if (samplepdfs->At(i) == NULL || NoOverflowIntegral(((TH2Poly*)(samplepdfs->At(i)))) == 0) continue;
     std::cout << SampleName[i] << ":" << std::endl;
     std::cout << "  MC only:             " << NoOverflowIntegral(MCOnly[i]) << std::endl;
     MCOnlyBin += NoOverflowIntegral(MCOnly[i]) ;
     std::cout << "  MC only nobin:       " << OverflowIntegral(MCOnly[i]) << std::endl;
     MCOnlyUn += OverflowIntegral(MCOnly[i]) ;

     std::cout << "  POT only:            " << NoOverflowIntegral(POTOnly[i]) << std::endl;
     POTOnlyBin += NoOverflowIntegral(POTOnly[i]) ;
     std::cout << "  POT only nobin:      " << OverflowIntegral(POTOnly[i]) << std::endl;
     POTOnlyUn += OverflowIntegral(POTOnly[i]) ;

     std::cout << "  Flux only:           " << NoOverflowIntegral(FluxOnly[i]) << std::endl;
     FluxOnlyBin += NoOverflowIntegral(FluxOnly[i]) ;
     std::cout << "  Flux only nobin:     " << OverflowIntegral(FluxOnly[i]) << std::endl;
     FluxOnlyUn += OverflowIntegral(FluxOnly[i]) ;

     std::cout << "  Xsec only:           " << NoOverflowIntegral(XsecOnly[i]) << std::endl;
     XsecOnlyBin += NoOverflowIntegral(XsecOnly[i]) ;
     std::cout << "  Xsec only nobin:     " << OverflowIntegral(XsecOnly[i]) << std::endl;
     XsecOnlyUn += OverflowIntegral(XsecOnly[i]) ;

     std::cout << "  Det only:            " << NoOverflowIntegral(DetOnly[i]) << std::endl;
     DetOnlyBin += NoOverflowIntegral(DetOnly[i]) ;
     std::cout << "  Det only nobin:      " << OverflowIntegral(DetOnly[i]) << std::endl;
     DetOnlyUn += OverflowIntegral(DetOnly[i]) ;

     std::cout << "  ND Cov only:         " << NoOverflowIntegral(NDCovOnly[i]) << std::endl;
     NDCovOnlyBin += NoOverflowIntegral(NDCovOnly[i]) ;
     std::cout << "  ND Cov nobin:        " << OverflowIntegral(NDCovOnly[i]) << std::endl;
     NDCovOnlyUn += OverflowIntegral(NDCovOnly[i]) ;

     std::cout << "  All comb:            " << NoOverflowIntegral(AllOnly[i]) << std::endl;
     AllOnlyBin += NoOverflowIntegral(AllOnly[i]) ;
     std::cout << "  All comb nobin:      " << OverflowIntegral(AllOnly[i]) << std::endl;
     AllOnlyUn += OverflowIntegral(AllOnly[i]) ;

     // Write
     DiagFile->cd();

     // Do this or else root doesn't draw these (it thinks entries == 0 so doesn't draw anything). This is because we're using AddBinContent rather than Fill
     MCOnly[i]->SetEntries(1);
     POTOnly[i]->SetEntries(1);
     FluxOnly[i]->SetEntries(1);
     XsecOnly[i]->SetEntries(1);
     DetOnly[i]->SetEntries(1);
     NDCovOnly[i]->SetEntries(1);
     AllOnly[i]->SetEntries(1);
     NuMuPmu[i]->SetEntries(1);

     MCOnly[i]->Write();
     POTOnly[i]->Write();
     FluxOnly[i]->Write();
     XsecOnly[i]->Write();
     DetOnly[i]->Write();
     NDCovOnly[i]->Write();
     AllOnly[i]->Write();
     NuMuPmu[i]->Write();
   }
   std::cout << "Totals: " << std::endl;
   std::cout << "  MCOnly bin:      " << MCOnlyBin << std::endl;
   std::cout << "  MCOnly un:       " << MCOnlyUn << std::endl;
   std::cout << "  POTOnly bin:     " << POTOnlyBin << std::endl;
   std::cout << "  POTOnly un:      " << POTOnlyUn << std::endl;
   std::cout << "  FluxOnly bin:    " << FluxOnlyBin << std::endl;
   std::cout << "  FluxOnly un:     " << FluxOnlyUn << std::endl;
   std::cout << "  XsecOnly bin:    " << XsecOnlyBin << std::endl;
   std::cout << "  XsecOnly un:     " << XsecOnlyUn << std::endl;
   std::cout << "  DetOnly bin:     " << DetOnlyBin << std::endl;
   std::cout << "  DetOnly un:      " << DetOnlyUn << std::endl;
   std::cout << "  NDCovOnly bin:   " << NDCovOnlyBin << std::endl;
   std::cout << "  NDCovOnly un:    " << NDCovOnlyUn << std::endl;
   std::cout << "  AllOnly bin:     " << AllOnlyBin << std::endl;
   std::cout << "  AllOnly un:      " << AllOnlyUn << std::endl;

   DiagFile->Write();
   DiagFile->Close();
   delete DiagFile;
   std::cout << "With xsec systematics: " << std::endl;
   XsecCov->printNominalCurrProp();
 #endif
   //afterclock.Stop();
   //std::cout << "afterclock took " << afterclock.RealTime() << std::endl;
} // end function
#endif


// *******************************************
// Calculate the cross-section weight
double samplePDFND::CalcXsecWeight(const int i) {
// *******************************************

   double xsecw = 1.;

   // Should apply the "once only" weight
   // These are weights that apply to the nominal MC only once
   xsecw *= NDEve[i].weight;
 #if DEBUG > 0
   if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
     std::cerr << "NEGATIVE XSEC ONE-TIME WEIGHT" << std::endl;
     std::cerr << "one-time weight = " << NDEve[i].weight << std::endl;
     std::cerr << "Setting to zero!" << std::endl;
     XsecCov->printNominalCurrProp();
   }
 #if DEBUG == 2
   std::cout << "===============\nPRINTING EVENT BY EVENT WEIGHTS\n===============" << std::endl;
   std::cout << "EVENT " << i << std::endl;
   std::cout << "One-time weight = " << xsecw << std::endl;
 #endif
 #endif

   // Get the spline weights
   // These are different for CPU and GPU code because we don't call TSpline3->Eval for GPU
   xsecw *= CalcXsecWeight_Spline(i);
 #if DEBUG > 0
   if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
     std::cerr << "NEGATIVE XSEC SPLINE WEIGHT" << std::endl;
     std::cerr << "spline weight = " << CalcXsecWeight_Spline(i) << std::endl;
     std::cerr << "Setting to zero!" << std::endl;
     XsecCov->printNominalCurrProp();
   }
 #if DEBUG == 2
   std::cout << "Spline weight = " <<  CalcXsecWeight_Spline(i) << std::endl;
 #endif

 #endif

   // Get the normalisatiSon weights
   // These are the same for CPU and GPU code
   //std::cout<<"CalcXsecWeight at "<<i<<" = "<<CalcXsecWeight_Norm(i)<<std::endl;
   xsecw *= CalcXsecWeight_Norm(i);
   //std::cout<<"xsecw at "<<i<<" = "<<xsecw<<std::endl;
#if DEBUG > 0
   if (xsecw < 0 || xsecw > __LARGE_WEIGHT__) {
     std::cerr << "NEGATIVE XSEC NORM WEIGHT" << std::endl;
     std::cerr << "norm weight = " << CalcXsecWeight_Norm(i) << std::endl;
     std::cerr << "Setting to zero!" << std::endl;
     XsecCov->printNominalCurrProp();
   }
 #endif

   // Do a check on the cross-section weight
   // If negative (maybe we've gone outside the spline interpolation), set to zero
   if (xsecw < 0) {
     xsecw = 0;
     // If large, cap at 100
   } else if (xsecw > __LARGE_WEIGHT__) {
     xsecw = __LARGE_WEIGHT__;
   }

 #if DEBUG == 2
   std::cout << "Total xsec weight = " << xsecw << std::endl;
 #endif

   return xsecw;
 }



#ifdef CUDA
// *********************************************
// Calculate cross-section weights for using GPU
double samplePDFNDGPU::CalcXsecWeight_Spline(const int i) {
  // *********************************************
  double xsecw = 1.;
  // Loop over the spline parameters and get their responses
  for (int id = 0; id < nSplineParams; id++) {
#ifdef DEBUG
    if (std::isnan(xsecw)) {
      std::cerr << "Found nan xsecw, event " << i << ", spline " << splineParsNames[id] << std::endl;
      throw;
    }
#endif
    xsecw *= splineMonolith->cpu_weights[i*nSplineParams+id];
    //KS: if weight is less then 0 there is no need to make spline any further
    //we can just return 0 and save some precious CPU time
    if (xsecw <= 0)
    {
        return 0;
    //std::cerr << "Found negative xsecw, event " << i << ", spline " << splineParsNames[id] << " = " << xsecw << std::endl;
    }
  } // end the k loop

  // If we're debugging check CPU and GPU weight
#ifdef DEBUG_DUMP
  CompareCPU_GPU_Splines(i);
#endif
  return xsecw;
} // end CalcXsecWeight_Spline

#else
  #if USE_SPLINE >= USE_TF1
// ***************************************************************************
// Calculate the TF1 weight for one event i
double samplePDFND::CalcXsecWeight_Spline(const int i) {
// ***************************************************************************
// The returned cross-section weight
   double xsecw = 1.0;
   // Loop over the spline parameters and get their responses
   for (int id = 0; id < nSplineParams; ++id) {
     // Check that the TF1 exists
     if (xsecInfo[i].GetFunc(id)) {
       // Get the global index in cross-section parameterisation for this spline parameter
       int GlobalIndex = splineParsIndex[id];
       // Then get the variation we want to reweight to
       double xvar = XsecCov->calcReWeight(GlobalIndex);
       // Calculate xsecw by evaluating the TF1
       xsecw *= xsecInfo[i].Eval(id, xvar);
  #ifdef DEBUG
       if (std::isnan(xsecw)) {
         std::cerr << "Found nan xsecw, event " << i << ", spline " << splineParsNames[id] << std::endl;
         throw;
       }
  #endif
       //KS: if weight is less then 0 there is no need to make spline any further
       //we can just return 0 and save some precious CPU time
       if (xsecw <= 0) 
       {
         return 0;
         //std::cerr << "Found negative xsecw, event " << i << ", spline " << splineParsNames[id] << " = " << xsecw << std::endl;
       }
     }
   } // End the for loop
   return xsecw;
 }

  #else
// ***************************************************************************
// Calculate the spline weight for one event i
double samplePDFND::CalcXsecWeight_Spline(const int i) {
// ***************************************************************************
   double xsecw = 1.0;
   // Loop over the spline parameters and get their responses
   for (int spline = 0; spline < nSplineParams; ++spline) {
     xsecw *= FastSplineEval((xsecInfo[i].GetFunc(spline)), spline);
 #ifdef DEBUG
     if (std::isnan(xsecw)) {
       std::cerr << "Found nan xsecw, event " << i << ", spline " << splineParsNames[spline] << std::endl;
       throw;
     }
 #endif
       //KS: if weight is less then 0 there is no need to make spline any further
       //we can just return 0 and save some precious CPU time
       if (xsecw <= 0) 
       {
         return 0;
         //std::cerr << "Found negative xsecw, event " << i << ", spline " << splineParsNames[id] << " = " << xsecw << std::endl;
       }
   } // End the for loop
   return xsecw;
 }
  #endif
#endif

// ***************************************************************************
// Calculate the normalisation weight for one event
double samplePDFND::CalcXsecWeight_Norm(const int i) {
// ***************************************************************************

   double xsecw = 1.0;

    for(int iParam = 0; iParam < NDEve[i].nxsec_norm_pointers; ++iParam)
    {
       xsecw *= *(NDEve[i].xsec_norm_pointers[iParam]);
    #ifdef DEBUG
       if (TMath::IsNaN(xsecw)) std::cout << "iParam=" << iParam << "xsecweight=nan from norms" << std::endl;
    #endif
    }
   return xsecw;
}



// ***************************************************************************
// Here we set the cross-section covariance matrix for the samplePDF
// Essentially involves: Checking the input root file for the covarianceXsec class
//                       Finding what parameters are normalisation parameters
//                       Finding what parameters are spline parameters (actually doesn't happen here in the code, see instead fillReweightingBins()
//                       Once all the normalisation parameters are found, set their modes (e.g. CC coherent normalisation should only apply to coherent interaction modes)
//                       Most of this is done automatically now, other than setting what xsec model it is by looking at the size (could probably be made even better by string comparison on the input root file name!) and the modes (could also probably be done by string comparison, e.g. if (string.find("CCQE") != std::string::npos) norm_mode = kMaCh3_CCQE)
void samplePDFND::setXsecCov(covarianceXsec * const xsec_cov) {
// ***************************************************************************

  // Set the XsecCov var
  XsecCov = xsec_cov;

  // Get the number of normalisation parameters
  nxsec_norm_modes = XsecCov->GetNumNearNormParams();

  // Get the infomation for the normalistion parameters
  xsec_norms = XsecCov->GetNearNormPars();

  // Number of spline parameters
  nSplineParams   = XsecCov->GetNumNearSplineParams();
  splineParsIndex = XsecCov->GetNearSplineParsIndex();
  splineParsNames = XsecCov->GetNearSplineParsNames();
  splineFileParsNames = XsecCov->GetNearSplineFileParsNames();
  // Number of unique splines (neutrino and anti-neutrino can share spline name but have a different parameter in covarianceXsec2015)
  nSplineParamsUniq   = XsecCov->GetNumSplineParamsUniq();
  splineParsUniqIndex = XsecCov->GetSplineParsUniqIndex();
  splineParsUniqNames = XsecCov->GetSplineParsUniqNames();

  // Search through which parameters exist in splinePars and don't exist in splineParsUniq
  int nSplineParamsShare    = XsecCov->GetNumSplineParamsShare();
  std::vector<int> splineParsShareIndex  = XsecCov->GetSplineParsShareIndex();
  std::vector<int> splineParsShareToUniq = XsecCov->GetSplineParsShareToUniq();
  std::vector<std::string> splineParsShareNames  = XsecCov->GetSplineParsShareNames();

  // Find the functional parameters (Eb and Alpha_q3)
  nFuncParams = XsecCov->GetNumNearFuncParams();
  funcParsNames = XsecCov->GetNearFuncParsNames();
  funcParsIndex = XsecCov->GetNearFuncParsIndex();

  // Keep the bit-field first entry
  for (int i = 0; i < XsecCov->GetNumParams(); i++) {
    xsecBitField.push_back(XsecCov->GetXSecParamID(i,1));
  }

  // Output the normalisation parameters as a sanity check!
  std::cout << "Normalisation parameters: " << std::endl;
  for (int i = 0; i < nxsec_norm_modes; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << xsec_norms[i].index << "|" << std::setw(30) << xsec_norms[i].name << "|" <<  std::setw(5);
    for(unsigned j = 0; j < xsec_norms[i].modes.size(); j++){
      std::cout<< xsec_norms[i].modes[j];
    }
    std::cout<< std::endl;

  }

  std::cout << "Total spline parameters: " << std::endl;
  for (int i = 0; i < nSplineParams; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << splineParsIndex.at(i) << "|" << std::setw(30) << splineParsNames.at(i) << std::endl;
  }

  std::cout << "Unique spline parameters: " << std::endl;
  for (int i = 0; i < nSplineParamsUniq; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << splineParsUniqIndex.at(i) << "|" << std::setw(30) << splineParsUniqNames.at(i) << std::endl;
  }

  std::cout << "Shared spline parameters: " << std::endl;
  for (int i = 0; i < nSplineParamsShare; ++i) {
    std::cout << std::setw(5) << i << "|" << std::setw(5) << splineParsShareIndex.at(i) << "|" << std::setw(30) << splineParsShareNames.at(i) << std::endl;
  }

#if USE_SPLINE < USE_TF1
  // Initialise
  SplineInfoArray = new FastSplineInfo[nSplineParams];
  for (int i = 0; i < nSplineParams; ++i) {
    SplineInfoArray[i].nPts = -999;
    SplineInfoArray[i].xPts = NULL;
    SplineInfoArray[i].CurrSegment = -999;
  }
#endif
} // End setting of setXsecCov for normalisation parameters


// ***************************************************************************
// Get individual sample likelihood
double samplePDFND::getSampleLikelihood(int isample) {
// ***************************************************************************

  double negLogLsample = 0.;

  // Skip disabled samples
  if (datapdfs->At(isample) == NULL || samplepdfs->At(isample) == NULL) return 0.0;

  for (int j = 1; j < maxBins[isample] - __TH2PolyOverflowBins__; ++j)
  {
      double mc = samplePDF_array[isample][j];
      double dat = samplePDF_data_array[isample][j];
      double w2 = samplePDF_w2_array[isample][j];
      negLogLsample += getTestStatLLH(dat, mc, w2);
  }
  return negLogLsample;
}

// ***************************************************************************
// Get the likelihood
double samplePDFND::getLikelihood() {
// ***************************************************************************

  double negLogL = 0.;

  int isample = 0;
  int j       = 0;

  // This only brings speed increase of 2 for 8 processors but might as well leave in
  // If DEBUG > 0 we also verify that the multi-threaded version equates to single thread (see below)
#ifdef MULTITHREAD
#pragma omp parallel for reduction(+:negLogL) private(isample, j)
#endif
  for (isample = 0; isample < nSamples; ++isample) 
  {
    double negLogLsample = 0.;

    if (datapdfs->At(isample) == NULL) continue;

    //KS: "maxBins[isample] - __TH2PolyOverflowBins__" is identical to "((TH2Poly*)datapdfs->At(isample))->GetNumberOfBins()+1", but make as less relying on root, and could be tiny bit faster
    for (j = 1; j < maxBins[isample] - __TH2PolyOverflowBins__; ++j)
    {
      double mc = samplePDF_array[isample][j];
      double dat = samplePDF_data_array[isample][j];
      double w2 = samplePDF_w2_array[isample][j];
      negLogLsample += getTestStatLLH(dat, mc, w2);

#if DEBUG == 3
      //KS: TH2Poly has excesive bins to put it liglty GetBinContent(j) give (j-1)-th that's why with GetBins() we have to acces (j-1).
      std::cout<< "Sample "<<SampleName[isample]<<std::endl;
      TH2PolyBin* polybinTemp = (TH2PolyBin*)(((TH2Poly*)samplepdfs->At(isample))->GetBins()->At(j-1)->Clone());
      std::cout<< "Bin "<< j<<" GetXMin "<< polybinTemp->GetXMin()<<" GetXMax "<< polybinTemp->GetXMax() <<std::endl;
      std::cout<< "Bin "<< j<<" GetYMin "<< polybinTemp->GetYMin()<<" GetYMax "<< polybinTemp->GetYMax() <<std::endl;
      std::cout<< "mc "<<mc<<" data "<<dat<<" w2 "<<w2<<" negLogLsample "<< getTestStatLLH(dat, mc, w2)<<std::endl;
#endif
    }
    negLogL += negLogLsample;
  } // end sample for loop


  // *************************************************************
  // THREAD SAFE VERIFICATION
  // TESTS IF MULTI-THREADED MODE GIVES DIFFERENT LOGL TO SINGLE-THREADED MODE
  // *************************************************************
#ifdef DEBUG
  // Test single-thread version
  double negLogLsingle = 0.;
  for (isample = 0; isample < nSamples; isample++)
  {
    if (datapdfs->At(isample) == NULL || samplepdfs->At(isample) == NULL) continue;
    for(j = 1; j < ((TH2Poly*)datapdfs->At(isample))->GetNumberOfBins()+1; j++)
    {
      double mc = samplePDF_array[isample][j];
      double dat = samplePDF_data_array[isample][j];
      // Can use normal Poisson, Barlow-Beeston or "IceCube" LLH
      double w2 = samplePDF_w2_array[isample][j];
      negLogLsingle += getTestStatLLH(dat, mc, w2);
    }
  } // end sample for loop

  // Arbitrarily defined 1.E-6, but should be sufficient
  if (fabs(negLogL - negLogLsingle) > 1.E-6) {
    std::cout << "negLogL_MP - negLogL_SP > 1.E-6!" << std::endl;
    std::cout << fabs(negLogL - negLogLsingle) << std::endl;
  }
  // *************************************************************
  // END OF THREAD SAFE VERIFICATION
  // *************************************************************
#endif

  return negLogL;
}

// *************************
// Calculate the Barlow-Beeston likelhood contribution from MC statistics
// Assumes the beta scaling parameters are Gaussian distributed
// Follows arXiv:1103.0354 section 5 and equation 8, 9, 10, 11 on page 4/5
// Essentially solves equation 11
// data is data, mc is mc, w2 is Sum(w_{i}^2) (sum of weights squared), which is sigma^2_{MC stats}
double samplePDFND::getTestStatLLH(double data, double mc, double w2) {
// *************************

  // Need some MC
  if (mc == 0) return 0.0;

  // The MC used in the likeliihood calculation
  // Is allowed to be changed by Barlow Beeston beta parameters
  double newmc = mc;

  // Not full Barlow-Beeston or what is referred to as "light": we're not introducing any more parameters
  // Assume the MC has a Gaussian distribution around generated
  // As in https://arxiv.org/abs/1103.0354 eq 10, 11

  // The penalty from MC statistics using Barlow-Beeston
  double penalty = 0;
  if (fTestStatistic == kBarlowBeeston) {
    // Barlow-Beeston uses fractional uncertainty on MC, so sqrt(sum[w^2])/mc
    double fractional = sqrt(w2)/mc;
    // -b/2a in quadratic equation
    double temp = mc*fractional*fractional-1;
    // b^2 - 4ac in quadratic equation
    double temp2 = temp*temp + 4*data*fractional*fractional;
    if (temp2 < 0) {
      std::cerr << "Negative square root in Barlow Beeston coefficient calculation!" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }  
    // Solve for the positive beta
    double beta = (-1*temp+sqrt(temp2))/2.;
    newmc = mc*beta;
    // And penalise the movement in beta relative the mc uncertainty
    if (fractional > 0) penalty = (beta-1)*(beta-1)/(2*fractional*fractional);
    else penalty = 0;
  }

  // Calculate the new Poisson likelihood
  // For Barlow-Beeston newmc is modified, so can only calculate Poisson likelihood after Barlow-Beeston
  // For the Poisson likelihood, this is just the usual calculation
  // For IceCube likelihood, we calculate it later
  double stat = 0;
  // All likelihood calculations may use the bare Poisson likelihood, so calculate here
  if (data == 0) stat = newmc;
  else if (newmc > 0) stat = newmc-data+data*TMath::Log(data/newmc);

  // Also try the IceCube likelihood
  // It does not modify the MC content
  // https://arxiv.org/abs/1901.04645
  // Argüelles, C.A., Schneider, A. & Yuan, T. J. High Energ. Phys. (2019) 2019: 30. https://doi.org/10.1007/JHEP06(2019)030
  // We essentially construct eq 3.16 and take the logarithm
  // in eq 3.16, mu is MC, sigma2 is w2, k is data
  if (fTestStatistic == kIceCube) {
    // If there for some reason is 0 mc uncertainty, return the Poisson LLH
    if (w2 == 0) return stat;

    // Reset the penalties if there is mc uncertainty
     stat = 0.0;
     penalty = 0.0;
     // Auxillary variables
     long double b = mc/w2;
     long double a = mc*b+1;
     long double k = data;
     // Use C99's implementation of log of gamma function to not be C++11 dependent
     stat = -1*(a * logl(b) + lgammal(k+a) - lgammal(k+(long double)1) - ((k+a)*log1pl(b)) - lgammal(a));
   }

   // Return the statistical contribution and penalty
   return stat+penalty;
 }

// *************************
TH1* samplePDFND::getPDF(int Selection) {
// *************************

  if (samplepdfs->At(Selection) == NULL) return NULL;

  //First Reset Histogram
  if (ndims[Selection] == 1) {
    ((TH1*)samplepdfs->At(Selection))->Reset("");
    ((TH1*)samplepdfs->At(Selection))->Fill(0.0, 0.0);
  } else if (ndims[Selection] == 2) {
    ((TH2Poly*)samplepdfs->At(Selection))->Reset("");
    ((TH2Poly*)samplepdfs->At(Selection))->Fill(0.0, 0.0, 0.0);
  }

  //Fill histogram from master array
  for (__int__ j = 0; j < maxBins[Selection]-__TH2PolyOverflowBins__; ++j)
  {
    ((TH2Poly*)samplepdfs->At(Selection))->SetBinContent(j, samplePDF_array[Selection][j]);
  }
  for (__int__ j =-__TH2PolyOverflowBins__; j < 0; ++j)
  {
    ((TH2Poly*)samplepdfs->At(Selection))->SetBinContent(j, samplePDF_array[Selection][maxBins[Selection]+j]);
  }

  return (TH1*)samplepdfs->At(Selection);
}

// *************************
TH2Poly* samplePDFND::getW2(int Selection) {
// *************************

  if(W2Hist.at(Selection) == NULL ) return NULL;

  //First Reset Histogram
  W2Hist[Selection]->Reset("");
  W2Hist[Selection]->Fill(0.0, 0.0, 0.0);
  
  //Fill histogram from master array
  for (__int__ j = 0; j < maxBins[Selection]-__TH2PolyOverflowBins__; ++j)
  {
    W2Hist[Selection]->SetBinContent(j, samplePDF_w2_array[Selection][j]);
  }
  for (__int__ j =-__TH2PolyOverflowBins__; j < 0; ++j)
  {
    W2Hist[Selection]->SetBinContent(j, samplePDF_w2_array[Selection][maxBins[Selection]+j]);
  }

  return (TH2Poly*)W2Hist.at(Selection);
}

// *************************
// Add data from a vector of doubles
// Only OK for TH1D
void samplePDFND::addData(std::vector<double> &dat, int index) {
// *************************

   if (datapdfs->At(index) == NULL) {
     std::cerr << "Failed to replace data because selection " << SampleName[index] << "(" << index << ")" <<  " is not enabled" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   // Check our data is TH1D
   if (ndims[index] != 1) {
     std::cerr << "Failed to replace data because datapdf is not TH1D" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   // Reset the current data
   ((TH1*)datapdfs->At(index))->Reset();

   // Replace the data with new data
   for (int i = 0; i < int(dat.size()); i++) {
     ((TH1*)datapdfs->At(index))->Fill(dat[i]);
   }

   return;
 }

// *************************
// Add data from a vector of vector of doubles
void samplePDFND::addData(std::vector< vector <double> > &dat, int index) {
// *************************

   if (datapdfs->At(index) == NULL) {
     std::cerr << "Failed to replace data because selection " << SampleName[index] << "(" << index << ")" <<  " is not enabled" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   // Check our data is TH2D
   if (ndims[index] != 2) {
     std::cerr << "Failed to replace data because datapdf is not TH1D" << std::endl;
     std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
     throw;
   }

   ((TH2Poly*)datapdfs->At(index))->Reset("");//An option is required here in the reset function for TH2Polys so have to have ""

   for(int i=0; i<int(dat.size()); i++) {
     ((TH2Poly*)datapdfs->At(index))->Fill(dat[0][i], dat[1][i]);
   }

   return;
}

// *************************
// Add data from TH1D
void samplePDFND::addData(TH1D* binneddata, int index){
  // *************************

   if (datapdfs->At(index)) {
     datapdfs->RemoveAt(index);
     datapdfs->AddAt(binneddata,index);
   } else {
     datapdfs->AddAt(binneddata,index);
   }

   return;
}

// *************************
// Add data from TH2D
void samplePDFND::addData(TH2Poly* binneddata, int index){
// *************************

   if(datapdfs->At(index)) {
     datapdfs->RemoveAt(index);
     datapdfs->AddAt(binneddata,index);
   } else {
     datapdfs->AddAt(binneddata,index);
   }

   return;
}

// *************************
// Get the event rate
double samplePDFND::getEventRate(int index) {
// *************************

   if(samplepdfs->At(index)) {
     return NoOverflowIntegral((TH2Poly*)getPDF(index));
   } else {
     return 0;
   }
 }



#if USE_SPLINE >= USE_TF1
// ***************************************************************************
// TGraph** because each xsecgraph points to a TObjArray and each collection of xsecgraph has nSplines entries
//
// #########################
// WARNING THIS NEEDS TO MATCH THE ENUMERATION IN chain->SetBranchAddress IN samplePDFND::fillReweightingBins()
void samplePDFND::SetSplines(TGraph** &xsecgraph, const int EventNumber) {
// ***************************************************************************

   xsecInfo[EventNumber].SetSplineNumber(nSplineParams);

 #ifdef DEBUG_TF1
   std::stringstream ss;
   ss << EventNumber;
   dumpfile->cd();
 #endif

   // Now loop over and set the splines that we have
   // This is definitely the set of unique parameters
   // When we do TSpline3->Eval(paramVal) we don't want paramVal to be the same though!
   for (int j = 0; j < nSplineParams; ++j) {

     // The TF1 object we build from fitting the TGraph
     TF1 *Fitter = NULL;
 #if USE_SPLINE == USE_TF1_red
     TF1_red *tf1_red = NULL;
 #endif

     // For a valid spline we require it not be null and have more than one point
     if (xsecgraph[j] && xsecgraph[j]->GetN() > 1) {

       // Set a unique name
       std::stringstream SplineTitle;
       SplineTitle << EventNumber << "_" << xsecgraph[j]->GetName();

       // For 2p2h shape C and O we can't fit a polynomial: try a linear combination of two linear functions around 0
       if (j == 3 || j == 4) {
        Fitter = new TF1(SplineTitle.str().c_str(), "(x<=0)*(1+[0]*x)+(x>0)*([1]*x+1)", xsecgraph[j]->GetX()[0], xsecgraph[j]->GetX()[xsecgraph[j]->GetN()-1]);
        // Fit 5hd order polynomial for all other parameters
       } else {
        Fitter = new TF1(SplineTitle.str().c_str(), "1+[0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x+[4]*x*x*x*x*x", xsecgraph[j]->GetX()[0], xsecgraph[j]->GetX()[xsecgraph[j]->GetN()-1]);
       }
       // Fit the TF1 to the graph
       xsecgraph[j]->Fit(Fitter, "Q0");

 #if USE_SPLINE == USE_TF1_red
       // Make the reduced TF1 if we want
       tf1_red = new TF1_red(Fitter, j);
 #endif

       // Debug the TF1 by comparing it to the TSpline3
 #ifdef DEBUG_TF1
       // Make the comparison spline
       std::string name = xsecgraph[j]->GetName();
       name = ss.str()+"_"+name;
       TSpline3 *spline = new TSpline3(xsecgraph[j]->GetName(), xsecgraph[j]);
       xsecgraph[j]->SetNameTitle((name+"_gr").c_str(), (name+"_gr").c_str());
       Fitter->SetNameTitle((name+"_tf").c_str(), (name+"_tf").c_str());
       spline->SetNameTitle((name+"_sp").c_str(), (name+"_sp").c_str());

       // The total xrange
       double xrange = xsecgraph[j]->GetX()[xsecgraph[j]->GetN()-1] - xsecgraph[j]->GetX()[0];
       bool checkit = false;
       double point = -999;
       // Scan 100 points and look for differences in weights
       const int npoints = 100;
       for (int z = 0; z < npoints; ++z) {
        double val = xsecgraph[j]->GetX()[0] + z*xrange/double(npoints);
        double fitted = Fitter->Eval(val);
        double splined = spline->Eval(val);
        if (fabs(fitted/splined-1) > 0.05) {
          checkit = true;
          point = val;
          break;
        }
       }

      // If we've found a dodgy event write to file
      if (checkit == true) {
      TCanvas *temp = new TCanvas("canv","canv", 800,800);
      temp->SetName((name+"_canv").c_str());
      temp->SetTitle((name+"_canv").c_str());
      xsecgraph[j]->SetLineColor(kBlack);
      xsecgraph[j]->SetLineWidth(2);
      xsecgraph[j]->SetLineStyle(kSolid);
      xsecgraph[j]->SetFillStyle(0);
      xsecgraph[j]->SetFillColor(0);
      Fitter->SetLineColor(kRed);
      Fitter->SetLineStyle(kDashed);
      Fitter->SetLineWidth(2);
      spline->SetLineColor(kOrange);
      spline->SetLineStyle(kSolid);
      spline->SetLineWidth(2);

      temp->cd();
      xsecgraph[j]->Draw();
      TLine *line = new TLine(point, 0, point, 20);
      line->SetLineColor(kRed);
      line->SetLineWidth(2);
      xsecgraph[j]->GetXaxis()->SetTitle("Variation (rel nominal)");
      xsecgraph[j]->GetYaxis()->SetTitle("Weight");
      spline->Draw("same");
      Fitter->Draw("same");
      line->Draw("same");
      TLegend *leg = new TLegend(0.65, 0.85, 1.0, 1.0);
      leg->SetHeader(name.c_str());
      leg->AddEntry(xsecgraph[j], "Graph");
      leg->AddEntry(spline, "TSpline3");
      leg->AddEntry(Fitter, "TF1");
      leg->Draw("same");

      dumpfile->cd();
      xsecgraph[j]->Write();
      Fitter->Write();
      spline->Write();
      temp->Write();

      delete temp;
      delete leg;
      delete line;
       }
       delete spline;
 #endif

       // For events without the jth TGraph, set to NULL
     } else {
       //Fitter = dummy;
       Fitter = NULL;
 #if USE_SPLINE == USE_TF1_red
       tf1_red = NULL;
 #endif
     }

     // Now TF1* Fitter is either a TF1 or a pointer to NULL for the jth parameter
     // Save it to the struct
 #if USE_SPLINE == USE_TF1_red
     xsecInfo[EventNumber].SetFunc(j, tf1_red);
 #elif USE_SPLINE == USE_TF1
     xsecInfo[EventNumber].SetFunc(j, Fitter);
 #endif
   } // end for j loop

   // Delete the TGraph pointers
   for (int i = 0; i < nSplineParams; ++i) {
     delete xsecgraph[i];
   }
 } // end function

 #else
// ***************************************************************************
// TGraph** because each xsecgraph points to a TObjArray and each collection of xsecgraph has nSplines entries
//
// #########################
// WARNING THIS NEEDS TO MATCH THE ENUMERATION IN chain->SetBranchAddress IN samplePDFND::fillReweightingBins()
void samplePDFND::SetSplines(TGraph** &xsecgraph, const int EventNumber) {
// ***************************************************************************

   xsecInfo[EventNumber].SetSplineNumber(nSplineParams);

   // Now loop over and set the splines that we have
   // This is definitely the set of unique parameters
   // When we do TSpline3->Eval(paramVal) we don't want paramVal to be the same though!
   //KS: Tried multithreading here in with one thread it is faster
   for (int j = 0; j < nSplineParams; ++j) {

     // Here's the TSpline3
     TSpline3* spline = NULL;
 #if USE_SPLINE < USE_TSpline3
     __SPLINE_TYPE__ *spline_red = NULL;
 #endif
     // For a valid spline we require it not be null and have more than one point
     if (xsecgraph[j] && xsecgraph[j]->GetN() > 1) {

       // Create the TSpline3* from the TGraph* and build the coefficients
       spline = new TSpline3(xsecgraph[j]->GetName(), xsecgraph[j]);
       spline->SetNameTitle(xsecgraph[j]->GetName(), xsecgraph[j]->GetName());

       // Reduce to use linear spline interpolation for certain parameters
       // Not the most elegant way: use TSpline3 object but set coefficients to zero and recalculate spline points; the smart way (but more human intensive) would be to save memory here and simply not store the zeros at all
       // Get which parameters should be linear from the fit manager
       // Convert the spline number to global xsec parameter
       int globalxsec = splineParsIndex[j];
       for (std::vector<int>::iterator it = linearsplines.begin(); it != linearsplines.end(); ++it) {
         if ((*it) == globalxsec) {
           // Loop over the splines points
           for (int k = 0; k < spline->GetNp(); ++k) {
             Double_t x1, y1, b1, c1, d1, x2, y2, b2, c2, d2 = 0;
             spline->GetCoeff(k, x1, y1, b1, c1, d1);
             spline->GetCoeff(k+1, x2, y2, b2, c2, d2);
             double tempb = (y2-y1)/(x2-x1);
             spline->SetPointCoeff(k, tempb, 0, 0); 
           }
         }
       }

#if USE_SPLINE < USE_TSpline3
       // Make the reduced TSpline3 format and delete the old spline
       spline_red = new __SPLINE_TYPE__(spline, j);
#endif

       // Fill the SplineInfoArray entries with information on each splinified parameter
       if (SplineInfoArray[j].xPts == NULL) {
         // Fill the number of points
#if USE_SPLINE < USE_TSpline3
         SplineInfoArray[j].nPts = spline_red->GetNp();
#else
         SplineInfoArray[j].nPts = spline->GetNp();
#endif

         // Fill the x points
         SplineInfoArray[j].xPts = new double[SplineInfoArray[j].nPts];
         for (int k = 0; k < SplineInfoArray[j].nPts; ++k) {
           double xtemp = -999.99;
           double ytemp = -999.99;
#if USE_SPLINE < USE_TSpline3
           spline_red->GetKnot(k, xtemp, ytemp);
#else
           spline->GetKnot(k, xtemp, ytemp);
#endif
           SplineInfoArray[j].xPts[k] = xtemp;
         }
       }
       // For events without the jth spline, set to NULL
     } else {
       spline = NULL;
#if USE_SPLINE < USE_TSpline3
       spline_red = NULL;
#endif
     }
     // Save the TSpline3* into the struct
#if USE_SPLINE < USE_TSpline3
     xsecInfo[EventNumber].SetFunc(j, spline_red);
#else
     xsecInfo[EventNumber].SetFunc(j, spline);
#endif
   } // end for j loop

   // Delete the TGraph pointers
   for (int i = 0; i < nSplineParams; ++i) {
     delete xsecgraph[i];
   }

 } // end function

// ***************************************************************************
// TGraph** because each xsecgraph points to a TObjArray and each collection of xsecgraph has nSplines entries
//
// #########################
// WARNING THIS NEEDS TO MATCH THE ENUMERATION IN chain->SetBranchAddress IN samplePDFND::fillReweightingBins()
// Reduce number of points in TGraph to 3
void samplePDFND::SetSplines_Reduced(TGraph** &xsecgraph, const int EventNumber) {
// ***************************************************************************

  xsecInfo[EventNumber].SetSplineNumber(nSplineParams);
  // First look at the TGraphs
  for (int i = 0; i < nSplineParams; ++i) {

    // Here's the reduced graph
    TGraph* ReducedGraph = NULL;
    // And the TSpline3 we build from that reduced graph
    TSpline3* Spline = NULL;
#if USE_SPLINE < USE_TSpline3
    __SPLINE_TYPE__ *spline_red = NULL;
#endif

    // If this criteria is true then the spline response is not flat
    if (xsecgraph[i] && xsecgraph[i]->GetN() > 1) {
      // Get the TGraph for this spline parameter
      TGraph *TempGraph = xsecgraph[i];
      // Now look at the x coordinates
      double *xaxis = TempGraph->GetX();
      // And the responses in y

      double *yresp = TempGraph->GetY();
      // The number of points in the original graph
      int nPoints = TempGraph->GetN();
      // The reduced x and y response (set to be half of original)
      double *x_red = new double[nPoints/2+1];
      double *y_red = new double[nPoints/2+1];
      if (nPoints % 2 != 1) {
        std::cerr << "Reduced spline points method only works for odd number of points, sorry" << std::endl;
        std::cerr << "Change the code in " << __FILE__ << ":" << __LINE__ << std::endl;
        std::cerr << "Or revert to using _NOT_ reduced spline code at ND (see SetSplines function in " << __FILE__ << ")" << std::endl;
        throw;
      }
      // Skip every second point
      for (int j = 0; j < (nPoints/2)+1; ++j) {
        // If we're on the last point in the reduced we want to pick out the last point in the total TGraph
        x_red[j] = xaxis[j*2];
        y_red[j] = yresp[j*2];
      }
      // Now make the new TGraphs from these x and y
      ReducedGraph = new TGraph(nPoints/2+1, x_red, y_red);
      ReducedGraph->SetName(xsecgraph[i]->GetName());

      // Now build the TSpline3
      Spline = new TSpline3(ReducedGraph->GetName(), ReducedGraph);
      Spline->SetNameTitle(ReducedGraph->GetName(), ReducedGraph->GetName());
      // Make the reduced TSpline3 format and delete the old spline
#if USE_SPLINE < USE_TSpline3
      spline_red = new __SPLINE_TYPE__(Spline, i);
#endif

      // Fill the SplineInfoArray entries with information on each splinified parameter
      if (SplineInfoArray[i].xPts == NULL) {
        // Fill the number of points
#if USE_SPLINE < USE_TSpline3
        SplineInfoArray[i].nPts = spline_red->GetNp();
#else
        SplineInfoArray[i].nPts = Spline->GetNp();
#endif
        // Fill the x points
        SplineInfoArray[i].xPts = new double[SplineInfoArray[i].nPts];
        for (int k = 0; k < SplineInfoArray[i].nPts; ++k) {
          double xtemp = -999.99;
          double ytemp = -999.99;
#if USE_SPLINE < USE_TSpline3
          spline_red->GetKnot(k, xtemp, ytemp);
#else
          Spline->GetKnot(k, xtemp, ytemp);
#endif
          SplineInfoArray[i].xPts[k] = xtemp;
        }
      }

      // Delete the temporary memory
      delete[] x_red;
      delete[] y_red;
      delete ReducedGraph;
    } else {
      ReducedGraph = NULL;
      Spline = NULL;
    }
    // Finall save the pointer to the reduced spline in the struct
#if USE_SPLINE < USE_TSpline3
    xsecInfo[EventNumber].SetFunc(i, spline_red);
#else
    xsecInfo[EventNumber].SetFunc(i, Spline);
#endif
  } // end loop over xsec parameter

  // Delete the TGraph pointers
  for (int i = 0; i < nSplineParams; ++i) {
    delete xsecgraph[i];
  }
} // end function
#endif

// ***************************************************************************
// Set up our experiment, loop over the events (data and MC), associate cross-section splines with events
// Needs to have addSelection called!
void samplePDFND::fillReweightingBins() {
// ***************************************************************************

  std::cout << "starting fillReweightingBins" <<  std::endl;
  // Check that the covariances are filled
  CheckCovariances();

  FindAdditionalInfo();

  FindNormBins();

  FindNormPointer();

  InitialisePDF();

  //KS: Shrink memory to accepted number of events
  NDEve.resize(nEvents);
#ifndef DEBUG
  //KS: We only needed this to set some varialbes but now we can delete it to save RAM, however do this in standard fit, in debug mode we might need thos varaibles for verbose
  NDEve_Aux.clear();
  std::cout << "Removing auxilary variables not needed during fit, relase " << nEvents*(sizeof(__float__)*2 + sizeof(__int__)*2 + sizeof(bool))/1.E6 << " MB of RAM" << std::endl;
#endif

  std::cerr<<"Function fillReweightingBins is experiment specific however core code uses it"<<std::endl;
  std::cerr<<"Since you haven't implemented it I have to stop it"<<std::endl;
  throw;

} // End fillReweightingBins() function


// ***************************************************************************
// Function to reserve memory for ND objects
void samplePDFND::ReserveMemory(int nEve) {
// ***************************************************************************

  // Number of events which passed the above selection
  nEvents = nEve;
  //KS: nEve is maximal number of events we have in splines, allocate more memro for now, but later we will resize agian to free memory
  NDEve.resize(nEve);
  NDEve_Aux.resize(nEve);

  if (XsecCov) xsecInfo = new XSecStruct<__SPLINE_TYPE__*>[nEve];
}

// ***************************************************************************
// Helper function to find normalisation parameter bins
// This is to avoid hard-coding the normalisations and instead use the input to determine where each normalisation parameter lives in the global parameter index
// e.g. 2p2h normalisation for C/O lives in position x (which we could just hard-code); I prefer to use a search instead in case parameterisations change. It's really up to the end-user!
void samplePDFND::FindNormBins() {
// ***************************************************************************

  std::cout << "nEvents = "<<nEvents<<std::endl; 
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(unsigned int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    std::vector< int > XsecBins;
    NDEve_Aux[iEvent].xsec_norms_bins = XsecBins;

    if (XsecCov) 
    {
      for(std::vector<XsecNorms3>::iterator it = xsec_norms.begin(); it != xsec_norms.end(); ++it)
      {
        bool noH=false;
        bool targetmatch = false;
        if((*it).targets.size() == 0) targetmatch = true;
        for(unsigned itarget = 0; itarget < (*it).targets.size(); itarget++)
        {
          if(fabs((*it).targets.at(itarget) == NDEve_Aux[iEvent].target)) targetmatch = true;
          
          //KS: Target -1 mean all but not hydrogen, this is moslty used by SK
          //if target is -1 and this is only target we apply it, if there are other we 
          //apply only to defined targets
          //Example: element = [ -1 ] - all but not hydrogen
          //element = [ -1 12 16] - only to 12 and 16
          if((*it).targets.at(itarget) == -1 && (*it).targets.size() == 1)
          {
            noH=true;
            if( NDEve_Aux[iEvent].target != 1 ) targetmatch = true;
          }
        }

        if(!targetmatch)continue;
        if(noH && NDEve_Aux[iEvent].target == 1) continue;

        bool horncurrentmatch=false;
        if((*it).horncurrents.size() == 0) horncurrentmatch = true;
        else{
          for(unsigned ihorncurrent=0; ihorncurrent < (*it).horncurrents.size(); ihorncurrent++)//Check if is right horncurrent
          {
            if((*it).horncurrents.at(ihorncurrent) == 1 && !NDEve_Aux[iEvent].isRHC) horncurrentmatch = true;
            if((*it).horncurrents.at(ihorncurrent) == -1 && NDEve_Aux[iEvent].isRHC) horncurrentmatch = true;
          }
        }

        if(!horncurrentmatch) continue;

        bool flavmatch=false;
        if((*it).pdgs.size() == 0) flavmatch = true;
        else{
          for(unsigned ipdg=0; ipdg < (*it).pdgs.size(); ipdg++) //Check if is right neutrino flavour
          {       
            if((*it).pdgs.at(ipdg) == kNue && NDEve_Aux[iEvent].species == kNue) flavmatch = true;
            if((*it).pdgs.at(ipdg) == kNumu && NDEve_Aux[iEvent].species == kNumu) flavmatch = true;
            if((*it).pdgs.at(ipdg) == kNue_bar && NDEve_Aux[iEvent].species == kNue_bar) flavmatch = true;
            if((*it).pdgs.at(ipdg) == kNumu_bar && NDEve_Aux[iEvent].species == kNumu_bar) flavmatch = true;
          }
        }

        if(!flavmatch)continue;

        bool prodflavmatch=false;
        if((*it).preoscpdgs.size() == 0) prodflavmatch = true;
        else{
          for(unsigned ipdg=0; ipdg<(*it).preoscpdgs.size(); ipdg++) //Check if is right neutrino flavour
          {         
            if((*it).preoscpdgs.at(ipdg) == kNue && NDEve_Aux[iEvent].species == kNue) prodflavmatch=true;
            if((*it).preoscpdgs.at(ipdg) == kNumu && NDEve_Aux[iEvent].species == kNumu) prodflavmatch=true;
            if((*it).preoscpdgs.at(ipdg) == kNue_bar && NDEve_Aux[iEvent].species == kNue_bar) prodflavmatch=true;
            if((*it).preoscpdgs.at(ipdg) == kNumu_bar && NDEve_Aux[iEvent].species == kNumu_bar) prodflavmatch=true;
          }
        }

        if(!prodflavmatch) continue;

        bool modematch=false;
        if((*it).modes.size() == 0) modematch = true;
        else{
          for(unsigned imode=0; imode < (*it).modes.size(); imode++)
          {
            if((*it).modes.at(imode) == NDEve[iEvent].mode) modematch = true;
          }
        }

        if(!modematch) continue;

        if((*it).hasKinBounds)
        {
          if((*it).etru_bnd_low >= NDEve_Aux[iEvent].Enu && (*it).etru_bnd_low != -999) continue;
          if((*it).etru_bnd_high < NDEve_Aux[iEvent].Enu && (*it).etru_bnd_high != -999) continue;
          if((*it).q2_true_bnd_low >= NDEve_Aux[iEvent].Q2 && (*it).q2_true_bnd_low != -999) continue;
          if((*it).q2_true_bnd_high < NDEve_Aux[iEvent].Q2 && (*it).q2_true_bnd_high != -999) continue;
        }
        //Now set 'index bin' for each normalisation parameter
        //All normalisations are just 1 bin, so bin = index (where index is just the bin for that normalisation)
        int bin = (*it).index;
        //Check if parameter applies to SK detector
        if(XsecCov->GetXSecParamID(bin,1) & 1) //1 means parameter applies to ND
        {
          NDEve_Aux[iEvent].xsec_norms_bins.push_back(bin);
        }
      } // end iteration over xsec_norms
    } // end if (XsecCov)
  }//end loop over events
  return;
}

// ************************
//KS: Find pointer for each norm dial to reduce impact of covarianceXsec::calcReweight()
void samplePDFND::FindNormPointer() {
// ******************************************
    #ifdef MULTITHREAD
    #pragma omp parallel for
    #endif
    for (unsigned iEvent = 0; iEvent < nEvents; ++iEvent) 
    {
      int counter = 0;
      NDEve[iEvent].nxsec_norm_pointers = NDEve_Aux[iEvent].xsec_norms_bins.size();
      
      const double* xsecPointer[ NDEve[iEvent].nxsec_norm_pointers ];
      for( std::vector < __int__ > ::iterator lit = NDEve_Aux[iEvent].xsec_norms_bins.begin(); lit != NDEve_Aux[iEvent].xsec_norms_bins.end(); lit++)
      {
           xsecPointer[counter]= XsecCov->retPointer(*lit);
           counter += 1;
      }
      std::vector < const double* > Temp_norm_pointers;
      //KS: Terrible code convertng const double array to std::vector
      for(int iParam=0; iParam < NDEve[iEvent].nxsec_norm_pointers; iParam++)
      {
        Temp_norm_pointers.push_back(xsecPointer[iParam]);
      }
      NDEve[iEvent].xsec_norm_pointers =  Temp_norm_pointers;
     }
}


// ******************************************
//KS:  Find Detector bin, hist bin and other useful information using multithreading
void samplePDFND::FindAdditionalInfo() {
// ******************************************

  std::cerr<<"Function LoadSamples is experiment specific however core code uses it"<<std::endl;
  std::cerr<<"Since you haven't implemented it I have to stop it"<<std::endl;
  throw;

}

// ***************************************************************************
// Helper function to check if the covariances have been set
void samplePDFND::CheckCovariances() {
// ***************************************************************************

  // Check if we're doing smart processing
  // This shouldn't happen because (to my knowledge) there is no such root file yet

  if (XsecCov == NULL) {
    std::cerr << "XsecCov covariance not set, can't setup samplePDFND in " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  if (NDDetCov == NULL) {
    std::cerr << "detector covariance not set, can't setup samplePDFND in " << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

}


#if USE_SPLINE < USE_TF1
// *************************
// Only need to do the binary search once per parameter, not once per event!
// Takes down the number of binary searches from 1.2M to 17, ha!
// For ROOT version see root/hist/hist/src/TSpline3.cxx TSpline3::FindX(double)
void samplePDFND::FindSplineSegment() {
// *************************

  // Loop over the splines
  //KS: Tried multithreading here with 48 splines and it is faster with one thread, maybe in future multithreading will be worth revisiting
  for (int i = 0; i < nSplineParams; ++i) 
  {
    const int nPoints = SplineInfoArray[i].nPts;
    const double* xArray = SplineInfoArray[i].xPts;

    if (nPoints == -999 || xArray == NULL) {
      std::cerr << "ERROR" << std::endl;
      std::cerr << "SplineInfoArray[" << i << "] isn't set yet" << std::endl;
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      continue;
      //throw;
    }

    // Get the variation for this reconfigure for the ith parameter
    const double xvar = XsecCov->calcReWeight(splineParsIndex[i]);

#if DEBUG > 0
    std::cout << "Finding spline segment for xsec parameter " << splineParsNames[i] << " at " << xvar << std::endl;
#endif

    // The segment we're interested in (klow in ROOT code)
    int segment = 0;
    int kHigh = nPoints-1;

    //KS: Consider searching first in previous segment and +-1 segment, initial test show it is faster, howver this part of the code takes 10^-6s the improvment is really negligible, but for the bigger number of dials worht considering.
    
    // If the variation is below the lowest saved spline point
    if (xvar <= xArray[0]) {
      segment = 0;
      // If the variation is above the highest saved spline point
    } else if (xvar >= xArray[nPoints-1]) {
      // Yes, the -2 is indeed correct, see TSpline.cxx:814 and //see: https://savannah.cern.ch/bugs/?71651
      segment = kHigh;
      // If the variation is between the maximum and minimum, perform a binary search
    } else {
      // The top point we've got
      int kHalf = 0;
      // While there is still a difference in the points (we haven't yet found the segment)
      // This is a binary search, incrementing segment and decrementing kHalf until we've found the segment
      while (kHigh - segment > 1) {
        // Increment the half-step 
        kHalf = (segment + kHigh)/2;
        // If our variation is above the kHalf, set the segment to kHalf
        if (xvar > xArray[kHalf]) {
          segment = kHalf;
          // Else move kHigh down
        } else {
          kHigh = kHalf;
        }
      } // End the while: we've now done our binary search
    } // End the else: we've now found our point

    if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;

    // Save the segment for the ith parameter
    SplineInfoArray[i].CurrSegment = segment;

#ifdef DEBUG
    if (SplineInfoArray[i].xPts[segment] > xvar && segment != 0) {
      std::cerr << "Found a segment which is _ABOVE_ the variation!" << std::endl;
      std::cerr << "IT SHOULD ALWAYS BE BELOW! (except when segment 0)" << std::endl;
      std::cerr << splineParsNames[i] << std::endl;

      std::cerr << "Found segment   = " << segment << std::endl;
      std::cerr << "Doing variation = " << xvar << std::endl;
      std::cerr << "x in spline     = " << SplineInfoArray[i].xPts[segment] << std::endl;
      for (int j = 0; j < SplineInfoArray[j].nPts; ++j) {
        std::cerr << "    " << j << " = " << SplineInfoArray[i].xPts[j] << std::endl;
      }
      std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
      throw;
    }
#endif
  }
}
#endif

#ifdef CUDA
// *************************************************
// Prepare the weights before the reweight loop starts
// For GPU accelerated code this means passing variations to the GPU, calculating the weights on the GPU, and passing them back to memory
void samplePDFND::PrepareWeights() {
  // *************************************************

  #if USE_SPLINE < USE_TF1
  // With NIWG 2018a the parameters there's a parameter mapping that goes from spline parameter to a global parameter index
  // Find the spline segments
  FindSplineSegment();

  // This way we avoid doing 1.2M+ binary searches on the GPU
  // and literally just multiply lots of numbers together on the GPU without any algorithm
  for (int i = 0; i < nSplineParams; ++i) {
    // Update the values and which segment it belongs to
    vals[i] = XsecCov->calcReWeight(splineParsIndex[i]);
    segments[i] = SplineInfoArray[i].CurrSegment;
  }
  // Call the GPU code through the SplineMonolith class
  splineMonolith->EvalGPU_SepMany_seg(vals, segments);
  #else
  // Feed the parameter variations
  for (int i = 0; i < nSplineParams; ++i) {
    // Update the values and which segment it belongs to
    vals[i] = XsecCov->calcReWeight(splineParsIndex[i]);
  }
  splineMonolith->EvalGPU_TF1(vals);
  #endif

  #ifdef DEBUG_DUMP
  // Write to file
  std::stringstream ss;
  ss << "Reweight_" << nReconf;
  DebugFile->cd();
  TDirectory *NewDir = DebugFile->mkdir(ss.str().c_str());
  NewDir->cd();

  for (int j = 0; j < nSplineParams; ++j) {

    // Save the parameter variation in the title
    std::stringstream paramval;
    paramval << "_" << vals[j];
    std::string restoregpuname = std::string(gpu_weights_plot[j]->GetName());
    std::string restorecpuname = std::string(cpu_weights_plot[j]->GetName());
    std::string restorediffname = std::string(diff_weights_plot[j]->GetName());

    gpu_weights_plot[j]->SetTitle((std::string(gpu_weights_plot[j]->GetTitle())+"_"+paramval.str()).c_str());
    cpu_weights_plot[j]->SetTitle((std::string(cpu_weights_plot[j]->GetTitle())+"_"+paramval.str()).c_str());
    diff_weights_plot[j]->SetTitle((std::string(diff_weights_plot[j]->GetTitle())+"_"+paramval.str()).c_str());

    gpu_weights_plot[j]->Write();
    cpu_weights_plot[j]->Write();
    diff_weights_plot[j]->Write();

    gpu_weights_plot[j]->Reset();
    cpu_weights_plot[j]->Reset();
    diff_weights_plot[j]->Reset();

    gpu_weights_plot[j]->SetTitle(restoregpuname.c_str());
    cpu_weights_plot[j]->SetTitle(restorecpuname.c_str());
    diff_weights_plot[j]->SetTitle(restorediffname.c_str());
  }

  nReconf++;
  if (nReconf > 50) {
    std::cerr << "I've counted " << nReconf << " reconfigures" << std::endl;
    std::cerr << "Since DEBUG_DUMP saves ROOT files for every reconfigure I doubt you wan't more!" << std::endl;
    std::cerr << "Found " << badWeight << " bad weights" << std::endl;
    DebugFile->Close();
    throw;
  }
  #endif
}
#else
// *************************
// Prepare the weights before the reweight loop starts
void samplePDFND::PrepareWeights() {
// *************************
// Find the relevant spline segments for these parameter variations
  #if USE_SPLINE < USE_TF1
  FindSplineSegment();
  #endif
}
#endif

// *************************
void samplePDFND::GetKinVars(int Sample, KinematicTypes &TypeX, KinematicTypes &TypeY) {
// *************************
  //Lepton kinematic is default
  TypeX = kLeptonMomentum;
  TypeY = kLeptonCosTheta;

  if (samplepdfs->At(Sample) != NULL)
  {
    TypeX = kinvars[Sample][0];
    TypeY = kinvars[Sample][1];
  }
}

// **************************************************
//KS: Helper which initlaise PDF
void samplePDFND::InitialisePDF() {
// **************************************************

  // Take the overflow into account
  nSamples += 1;
  //KS: Learn how much bin we need for reweighting LLH calcaualtion etc.
  int GlobalNumberOfBin = 0;

  //For th2polys, the overflow bins are -1 to -9, depending on which direction your overflowing. For multithreading, we first add to a private array and then set this as the bin content at the end of the loop over events where array index = bin number. You can't have an array with negative indices so we'll use the final 9 bins as overflow. The last array element will equate to bin -1, the penultimate to bin -2 and so.
  for (__int__ i = 0; i < nSamples; i++)
  {
      if (samplepdfs->At(i) == NULL) continue;
      maxBins[i] += 10;
      GlobalNumberOfBin += maxBins[i];
  }
  std::cout<<"Using in total "<<GlobalNumberOfBin<<" bins for all samples"<<std::endl;

   //KS: Initialsie master arrays which keep number of events for each sample. Those are used in reweighting and LLH calcualtion but we still retain option to convert it back to TH2Poly
   // The master array for data
   samplePDF_data_array = new double*[nSamples]();
   // The master array for MC
   samplePDF_array = new double*[nSamples]();
   // The master array of weights^2, needed for Barlow Beeston
   samplePDF_w2_array = new double*[nSamples]();
   for (__int__ i = 0; i < nSamples; i++) {
     samplePDF_data_array[i] = new double[maxBins[i]]();
     samplePDF_array[i] = new double[maxBins[i]]();
     samplePDF_w2_array[i] = new double[maxBins[i]]();
     for (__int__ j = 0; j < maxBins[i]; ++j) {
       samplePDF_data_array[i][j] = 0.0;
       samplePDF_array[i][j] = 0.0;
       samplePDF_w2_array[i][j] = 0.0;
     }
   }

   //Update Data PDF now that we initalised
   UpdateDataPDF();
#ifdef MULTITHREAD
   // Allocate only if modepdf
   if (modepdf) {
   // The master mode array
   // Has an extra index to denote the mode
     samplePDF_mode_array = new double**[nSamples]();
     for (__int__ i = 0; i < nSamples; i++) {
       samplePDF_mode_array[i] = new double*[nModes+1]();
       for (__int__ j = 0; j < nModes+1; j++) {
         samplePDF_mode_array[i][j] = new double[maxBins[i]]();
         for (__int__ k = 0; k < maxBins[i]; ++k) {
           samplePDF_mode_array[i][j][k] = 0.0;
         }
       }
     }
   }
#endif
}

// **************************************************
//Helper which udpate data arrays from histogram
void samplePDFND::UpdateDataPDF() {
// **************************************************

    for (__int__ i = 0; i < nSamples; i++)
    {
      if (datapdfs->At(i) == NULL) continue;
      //Fill master array from histogram
      for (__int__ j = 0; j < maxBins[i]-__TH2PolyOverflowBins__; ++j)
      {
        samplePDF_data_array[i][j] = ((TH2Poly*)datapdfs->At(i))->GetBinContent(j);
      }
      for (__int__ j = -__TH2PolyOverflowBins__; j < 0; ++j)
      {
        samplePDF_data_array[i][maxBins[i]+j] = ((TH2Poly*)datapdfs->At(i))->GetBinContent(j);
      }
    }
}

// **************************************************
// Helper function to reset the data and MC histograms
void samplePDFND::ResetHistograms() {
// **************************************************

  // Loop over the samples and reset the and fill with zeros
  // Don't openMP this; no signficiant gain
  for (__int__ i = 0; i < nSamples; ++i)  
  {
    if (samplepdfs->At(i) == NULL) continue;
    for (__int__ j = 0; j < maxBins[i]; ++j)
    {
      //reset arrays
      samplePDF_array[i][j] = 0.0;
      if(firsttime) samplePDF_w2_array[i][j] = 0.0;
    }
    #if DEBUG > 0
    // Only reset the histograms that get refilled on event-by-event basis
    // i.e. don't need to reset POT scaled because they're fixed once we've looped once
    XsecOnly[i]->Reset("");
    NDCovOnly[i]->Reset("");
    AllOnly[i]->Reset("");
    #endif
  } // end loop over samples
  
  // If we want to plot according to mode
  if (modepdf) {
    for (__int__ i = 0; i < nSamples; ++i) {
      for (int j = 0; j < nModes+1; ++j) {
          // If 1D
          if (ndims[i] == 1) {
          ((TH1*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Reset();
          ((TH1*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Fill(0.0, 0.0);
          // If 2D
          } else if (ndims[i] == 2) {
          ((TH2Poly*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Reset("");
          ((TH2Poly*)((TObjArray*)samplemodepdfs->At(i))->At(j))->Fill(0.0, 0.0, 0.0);
          }
      } // end the mode loop
    } // end loop over samples
  } // end the if modepdf
} // end function


#if DEBUG > 0
// ***************************************************************************
// Dump the spline draw for a given event to an output file
// This is useful to check if the splines are bad
void samplePDFND::DumpSplines(double xsecw, double detw, int i) {
// ***************************************************************************

  std::cerr << "===============\nFOUND STRANGE EVENT IN SAMPLE\n===============" << std::endl;

  std::cerr << "Current XsecCov configuration:" << std::endl;
  XsecCov->printNominalCurrProp();

  PrintStructs(xsecw, detw, i);

  // Make a string describing the run
  std::string OutputName = std::string((*prod).Data());
  while (OutputName.find(" ") != std::string::npos) {
    OutputName.replace(OutputName.find(" "), 1, std::string("_"));
  }
  std::stringstream ss;
  ss << i;
  OutputName += "_event_" + ss.str() + ".root";

#if USE_SPLINE < USE_TF1
  // Make a TFile
  TFile *file = new TFile(OutputName.c_str(), "recreate");
  TCanvas *canv = new TCanvas("canv", "canv", 1080, 1080);
  file->cd();
  canv->cd();
  canv->Draw();
  // Loop over the spline parameters and write them to file
  for (int j = 0; j < nSplineParams; ++j) {
    if (xsecInfo[i].GetFunc(j) == NULL) continue;
    TSpline3 *Spline = NULL;
#if USE_SPLINE < USE_TSpline3
    Spline = xsecInfo[i].GetFunc(j)->ConstructTSpline3();
#else
    Spline = xsecInfo[i].GetFunc(j);
#endif
    Spline->Draw("LP*");
    canv->Write(Spline->GetTitle());
  }

  file->Write();
  file->Close();
  delete file;
  delete canv;
  std::cerr << "Wrote collected spline response to " << OutputName << std::endl;
#endif
}
#endif

#ifdef DEBUG
// ***************************************************************************
// Prints the relevant struct information and weights for one event
void samplePDFND::PrintStructs(double xsecw, double detw, const int i) {
// ***************************************************************************
  std::cout << "===============\nPRINTING EVENT BY EVENT WEIGHTS\n===============" << std::endl;
  std::cout << "EVENT " << i << std::endl;
  NDEve[i].Print();
  NDEve_Aux[i].Print();
  std::cout << std::setw(20) << "detw     = " << detw << std::endl;
  std::cout << std::setw(20) << "xsecw    = " << xsecw << std::endl;
  std::cout << std::setw(20) << "totalw (xsec*det*flux)  = " << xsecw*detw*NDEve[i].flux_w << std::endl;
}
#endif

#ifdef CUDA
// *********************************************
// Fill the GPU with splines
void samplePDFND::fillGPUSplines() {
// *********************************************

  // Can pass the spline segments to the GPU instead of the values
  // Make these here and only refill them for each loop, avoiding unnecessary new/delete on each reconfigure
  #if USE_SPLINE < USE_TF1
  segments = new int[nSplineParams]();
  #endif
  vals = new float[nSplineParams]();

  // Make a vector of all the vectors of TSpline3 pointers
  // Makes life slightly easier than passing 17 arguments and having different SMonolith constructors
  std::vector< std::vector<__SPLINE_TYPE__*> > MasterSpline;

  // The order we want in the master spline is [event][splinenumber]
  // Push back all the sweet splines for each event into one big vector
  for (unsigned int i = 0; i < nEvents; ++i) {
    // Make a vector of the pointers for each event
    std::vector<__SPLINE_TYPE__*> TempVector;
    for (int j = 0; j < nSplineParams; ++j) {
      TempVector.push_back(xsecInfo[i].GetFunc(j));
    } // End the for loop over splines
    MasterSpline.push_back(TempVector);
  }

  // Now pass the master spline to the GPU preparer
  splineMonolith = new SMonolith(MasterSpline);

  // Now need to reset the xsecInfo to match the MasterSpline
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int i = 0; i < nEvents; ++i) {
    // Make a vector of the pointers for each event
    for (int j = 0; j < nSplineParams; ++j) {
      xsecInfo[i].SetFunc(j, MasterSpline[i][j]);
    } // End the for loop over splines
  }

  // Make TH1D which logs differences in CPU and GPU weights for each reconfigure
  #ifdef DEBUG_DUMP
  std::string title = FitManager->getOutputFilename();
  title += "_GPUweights.root";
  DebugFile = new TFile(title.c_str(), "RECREATE");

  gpu_weights_plot = new TH1D*[nSplineParams];
  cpu_weights_plot = new TH1D*[nSplineParams];
  diff_weights_plot = new TH1D*[nSplineParams];
  for (int i = 0; i < nSplineParams; ++i) {
    gpu_weights_plot[i] = new TH1D((splineParsNames[i]+"_GPU").c_str(), (splineParsNames[i]+"_GPU").c_str(), nEvents, 0, nEvents);
    cpu_weights_plot[i] = new TH1D((splineParsNames[i]+"_CPU").c_str(), (splineParsNames[i]+"_CPU").c_str(), nEvents, 0, nEvents);
    diff_weights_plot[i] = new TH1D((splineParsNames[i]+"_diff").c_str(), (splineParsNames[i]+"_diff").c_str(), nEvents, 0, nEvents);

    gpu_weights_plot[i]->GetXaxis()->SetTitle("Event number");
    gpu_weights_plot[i]->GetYaxis()->SetTitle("Weight");
    cpu_weights_plot[i]->GetXaxis()->SetTitle(gpu_weights_plot[i]->GetXaxis()->GetTitle());
    cpu_weights_plot[i]->GetYaxis()->SetTitle(gpu_weights_plot[i]->GetYaxis()->GetTitle());
    diff_weights_plot[i]->GetXaxis()->SetTitle(gpu_weights_plot[i]->GetXaxis()->GetTitle());
    diff_weights_plot[i]->GetYaxis()->SetTitle("CPU - GPU weight");
  }
  // Clean up the memory if we aren't debugging
  #else
  CleanUpMemory();
  #endif
}

  #ifdef DEBUG_DUMP
// **************************************
// Compare the CPU and GPU weights
// Since we have both CPU and GPU splines in memory we might as well check they are the same!
void samplePDFND::CompareCPU_GPU_Splines(const int EventNumber) {
  // **************************************

  double CPU_weight = -1.0;
  double GPU_weight = -1.0;

  // Need to loop over the spline parameter
  for (int id = 0; id < nSplineParams; ++id) {
    // Get the splines
    __SPLINE_TYPE__ *spl = xsecInfo[EventNumber].GetFunc(id);
    // If this is a NULL pointer we don't have a spline so continue
    if (spl == NULL) continue;

    // Get the variation for this spline parameter ID
    const int GlobalIndex = splineParsIndex[id];
    const double var = XsecCov->calcReWeight(GlobalIndex);

    // Write the CPU and GPU weights
    GPU_weight = splineMonolith->cpu_weights[EventNumber*nSplineParams+id];
    CPU_weight = spl->Eval(var);
    double diff = (CPU_weight - GPU_weight);

    gpu_weights_plot[id]->SetBinContent(EventNumber, GPU_weight);
    cpu_weights_plot[id]->SetBinContent(EventNumber, CPU_weight);
    diff_weights_plot[id]->SetBinContent(EventNumber, diff);

    if (diff > 1.E-5) {
      badWeight++;
      std::cerr << "Found difference in splines greater than 1E-5!" << std::endl;

      std::cerr << "   Event no:      " << EventNumber << std::endl;
      std::cerr << "   Event mode:    " << Mode_ToString(xsecInfo[EventNumber].mode) << std::endl;
      std::cerr << "   Event species: " << xsecInfo[EventNumber].species << std::endl;

      std::cerr << "   Parameter:     " << splineParsNames[id] << std::endl;
      std::cerr << "   Variation:     " << var << std::endl;
      std::cerr << "   CPU_weight:    " << CPU_weight << std::endl;
      std::cerr << "   GPU_weight:    " << GPU_weight << std::endl;
      std::cerr << "   Diff:          " << diff << std::endl;

    #if USE_SPLINE < USE_TF1
      std::cerr << "   Segment:       " << SplineInfoArray[id].CurrSegment << std::endl;
      std::cerr << "   Segment spl:   " << spl->FindX(var) << std::endl;
    #else
      std::cerr << "   Coeffs:        " << std::endl;
      spl->Print();
    #endif
      std::cerr << "   " << badWeight << " weights bad so far" << std::endl;
      std::cerr << "   On " << nReconf << " reweights so far" << std::endl;
    }
  }
}
  #endif

// *************************************
// Helper function to clean up the splines which we've now got stores on the GPU
void samplePDFND::CleanUpMemory() {
// *************************************
  #ifndef DEBUG_DUMP
  // Delete all the TSpline3 memory allocated since this now lives on the GPU and is no longer needed on the GPU
  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (unsigned int i = 0; i < nEvents; ++i) {
    for (int j = 0; j < nSplineParams; ++j) {
      // Delete the allocated memory
      if (xsecInfo[i].GetFunc(j) != NULL) {
        delete xsecInfo[i].GetFunc(j);
      }
      // Set the memory to NULL
      __SPLINE_TYPE__* Temp = NULL;
      xsecInfo[i].SetFunc(j, Temp);
    }
  }
  #endif
}
#endif

