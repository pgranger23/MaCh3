#include "covarianceBase.h"

// ********************************************
covarianceBase::covarianceBase(const char *name, const char *file) : inputFile(std::string(file)), pca(false) {
// ********************************************

  MACH3LOG_INFO("Constructing instance of covarianceBase");
  init(name, file);
  FirstPCAdpar = -999;
  LastPCAdpar = -999;
}
// ********************************************
covarianceBase::covarianceBase(std::vector<std::string> YAMLFile, const char *name, double threshold, int FirstPCA, int LastPCA) : inputFile(YAMLFile[0].c_str()), matrixName(name), pca(true), eigen_threshold(threshold), FirstPCAdpar(FirstPCA), LastPCAdpar(LastPCA) {
// ********************************************

  MACH3LOG_INFO("Constructing instance of covarianceBase using ");
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    MACH3LOG_INFO("{}", YAMLFile[i]);
  }
  MACH3LOG_INFO("as an input");

  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("Principal component analysis but given the threshold for the principal components to be less than 0, or greater than (or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: ");
    MACH3LOG_INFO("Am instead calling the usual non-PCA constructor...");
    pca = false;
  }

  init(YAMLFile);
  // Call the innocent helper function
  if (pca) ConstructPCA();
}

// ********************************************
covarianceBase::covarianceBase(const char *name, const char *file, int seed) : inputFile(std::string(file)), pca(false) {
// ********************************************

  #ifdef MULTITHREAD
  if(seed != 0)
  {
    MACH3LOG_WARN("You have set seed to {}", seed);
    MACH3LOG_WARN("And you are running with MULTITHREAD");
    MACH3LOG_WARN("TRandom for each thread will have same seed");
    MACH3LOG_WARN("This is fine if this was your intention");
  }
  #endif

  init(name, file);
  FirstPCAdpar = -999;
  LastPCAdpar = -999;

}
// ********************************************
covarianceBase::covarianceBase(const char *name, const char *file, int seed, double threshold, int firstpcapar, int lastpcapar) : inputFile(std::string(file)), pca(true), eigen_threshold(threshold), FirstPCAdpar(firstpcapar), LastPCAdpar(lastpcapar) {
// ********************************************

  if (threshold < 0 || threshold >= 1) {
    MACH3LOG_INFO("NOTE: {} {}", name, file);
    MACH3LOG_INFO("Principal component analysis but given the threshold for the principal components to be less than 0, or greater than (or equal to) 1. This will not work");
    MACH3LOG_INFO("Please specify a number between 0 and 1");
    MACH3LOG_INFO("You specified: ");
    MACH3LOG_INFO("Am instead calling the usual non-PCA constructor...");
    pca = false;
  }
#ifdef MULTITHREAD
  if(seed != 0)
  {
    MACH3LOG_WARN("You have set seed to {}", seed);
    MACH3LOG_WARN("And you are running with MULTITHREAD");
    MACH3LOG_WARN("TRandom for each thread will have same seed");
    MACH3LOG_WARN("This is fine if this was your intention");
  }
#endif
  MACH3LOG_INFO("Constructing instance of covarianceBase");
  init(name, file);
  // Call the innocent helper function
  if (pca) ConstructPCA();
}

// ********************************************
//Destructor
covarianceBase::~covarianceBase(){
// ********************************************

  _fPreFitValue.clear();
  _fError.clear();
  _fCurrVal.clear();
  _fPropVal.clear();
  _fLowBound.clear();
  _fUpBound.clear();
  _fIndivStepScale.clear();
  _fFlatPrior.clear();

  delete[] randParams;
  delete[] corr_throw;

  if (covMatrix != NULL) delete covMatrix;
  if (invCovMatrix != NULL) delete invCovMatrix;
  if (throwMatrix_CholDecomp != NULL) delete throwMatrix_CholDecomp;

  for(int i = 0; i < _fNumPar; i++)
  {
    delete[] InvertCovMatrix[i];
    delete[] throwMatrixCholDecomp[i];
  }
  delete[] InvertCovMatrix;
  delete[] throwMatrixCholDecomp;
  
  const int nThreads = MaCh3Utils::GetNThreads();
  for (int iThread = 0;iThread < nThreads; iThread++)  delete random_number[iThread];
  delete[] random_number;
  if (throwMatrix != NULL) delete throwMatrix;
}

// ********************************************
void covarianceBase::ConstructPCA() {
// ********************************************

  // Check that covariance matrix exists
  if (covMatrix == NULL) {
    MACH3LOG_ERROR("Covariance matrix for {} has not yet been set", matrixName);
    MACH3LOG_ERROR("Can not construct PCA until it is set");
    throw;
  }

  //Check whether first and last pcadpar are set and if not just PCA everything
  if(FirstPCAdpar == -999 || LastPCAdpar == -999){
    if(FirstPCAdpar == -999 && LastPCAdpar == -999){
      FirstPCAdpar = 0;
      LastPCAdpar = covMatrix->GetNrows()-1;
    }
    else{
      MACH3LOG_ERROR("You must either leave FirstPCAdpar and LastPCAdpar at -999 or set them both to something");
      throw;
    }
  }
  if(FirstPCAdpar > covMatrix->GetNrows()-1 || LastPCAdpar>covMatrix->GetNrows()-1){
    MACH3LOG_ERROR("FirstPCAdpar and LastPCAdpar are higher than the number of parameters");
    MACH3LOG_ERROR("first: {} last: {}, params: {}", FirstPCAdpar, LastPCAdpar, covMatrix->GetNrows()-1);
    throw;
  }
  if(FirstPCAdpar < 0 || LastPCAdpar < 0){
    MACH3LOG_ERROR("FirstPCAdpar and LastPCAdpar are less than 0 but not default -999");
    MACH3LOG_ERROR("first: {} last: {}", FirstPCAdpar, LastPCAdpar);
    throw;
  }

  MACH3LOG_INFO("PCAing parameters {} through {} inclusive", FirstPCAdpar, LastPCAdpar);
  int numunpcadpars = covMatrix->GetNrows()-(LastPCAdpar-FirstPCAdpar+1);

  TMatrixDSym submat(covMatrix->GetSub(FirstPCAdpar,LastPCAdpar,FirstPCAdpar,LastPCAdpar));

  //CW: Calculate how many eigen values this threshold corresponds to
  TMatrixDSymEigen eigen(submat);
  eigen_values.ResizeTo(eigen.GetEigenValues());
  eigen_vectors.ResizeTo(eigen.GetEigenVectors());
  eigen_values = eigen.GetEigenValues();
  eigen_vectors = eigen.GetEigenVectors();
  double sum = 0;
  // Loop over eigen values and sum them up
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    sum += eigen_values(i);
  }
  nKeptPCApars = eigen_values.GetNrows();
  //CW: Now go through again and see how many eigen values correspond to threshold
  for (int i = 0; i < eigen_values.GetNrows(); ++i) {
    // Get the relative size of the eigen value
    double sig = eigen_values(i)/sum;
    // Check against the threshold
    if (sig < eigen_threshold) {
      nKeptPCApars = i;
      break;
    }
  }
  _fNumParPCA = numunpcadpars+nKeptPCApars;
  MACH3LOG_INFO("Threshold of {} on eigen values relative sum of eigen value ({}) generates {} eigen vectors, plus we have {} unpcad pars, for a total of {}", eigen_threshold, sum, nKeptPCApars, numunpcadpars, _fNumParPCA);

  //DB Create array of correct size so eigen_values can be used in CorrelateSteps
  eigen_values_master = std::vector<double>(_fNumParPCA, 1.0);
  for (int i = FirstPCAdpar; i < FirstPCAdpar+nKeptPCApars; ++i) {eigen_values_master[i] = eigen_values(i-FirstPCAdpar);}

  // Now construct the transfer matrices
  //These matrices will be as big as number of unPCAd pars plus number of eigenvalues kept
  TransferMat.ResizeTo(covMatrix->GetNrows(), _fNumParPCA);
  TransferMatT.ResizeTo(covMatrix->GetNrows(), _fNumParPCA);

  // Get a subset of the eigen vector matrix
  TMatrixD temp(eigen_vectors.GetSub(0, eigen_vectors.GetNrows()-1, 0, nKeptPCApars-1));
  
  //Make transfer matrix which is two blocks of identity with a block of the PCA transfer matrix in between
  TMatrixD temp2;
  temp2.ResizeTo(covMatrix->GetNrows(), _fNumParPCA);

  //First set the whole thing to 0
  for(int iRow = 0; iRow < covMatrix->GetNrows(); iRow++){
    for(int iCol = 0; iCol < _fNumParPCA; iCol++){
      temp2[iRow][iCol] = 0;
    }
  }
  //Set the first identity block
  if(FirstPCAdpar != 0){
    for(int iRow = 0; iRow < FirstPCAdpar; iRow++){
      temp2[iRow][iRow] = 1;
    }
  }

  //Set the transfer matrix block for the PCAd pars
  temp2.SetSub(FirstPCAdpar,FirstPCAdpar,temp);

  //Set the second identity block
  if(LastPCAdpar != covMatrix->GetNrows()-1){
    for(int iRow = 0;iRow < (covMatrix->GetNrows()-1)-LastPCAdpar; iRow++){
      temp2[LastPCAdpar+1+iRow][FirstPCAdpar+nKeptPCApars+iRow] = 1;
    }
  }
   
  TransferMat = temp2;
  // Copy the contents
  TransferMatT = TransferMat;
  // And then transpose
  TransferMatT.T();

  // Make a note that we have now done PCA
  pca = true;

  // Make the PCA parameter arrays
  fParCurr_PCA.ResizeTo(_fNumParPCA);
  fParProp_PCA.ResizeTo(_fNumParPCA);
  _fPreFitValue_PCA.resize(_fNumParPCA);

  //KS: make easy map so we could easily find un-decomposed parameters
  isDecomposed_PCA.resize(_fNumParPCA);
  fParSigma_PCA.resize(_fNumParPCA);
  for (int i = 0; i < _fNumParPCA; ++i)
  {
    fParSigma_PCA[i] = 1;
    isDecomposed_PCA[i] = -1;
  }
  for (int i = 0; i < FirstPCAdpar; ++i) isDecomposed_PCA[i] = i;
  
  for (int i = FirstPCAdpar+nKeptPCApars+1; i < _fNumParPCA; ++i) isDecomposed_PCA[i] = i+(size-_fNumParPCA);

  #ifdef DEBUG_PCA
  //KS: Let's dump all useful matrices to properly validate PCA
  DebugPCA(sum, temp, submat);
  #endif
}

void covarianceBase::init(const char *name, const char *file)
{
  // Set the covariance matrix from input ROOT file (e.g. flux, ND280, NIWG)
  TFile *infile = new TFile(file, "READ");
  if (infile->IsZombie()) {
    MACH3LOG_ERROR("Could not open input covariance ROOT file {} !!!", file);
    MACH3LOG_ERROR("Was about to retrieve matrix with name {}", name);
    throw;
  }

  // Should put in a 
  TMatrixDSym *CovMat = (TMatrixDSym*)(infile->Get(name));
  if (CovMat == NULL) {
    std::cerr << "Could not find covariance matrix name " << name << " in file " << file << std::endl;
    std::cerr << "Are you really sure " << name << " exists in the file?" << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__  << std::endl;
    throw;
  } 

  PrintLength = 35;

  const int nThreads = MaCh3Utils::GetNThreads();
  //KS: set Random numbers for each thread so each thread has different seed
  //or for one thread if without MULTITHREAD
  random_number = new TRandom3*[nThreads]();
  for (int iThread = 0; iThread < nThreads; iThread++) {
    random_number[iThread] = new TRandom3(0);
  }

  // Not using adaptive by default
  use_adaptive = false;
  total_steps = 0;
  lower_adapt = -999;
  upper_adapt = -999;

  // Set the covariance matrix
  size = CovMat->GetNrows();
  _fNumPar = size;
    
  InvertCovMatrix = new double*[_fNumPar]();
  throwMatrixCholDecomp = new double*[_fNumPar]();
  // Set the defaults to true
  for(int i = 0; i < _fNumPar; i++)
  {
    InvertCovMatrix[i] = new double[_fNumPar]();
    throwMatrixCholDecomp[i] = new double[_fNumPar]();
    for (int j = 0; j < _fNumPar; j++)
    {
      InvertCovMatrix[i][j] = 0.;
      throwMatrixCholDecomp[i][j] = 0.;
    }
  }

  setName(name);
  MakePosDef(CovMat);
  setCovMatrix(CovMat);

  if (_fNumPar <= 0) {
    MACH3LOG_CRITICAL("Covariance matrix {} has {} entries!", getName(), _fNumPar);
    throw;
  }
  _fNumParPCA = _fNumPar;

  ReserveMemory(_fNumPar);

  infile->Close();

  MACH3LOG_INFO("Created covariance matrix named: {}", getName());
  MACH3LOG_INFO("from file: {}", file);

  delete infile;
}

// ETA
// An init function for the YAML constructor
// All you really need from the YAML file is the number of Systematics
// Then get all the info from the YAML file in the covarianceXsec::ParseYAML function
void covarianceBase::init(std::vector<std::string> YAMLFile)
{
  _fYAMLDoc["Systematics"] = YAML::Node(YAML::NodeType::Sequence);
  for(unsigned int i = 0; i < YAMLFile.size(); i++)
  {
    YAML::Node YAMLDocTemp = YAML::LoadFile(YAMLFile[i]);
    for (const auto& item : YAMLDocTemp["Systematics"]) {
      _fYAMLDoc["Systematics"].push_back(item);
    }
  }

  const int nThreads = MaCh3Utils::GetNThreads();
  //KS: set Random numbers for each thread so each thread has different seed
  //or for one thread if without MULTITHREAD
  random_number = new TRandom3*[nThreads]();
  for (int iThread = 0; iThread < nThreads; iThread++) {
    random_number[iThread] = new TRandom3(0);
  }

  PrintLength = 35;

  // Not using adaptive by default
  use_adaptive = false;
  total_steps = 0;
  lower_adapt = -999;
  upper_adapt = -999;

  // Set the covariance matrix
  _fNumPar = _fYAMLDoc["Systematics"].size();
  size = _fNumPar;
    
  InvertCovMatrix = new double*[_fNumPar]();
  throwMatrixCholDecomp = new double*[_fNumPar]();
  // Set the defaults to true
  for(int i = 0; i < _fNumPar; i++)
  {
    InvertCovMatrix[i] = new double[_fNumPar]();
    throwMatrixCholDecomp[i] = new double[_fNumPar]();
    for (int j = 0; j < _fNumPar; j++)
    {
      InvertCovMatrix[i][j] = 0.;
      throwMatrixCholDecomp[i][j] = 0.;
    }
  }

  ReserveMemory(_fNumPar);

  TMatrixDSym* _fCovMatrix = new TMatrixDSym(_fNumPar);
  //_fDetString = std::vector<std::string>(_fNumPar);

  int i = 0;
  std::vector<std::map<std::string,double>> Correlations(_fNumPar);
  std::map<std::string, int> CorrNamesMap;

  //ETA - read in the systematics. Would be good to add in some checks to make sure
  //that there are the correct number of entries i.e. are the _fNumPar for Names,
  //PreFitValues etc etc.
  for (auto const &param : _fYAMLDoc["Systematics"])
  {
    _fFancyNames[i] = (param["Systematic"]["Names"]["FancyName"].as<std::string>());
    _fPreFitValue[i] = (param["Systematic"]["ParameterValues"]["PreFitValue"].as<double>());
    _fGenerated[i] = (param["Systematic"]["ParameterValues"]["Generated"].as<double>());
    _fIndivStepScale[i] = (param["Systematic"]["StepScale"]["MCMC"].as<double>());
    _fError[i] = (param["Systematic"]["Error"].as<double>());

    //ETA - a bit of a fudge but works
    std::vector<double> TempBoundsVec = param["Systematic"]["ParameterBounds"].as<std::vector<double>>();
    _fLowBound[i] = TempBoundsVec[0];
    _fUpBound[i] = TempBoundsVec[1];

    //ETA - now for parameters which are optional and have default values
    _fFlatPrior[i] = GetFromManager<bool>(param["Systematic"]["FlatPrior"], false);

    //Fill the map to get the correlations later as well
    CorrNamesMap[param["Systematic"]["Names"]["FancyName"].as<std::string>()]=i;

    //Also loop through the correlations
    if(param["Systematic"]["Correlations"]) {
      for(unsigned int Corr_i = 0; Corr_i < param["Systematic"]["Correlations"].size(); ++Corr_i){
        for (YAML::const_iterator it=param["Systematic"]["Correlations"][Corr_i].begin();it!=param["Systematic"]["Correlations"][Corr_i].end();++it) {
          Correlations[i][it->first.as<std::string>()] = it->second.as<double>();
        }
      }
    }
    i++;
  }

  //ETA
  //Now that we've been through all systematic let's fill the covmatrix
  //This makes the root TCov from YAML
  for(int j = 0; j < _fNumPar;j++) {
    (*_fCovMatrix)(j, j)=_fError[j]*_fError[j];
    //Get the map of parameter name to correlation fomr the Correlations object
    for (auto const& [key, val] : Correlations[j]) {
      int index = -1;

      //If you found the parameter name then get the index
      if (CorrNamesMap.find(key) != CorrNamesMap.end()) {
        index = CorrNamesMap[key];
      }
      else {
        MACH3LOG_ERROR("Parameter {} not in list! Check your spelling?", key);
        exit(5);
      }

      //
      double Corr1 = val;
      double Corr2 = 0;
      if(Correlations[index].find(_fFancyNames[j]) != Correlations[index].end()) {
        Corr2 = Correlations[index][_fFancyNames[j]];
        //Do they agree to better than float precision?
        if(std::abs(Corr2 - Corr1) > FLT_EPSILON) {
          MACH3LOG_ERROR("Correlations are not equal between {} and {}", _fFancyNames[j], key);
          MACH3LOG_ERROR("Got : {} and {}", Corr2, Corr1);
          exit(5);
        }
      } else {
        MACH3LOG_ERROR("Correlation does not appear reciprocally between {} and {}", _fFancyNames[j], key);
        exit(5);
      }
      (*_fCovMatrix)(j, index)= (*_fCovMatrix)(index, j) = Corr1*_fError[j]*_fError[index];
    }
  }

  //Now make positive definite
  MakePosDef(_fCovMatrix);
  setCovMatrix(_fCovMatrix);

  if (_fNumPar <= 0) {
    MACH3LOG_ERROR("Covariance object has {} systematics!", _fNumPar);
    throw;
  }
  _fNumParPCA = _fNumPar;

  MACH3LOG_INFO("Created covariance matrix from files: ");
  for(const auto &file : YAMLFile){
    MACH3LOG_INFO("{} ", file);
  }
  MACH3LOG_INFO("----------------");
  MACH3LOG_INFO("Found {} systematics parameters in total", _fNumPar);
  MACH3LOG_INFO("----------------");

  return;
}

void covarianceBase::init(TMatrixDSym* covMat) {

  size = covMat->GetNrows();
  _fNumPar = size;
  InvertCovMatrix = new double*[_fNumPar]();
  throwMatrixCholDecomp = new double*[_fNumPar]();
  // Set the defaults to true
  for(int i = 0; i < _fNumPar; i++)
  {
    InvertCovMatrix[i] = new double[_fNumPar]();
    throwMatrixCholDecomp[i] = new double[_fNumPar]();
    for (int j = 0; j < _fNumPar; j++)
    {
      InvertCovMatrix[i][j] = 0.;
      throwMatrixCholDecomp[i][j] = 0.;
    }
  }
  
  setCovMatrix(covMat);

  ReserveMemory(_fNumPar);

  MACH3LOG_INFO("Created covariance matrix named: {}", getName());
}

// Set the covariance matrix for this class
void covarianceBase::setCovMatrix(TMatrixDSym *cov) {
  if (cov == NULL) {
    MACH3LOG_ERROR("Could not find covariance matrix you provided to setCovMatrix");
    MACH3LOG_ERROR("{}:{}", __FILE__, __LINE__);
    throw;
  }
  covMatrix = cov;
  invCovMatrix = (TMatrixDSym*)cov->Clone();
  invCovMatrix->Invert();
  //KS: ROOT has bad memory management, using standard double means we can decrease most operation by factor 2 simply due to cache hits
  for (int i = 0; i < _fNumPar; i++)
  {
    for (int j = 0; j < _fNumPar; ++j)
    {
      InvertCovMatrix[i][j] = (*invCovMatrix)(i,j);
    }
  }

  setThrowMatrix(cov);
}

void covarianceBase::ReserveMemory(const int SizeVec) {

  _fNames = std::vector<std::string>(SizeVec);
  _fFancyNames = std::vector<std::string>(SizeVec);
  _fGenerated = std::vector<double>(SizeVec);
  _fPreFitValue = std::vector<double>(SizeVec);
  _fError = std::vector<double>(SizeVec);
  _fCurrVal = std::vector<double>(SizeVec);
  _fPropVal = std::vector<double>(SizeVec);
  _fLowBound = std::vector<double>(SizeVec);
  _fUpBound = std::vector<double>(SizeVec);
  _fFlatPrior = std::vector<bool>(SizeVec);
  _fIndivStepScale = std::vector<double>(SizeVec);

  corr_throw = new double[SizeVec]();
  // set random parameter vector (for correlated steps)
  randParams = new double[SizeVec];

  // Set the defaults to true
  for(int i = 0; i < SizeVec; i++) {
    _fGenerated.at(i) = 1.;
    _fPreFitValue.at(i) = 1.;
    _fError.at(i) = 0.;
    _fCurrVal.at(i) = 0.;
    _fPropVal.at(i) = 0.;
    _fLowBound.at(i) = -999.99;
    _fUpBound.at(i) = 999.99;
    _fFlatPrior.at(i) = false;
    _fIndivStepScale.at(i) = 1.;
    corr_throw[i] = 0.0;
  }

  _fGlobalStepScale = 1.0;
}

// Set all the covariance matrix parameters to a user-defined value
// Might want to split this
void covarianceBase::setPar(int i , double val) {

  MACH3LOG_INFO("Over-riding {}: ", GetParName(i));
  MACH3LOG_INFO("_fPropVal ({}), _fCurrVal ({}), _fPreFitValue ({}) to ({})", _fPropVal[i], _fCurrVal[i], _fPreFitValue[i], val);

  _fPropVal[i] = val;
  _fCurrVal[i] = val;
  _fPreFitValue[i] = val;

  // Transfer the parameter values to the PCA basis
  if (pca) TransferToPCA();
}

// Transfer a parameter variation in the parameter basis to the eigen basis
void covarianceBase::TransferToPCA() {
  if (!pca) {
    MACH3LOG_ERROR("Can not transfer to PCA if PCA isn't enabled");
    MACH3LOG_ERROR("{}:{}", __FILE__, __LINE__);
    throw;
  }
  // Make the temporary vectors
  TVectorD fParCurr_vec(_fNumPar);
  TVectorD fParProp_vec(_fNumPar);
  for (int i = 0; i < _fNumPar; ++i) {
    fParCurr_vec(i) = _fCurrVal[i];
    fParProp_vec(i) = _fPropVal[i];
  }

  fParCurr_PCA = TransferMatT*fParCurr_vec;
  fParProp_PCA = TransferMatT*fParProp_vec;
}

// Transfer a parameter variation in the eigen basis to the parameter basis
void covarianceBase::TransferToParam() {
  if (!pca) {
    MACH3LOG_ERROR("Can not transfer to PCA if PCA isn't enabled");
    MACH3LOG_ERROR("{}:{}", __FILE__, __LINE__);
    throw;
  }

  // Make the temporary vectors
  TVectorD fParProp_vec = TransferMat*fParProp_PCA;
  TVectorD fParCurr_vec = TransferMat*fParCurr_PCA;
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < _fNumPar; ++i) {
    _fPropVal[i] = fParProp_vec(i);
    _fCurrVal[i] = fParCurr_vec(i);
  }
}

const std::vector<double> covarianceBase::getProposed() const {
  std::vector<double> props(_fNumPar);
  for (int i = 0; i < _fNumPar; ++i) props[i] = _fPropVal[i];
  return props;
}

// Throw nominal values
void covarianceBase::throwNominal(bool nomValues, int seed) {

  TVectorD* vec = new TVectorD(_fNumPar);
  for (int i = 0; i < _fNumPar; i++) {
    (*vec)(i) = 1.0;
  }

  ThrowParms* nom_throws = new ThrowParms(*vec, (*covMatrix));
  nom_throws->SetSeed(seed);
  std::vector<double> nominal = getNominalArray();
  nominal.clear();
  nominal.resize(_fNumPar);

  // If we want to put the nominals somewhere else than user specified
  // Don't fully understand this though: won't we have to reweight the MC somehow?
  // nominal[i] is used in GetLikelihood() as the penalty term, so we're essentially setting a random parameter penalty term?
  if (!nomValues)
  {
    bool throw_again = true;

    while(throw_again == true)
    {
      throw_again = false;
      MACH3LOG_INFO("Setting {} nominal values to random throws.", getName());
      nom_throws->ThrowSet(nominal);

      for (int i = 0; i < _fNumPar; i++)
      {
      // if parameter is fixed, dont throw
        if (_fError[i] < 0) {
          nominal[i] = 1.0;
          continue;
        }

        if (nominal[i] < 0) {
          nominal[i] = 0.0;
          throw_again = true;
        }
      }
    }
  } else {

    // If we want nominal values, set all entries to 1 (defined as nominal in MaCh3)
    for (int i = 0; i < int(nominal.size()); i++) {
      nominal[i] = 1.0;
    }
  }

  delete nom_throws;
  delete vec;
}

// *************************************
// Throw the parameters according to the covariance matrix
// This shouldn't be used in MCMC code ase it can break Detailed Balance;
void covarianceBase::throwParameters() {
// *************************************

  // First draw new randParams
  randomize();

  if (!pca) {
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumPar; ++i) {
      // Check if parameter is fixed first: if so don't randomly throw
      if (isParameterFixed(i)) continue;

      _fPropVal[i] = _fPreFitValue[i] + corr_throw[i];
      int throws = 0;
      // Try again if we the initial parameter proposal falls outside of the range of the parameter
      while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
#ifdef MULTITHREAD
        randParams[i] = random_number[omp_get_thread_num()]->Gaus(0, 1);
#else
        randParams[i] = random_number[0]->Gaus(0,1);
#endif
        const double corr_throw_single = MatrixVectorMultiSingle(throwMatrixCholDecomp, randParams, _fNumPar, i);
        _fPropVal[i] = _fPreFitValue[i] + corr_throw_single;
        if (throws > 10000) 
        {
          //KS: Since we are multithreading there is danger that those messages
          //will be all over the place, small price to pay for faster code
          MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws, i);
          MACH3LOG_WARN("Matrix: {}", matrixName);
          MACH3LOG_WARN("Param: {}", _fNames[i]);
          MACH3LOG_WARN("Setting _fPropVal:  {} to {}", _fPropVal[i], _fPreFitValue[i]);
          MACH3LOG_WARN("I live at {}:{}", __FILE__, __LINE__);
          _fPropVal[i] = _fPreFitValue[i];
          //throw;
        }
        throws++;
      }
      _fCurrVal[i] = _fPropVal[i];
    }
  }
  else
  {
    MACH3LOG_CRITICAL("Hold on, you are trying to run Prior Predicitve Code with PCA, which is wrong");
    MACH3LOG_CRITICAL("Sorry I have to kill you, I mean your job");
    throw;
  }
}

// *************************************
// Throw each parameter within their 1 sigma range
// Used to start the chain in different states
void covarianceBase::RandomConfiguration() {
  // *************************************
  // Have the 1 sigma for each parameter in each covariance class, sweet!
  // Don't want to change the nominal array because that's what determines our likelihood
  // Want to change the fParProp, fParCurr, fParInit
  // fParInit and the others will already be set
  for (int i = 0; i < _fNumPar; ++i) {
    // Check if parameter is fixed first: if so don't randomly throw
    if (isParameterFixed(i)) continue;
    // Check that the sigma range is larger than the parameter range
    // If not, throw in the valid parameter range instead
    const double paramrange = _fUpBound[i] - _fLowBound[i];
    const double sigma = sqrt((*covMatrix)(i,i));
    double throwrange = sigma;
    if (paramrange < sigma) throwrange = paramrange;

    _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1)*throwrange;
    // Try again if we the initial parameter proposal falls outside of the range of the parameter
    // Really only relevant for the xsec parameters; the flux and ND280 have -999 and 999 set to the limits!
    int throws = 0;
    while (_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]) {
      if (throws > 1000) {
        MACH3LOG_WARN("Tried {} times to throw parameter {} but failed", throws, i);
        MACH3LOG_WARN("Matrix: {}", matrixName);
        MACH3LOG_WARN("Param: {}", _fNames[i]);
        MACH3LOG_WARN("I live at {}:{}", __FILE__, __LINE__);
        throw;
      }
      _fPropVal[i] = _fPreFitValue[i] + random_number[0]->Gaus(0, 1)*throwrange;
      throws++;
    }
    MACH3LOG_INFO("Setting current step in {} param {} = {} from {}", matrixName, i, _fPropVal[i], _fCurrVal[i]);
    _fCurrVal[i] = _fPropVal[i];
  }
  if (pca) TransferToPCA();

}

// *************************************
// Set a single parameter
void covarianceBase::setSingleParameter(const int parNo, const double parVal) {
  // *************************************
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  std::cout << "Setting " << GetParName(parNo) << "(parameter " << parNo << ") to " << parVal << std::endl;

  if (pca) TransferToPCA();
}

void covarianceBase::setParCurrProp(const int parNo, const double parVal) {
  _fPropVal[parNo] = parVal;
  _fCurrVal[parNo] = parVal;
  std::cout << "Setting " << GetParName(parNo) << "(parameter " << parNo << ") to " << parVal << std::endl;
  if (pca) TransferToPCA();
}

// Propose a step for the set of systematics parameters this covariance class holds
void covarianceBase::proposeStep() {
  // Make the random numbers for the step proposal
  randomize();
  CorrelateSteps();
  if(use_adaptive && total_steps < upper_adapt) updateAdaptiveCovariance();
}

// ************************************************
// "Randomize" the parameters in the covariance class for the proposed step
// Used the proposal kernel and the current parameter value to set proposed step
// Also get a new random number for the randParams
void covarianceBase::randomize() {
// ************************************************
  if (!pca) {
//KS: By multithreading here we gain at least factor 2 with 8 threads with ND only fit      
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumPar; ++i) {
      // If parameter isn't fixed
      if (_fError[i] > 0.0) {
#ifdef MULTITHREAD
        randParams[i] = random_number[omp_get_thread_num()]->Gaus(0, 1);
#else
        randParams[i] = random_number[0]->Gaus(0, 1);
#endif 
        // If parameter IS fixed
      } else {
        randParams[i] = 0.0;
      }
    } // end for
  // If we're in the PCA basis we instead throw parameters there (only _fNumParPCA parameter)
  } else {
    // Scale the random parameters by the sqrt of eigen values for the throw
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumPar; ++i)
    {
      if (fParSigma_PCA[i] > 0. && i < _fNumParPCA)
      {
#ifdef MULTITHREAD        
        randParams[i] = random_number[omp_get_thread_num()]->Gaus(0,1);
#else        
        randParams[i] = random_number[0]->Gaus(0,1);
#endif
      } else { // If parameter IS fixed or out od bounds
        randParams[i] = 0.0;
      }
    }
  }
}


// ************************************************
// Correlate the steps by setting the proposed step of a parameter to its current value + some correlated throw
void covarianceBase::CorrelateSteps() {
// ************************************************
  //KS: Using custom function compared to ROOT one with 8 threads we have almost factor 2 performance increase, by replacing TMatrix with just double we increase it even more
  MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);

  // If not doing PCA
  if (!pca) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumPar; ++i) {
      if (_fError[i] > 0.) {
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i]*_fGlobalStepScale*_fIndivStepScale[i];
      }
    }
    // If doing PCA throw uncorrelated in PCA basis (orthogonal basis by definition)
  } else { 
    // Throw around the current step
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumParPCA; ++i)
    {
      if (fParSigma_PCA[i] > 0.) 
      {
        double IndStepScale = 1.;
        //KS: If undecomposed parameter apply individual step scale and Cholesky for better acceptance rate
        if(isDecomposed_PCA[i] >= 0)
        {
          IndStepScale *= _fIndivStepScale[isDecomposed_PCA[i]];
          IndStepScale *= corr_throw[isDecomposed_PCA[i]];
        }
        //If decomposed apply only random number
        else
        {
         IndStepScale *= randParams[i];
         //KS: All PCA-ed parameters have the same step scale
         IndStepScale *= _fIndivStepScale[FirstPCAdpar];
        }
        fParProp_PCA(i) = fParCurr_PCA(i)+_fGlobalStepScale*IndStepScale*eigen_values_master[i];
      }
    }
    // Then update the parameter basis
    TransferToParam();
  }
}

// Update so that current step becomes the previously proposed step
void covarianceBase::acceptStep() {
  if (!pca) {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumPar; ++i) {
      // Update state so that current state is proposed state
      _fCurrVal[i] = _fPropVal[i];
    }
  } else {
  // Update the book-keeping for the output
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
    for (int i = 0; i < _fNumParPCA; ++i) {
      fParCurr_PCA(i) = fParProp_PCA(i);
    }
    // Then update the parameter basis
    TransferToParam();
  }
}

// Throw the proposed parameter by mag sigma
// Should really just have the user specify this throw by having argument double
void covarianceBase::throwParProp(const double mag) {
  randomize();
  if (!pca) {
    // Make the correlated throw
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);
    // Number of sigmas we throw
    for (int i = 0; i < _fNumPar; i++) {
      if (_fError[i] > 0.)
        _fPropVal[i] = _fCurrVal[i] + corr_throw[i]*mag;
    }
  } else {
    for (int i = 0; i < _fNumPar; i++) {
      fParProp_PCA(i) = fParCurr_PCA(i)+mag*randParams[i];
    }
    // And update the fParCurr in the basis
    // Then update the fParProp in the parameter basis using the transfer matrix, so likelihood is evaluated correctly
    TVectorD proposed = TransferMat*fParProp_PCA;
    for (int i = 0; i < _fNumPar; ++i) {
      if (fParSigma_PCA[i] > 0.) {
        _fPropVal[i] = proposed(i);
      }
    }
  }
}

// Helper function to throw the current parameter by mag sigmas
// Can study bias in MCMC with this; put different starting parameters
void covarianceBase::throwParCurr(const double mag)
{
  randomize();
  if (!pca) {
    // Get the correlated throw vector
    MatrixVectorMulti(corr_throw, throwMatrixCholDecomp, randParams, _fNumPar);
    // The number of sigmas to throw
    // Should probably have this as a default parameter input to the function instead
    for (int i = 0; i < _fNumPar; i++) {
      if (_fError[i] > 0.){
        _fCurrVal[i] = corr_throw[i]*mag;
      }
    }
  } else {
    for (int i = 0; i < _fNumPar; i++) {
      fParProp_PCA(i) = mag*randParams[i];
    }
    // And update the fParCurr in the basis
    // Then update the fParProp in the parameter basis using the transfer matrix, so likelihood is evaluated correctly
    TVectorD current = TransferMat*fParCurr_PCA;
    for (int i = 0; i < _fNumPar; ++i) {
      if (fParSigma_PCA[i] > 0.) {
        _fCurrVal[i] = current(i);
      }
    }
  }
}

// Function to print the nominal values
void covarianceBase::printNominal() {
  MACH3LOG_INFO("Prior values for {} covarianceBase:", getName());
  for (int i = 0; i < _fNumPar; i++) {
    MACH3LOG_INFO("    {}   {} ", GetParFancyName(i), getParInit(i));
  }
}

// Function to print the nominal, current and proposed values
void covarianceBase::printNominalCurrProp() {

  MACH3LOG_INFO("Printing parameters for {}", getName());
  // Dump out the PCA parameters too
  if (pca) {
    MACH3LOG_INFO("PCA:");
    for (int i = 0; i < _fNumParPCA; ++i) {
      std::cout << std::setw(PrintLength) << std::left << "PCA " << i << " Current: " << fParCurr_PCA(i) << " Proposed: " << fParProp_PCA(i) << std::endl;
    }
  }
  MACH3LOG_INFO("{:<30} {:<10} {:<10} {:<10}", "Name", "Prior", "Current", "Proposed");
  for (int i = 0; i < _fNumPar; ++i) {
    MACH3LOG_INFO("{:<30} {:<10.2f} {:<10.2f} {:<10.2f}", GetParFancyName(i), _fPreFitValue[i], _fCurrVal[i], _fPropVal[i]);
  }
  //KS: "\n" is faster performance wise, keep std::endl at the end to flush just in case, also looks pretty
  //std::cout << std::endl;
}

// ********************************************
// Get the likelihood in the case where we want to include priors on the parameters
// fParEvalLikelihood stores if we want to evaluate the likelihood for the given parameter
//                    true = evaluate likelihood (so run with a prior)
//                    false = don't evaluate likelihood (so run without a prior)
double covarianceBase::CalcLikelihood() {
// ********************************************

  double logL = 0.0;
  //TStopwatch clock;
  ///clock.Start();
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:logL)
  #endif
  for(int i = 0; i < _fNumPar; ++i){
    for (int j = 0; j <= i; ++j) {
      if (!_fFlatPrior[i] && !_fFlatPrior[j]) {
        //KS: Since matrix is symmetric we can calculate non diagonal elements only once and multiply by 2, can bring up to factor speed decrease.
        short int scale = (i != j) ? 2 : 1;
        logL += scale * 0.5*(_fPropVal[i] - _fPreFitValue[i])*(_fPropVal[j] - _fPreFitValue[j])*InvertCovMatrix[i][j];
      }
    }
  }
  //clock.Stop();
  //std::cout << __FILE__ << "::GetLikelihood took " << clock.RealTime() << "s" << std::endl;

  return logL;
}
// ********************************************
int covarianceBase::CheckBounds() {
// ********************************************

  int NOutside = 0;
  #ifdef MULTITHREAD
  #pragma omp parallel for reduction(+:NOutside)
  #endif
  for (int i = 0; i < _fNumPar; ++i){
    if(_fPropVal[i] > _fUpBound[i] || _fPropVal[i] < _fLowBound[i]){
      NOutside++;
    }
  }
  return NOutside;
}
// ********************************************
double covarianceBase::GetLikelihood() {
// ********************************************

  // Checkbounds and calclikelihood are virtual
  // Default behaviour is to reject negative values + do std llh calculation
  const int NOutside = CheckBounds();
  
  if(NOutside > 0)
    return NOutside*_LARGE_LOGL_;

  return CalcLikelihood();
}


void covarianceBase::printPars() {

  MACH3LOG_INFO("Number of pars: {}", _fNumPar);
  MACH3LOG_INFO("Current {} parameters:", matrixName);
  for(int i = 0; i < _fNumPar; i++) {
    std::cout << std::fixed << std::setprecision(5) << _fNames[i].c_str() << " current: \t" << _fCurrVal[i] << "   \tproposed: \t" << _fPropVal[i] << std::endl;
  }

  return;
}

// Sets the proposed parameters to the nominal values
void covarianceBase::setParameters(std::vector<double> pars) {

  // If empty, set the proposed to nominal
  if (pars.empty()) {
    // For xsec this means setting to the prior (because nominal is the prior)
    for (int i = 0; i < _fNumPar; i++) {
      _fPropVal[i] = _fPreFitValue[i];
    }
    // If not empty, set the parameters to the specified
  } else {

	if (pars.size() != size_t(_fNumPar)) {
      MACH3LOG_ERROR("Warning: parameter arrays of incompatible size! Not changing parameters! {} has size {} but was expecting {}", matrixName, pars.size(), _fNumPar);
      throw;
    }

    unsigned int parsSize = pars.size();
    for (unsigned int i = 0; i < parsSize; i++) {
	  //Make sure that you are actually passing a number to set the parameter to
	  if(std::isnan(pars[i])) {
		std::cerr << "Error: trying to set parameter value to a nan for parameter " << GetParName(i) << " in matrix " << matrixName << ". This will not go well!" << std::endl;
		throw;
	  } else {
		_fPropVal[i] = pars[i];
	  }
    }
  }

  // And if pca make the transfer
  if (pca) {
    TransferToPCA();
    TransferToParam();
  }

  return;
}

// ********************************************
void covarianceBase::SetBranches(TTree &tree, bool SaveProposal) {
// ********************************************

  // loop over parameters and set a branch
  for (int i = 0; i < _fNumPar; ++i) {
    tree.Branch(_fNames[i].c_str(), &_fCurrVal[i], Form("%s/D", _fNames[i].c_str()));
  }
  // When running PCA, also save PCA parameters
  if (pca) {
    for (int i = 0; i < _fNumParPCA; ++i) {
      tree.Branch(Form("%s_PCA", _fNames[i].c_str()), (double*)&(fParCurr_PCA.GetMatrixArray()[i]), Form("%s_PCA/D", _fNames[i].c_str()));
    }
  }
  if(SaveProposal)
  {
    // loop over parameters and set a branch
    for (int i = 0; i < _fNumPar; ++i) {
      tree.Branch(Form("%s_Prop", _fNames[i].c_str()), &_fPropVal[i], Form("%s_Prop/D", _fNames[i].c_str()));
    }
    // When running PCA, also save PCA parameters
    if (pca) {
      for (int i = 0; i < _fNumParPCA; ++i) {
        tree.Branch(Form("%s_PCA_Prop", _fNames[i].c_str()), (double*)&(fParProp_PCA.GetMatrixArray()[i]), Form("%s_PCA_Prop/D", _fNames[i].c_str()));
      }
    }
  }
}
// ********************************************
void covarianceBase::setStepScale(const double scale) {
// ********************************************
  if(scale == 0)
  {
    MACH3LOG_ERROR("You are trying so set StepScale to 0 this will not work");
    throw;
  }
  MACH3LOG_INFO("{} setStepScale() = {}", getName(), scale);
  _fGlobalStepScale = scale;
}

// ********************************************
void covarianceBase::toggleFixAllParameters() {
// ********************************************
  // fix or unfix all parameters by multiplying by -1
  if(!pca)
  {
    for (int i = 0; i < _fNumPar; i++) _fError[i] *= -1.0;
  } else{
     for (int i = 0; i < _fNumParPCA; i++) fParSigma_PCA[i] *= -1.0;
  }
  return;
}

// ********************************************
void covarianceBase::toggleFixParameter(const int i) {
// ********************************************
  if(!pca) {
	if (i > _fNumPar) {
      MACH3LOG_ERROR("Can't toggleFixParameter for parameter {} because size of covariance ={}", i, _fNumPar);
      MACH3LOG_ERROR("Fix this in your config file please!");
	  throw;
	} else {
	  _fError[i] *= -1.0;
      MACH3LOG_INFO("Setting {}(parameter {}) to fixed at {}", GetParName(i), i, _fCurrVal[i]);
	} 
  } else {
	int isDecom = -1;
	for (int im = 0; im < _fNumParPCA; ++im) {
	  if(isDecomposed_PCA[im] == i) {isDecom = im;}
	}
	if(isDecom < 0) {
      MACH3LOG_ERROR("Setting {}(parameter {}) to fixed at {}", GetParName(i), i, _fCurrVal[i]);
	  //throw; 
	} else {
	  fParSigma_PCA[isDecom] *= -1.0;
      MACH3LOG_INFO("Setting un-decomposed {}(parameter {}/{} in PCA base) to fixed at {}", GetParName(i), i, isDecom, _fCurrVal[i]);
	}
  }
  
  return;
}
// ********************************************
void covarianceBase::setEvalLikelihood(const int i, const bool eL) {
// ********************************************
  if (i > _fNumPar) {
    MACH3LOG_INFO("Can't setEvalLikelihood for Cov={}/Param={} because size of Covariance = {}", getName(), i, _fNumPar);
    MACH3LOG_ERROR("Fix this in your config file please!");
    throw;
  } else {

    if(eL){
      MACH3LOG_INFO("Setting {} (parameter {}) to flat prior", GetParName(i), i);
    }
    else{
      // HW :: This is useful
      MACH3LOG_INFO("Setting {} (parameter {}) to non-flat prior", GetParName(i), i);
    }
    _fFlatPrior[i] = eL;
  }
}

// ********************************************
//KS: Custom function to perform multiplication of matrix and vector with multithreading
void covarianceBase::MatrixVectorMulti(double* _restrict_ VecMulti, double** _restrict_ matrix, const double* _restrict_ vector, const int n) {
// ********************************************

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for (int i = 0; i < n; ++i)
  {
    double result = 0.0;
    for (int j = 0; j < n; ++j)
    {
      result += matrix[i][j]*vector[j];
    }
    VecMulti[i] = result;
  }
}

// ********************************************
double covarianceBase::MatrixVectorMultiSingle(double** _restrict_ matrix, const double* _restrict_ vector, const int Length, const int i) {
// ********************************************

  double Element = 0.0;
  for (int j = 0; j < Length; ++j) {
    Element += matrix[i][j]*vector[j];
  }
  return Element;
}

// ********************************************
void covarianceBase::setIndivStepScale(std::vector<double> stepscale) {
// ********************************************

  if ((int)stepscale.size() != _fNumPar)
  {
    MACH3LOG_WARN("Stepscale vector not equal to number of parameters. Quitting..");
    MACH3LOG_WARN("Size of argument vector: {}", stepscale.size());
    MACH3LOG_WARN("Expected size: {}", _fNumPar);
    return;
  }

  for (int iParam = 0 ; iParam < _fNumPar; iParam++) {
    _fIndivStepScale[iParam] = stepscale[iParam];
  }

  printIndivStepScale();

  return;
}


void covarianceBase::printIndivStepScale() {
  std::cout << "============================================================" << std::endl;
  std::cout << std::setw(PrintLength) << "Parameter:" << " | " << std::setw(11) << "Step scale:" << std::endl;
  for (int iParam = 0; iParam < _fNumPar; iParam++) {
    std::cout << std::setw(PrintLength) << _fNames[iParam].c_str() << " | " << std::setw(11) << _fIndivStepScale[iParam] << std::endl;
  }
  std::cout << "============================================================" << std::endl;
}

//Makes sure that matrix is positive-definite (so no error is thrown when
//throwNominal() is called) by adding a small number to on-diagonal elements
void covarianceBase::MakePosDef(TMatrixDSym *cov) {
  //DB Save original warning state and then increase it in this function to suppress 'matrix not positive definite' messages
  //Means we no longer need to overload
  if(cov == NULL){
    cov = &*covMatrix;
  }

  int originalErrorWarning = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  //DB Loop 1000 times adding 1e-9 which tops out at 1e-6 shift on the diagonal before throwing error
  int MaxAttempts = 1e5;
  int iAttempt = 0;
  bool CanDecomp = false;
  TDecompChol chdcmp;
  
  for (iAttempt = 0; iAttempt < MaxAttempts; iAttempt++) {
    chdcmp = TDecompChol(*cov);
    if (chdcmp.Decompose()) {
      CanDecomp = true;
      break;
    } else {
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
      for (int iVar = 0 ; iVar < _fNumPar; iVar++) {
        (*cov)(iVar,iVar) += pow(10,-9);
      }
    }
  }

  if (!CanDecomp) {
    MACH3LOG_ERROR("Tried {} times to shift diagonal but still can not decompose the matrix", MaxAttempts);
    MACH3LOG_ERROR("This indicates that something is wrong with the input matrix");
    throw;
  }
  if(total_steps < 2) {
    MACH3LOG_INFO("Had to shift diagonal {} time(s) to allow the covariance matrix to be decomposed", iAttempt);
  }
  //DB Resetting warning level
  gErrorIgnoreLevel = originalErrorWarning;

  return;
}

void covarianceBase::resetIndivStepScale() {
  std::vector<double> stepScales(_fNumPar);
  for (int i = 0; i <_fNumPar; i++) {
    stepScales[i] = 1.;
  }
  _fGlobalStepScale = 1.0;
  setIndivStepScale(stepScales);
}

// HW: Code for throwing from separate throw matrix, needs to be set after init to ensure pos-def
void covarianceBase::setThrowMatrix(TMatrixDSym *cov){
   if (cov == NULL) {
    MACH3LOG_ERROR("Could not find covariance matrix you provided to setThrowMatrix");
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  if (covMatrix->GetNrows() != cov->GetNrows()) {
    MACH3LOG_ERROR("Matrix given for throw Matrix is not the same size as the covariance matrix stored in object!");
    std::cerr << "Stored covariance matrix size:" << covMatrix->GetNrows() << std::endl;
    std::cerr << "Given matrix size:" << cov->GetNrows() << std::endl;
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  throwMatrix = (TMatrixDSym*)cov->Clone();
  if(total_steps <= lower_adapt) makeClosestPosDef(throwMatrix);
  else MakePosDef(throwMatrix);
  
  TDecompChol TDecompChol_throwMatrix(*throwMatrix);
  
  if(!TDecompChol_throwMatrix.Decompose()) {
    MACH3LOG_ERROR("Cholesky decomposition failed for {} trying to make positive definite", matrixName);
    std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
    throw;
  }

  throwMatrix_CholDecomp = new TMatrixD(TDecompChol_throwMatrix.GetU());
  throwMatrix_CholDecomp->T();

  //KS: ROOT has bad memory management, using standard double means we can decrease most operation by factor 2 simply due to cache hits
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for (int i = 0; i < _fNumPar; ++i)
  {
    for (int j = 0; j < _fNumPar; ++j)
    {
      throwMatrixCholDecomp[i][j] = (*throwMatrix_CholDecomp)(i,j);
    }
  }
}

void covarianceBase::updateThrowMatrix(TMatrixDSym *cov){
  delete throwMatrix;
  throwMatrix = NULL;
  delete throwMatrix_CholDecomp;
  throwMatrix_CholDecomp = NULL;
  setThrowMatrix(cov);
}

// The setter
void covarianceBase::useSeparateThrowMatrix(TString throwMatrixFileName, TString throwMatrixName, TString meansVectorName){
// Firstly let's check if the file exists
  TFile* throwMatrixFile = new TFile(throwMatrixFileName);
  resetIndivStepScale();
  if(throwMatrixFile->IsZombie()) {
    MACH3LOG_ERROR("Couldn't find throw Matrix file : {}", throwMatrixFileName);
    std::cerr<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  } //We're done for now

  TMatrixDSym* tmp_throwMatrix = (TMatrixDSym*)throwMatrixFile->Get(throwMatrixName);
  TVectorD* tmp_meansvec = (TVectorD*)throwMatrixFile->Get(meansVectorName);

  if(!tmp_throwMatrix){
    std::cerr<<"ERROR : Couldn't find throw matrix "<<throwMatrixName<<" in "<<throwMatrixFileName<<std::endl;
    std::cerr<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  if(!tmp_meansvec){
    std::cerr<<"ERROR : Couldn't find means vector "<<meansVectorName<<" in "<<throwMatrixFileName<<std::endl;
    std::cerr<<__FILE__<<" : "<<__LINE__<<std::endl;
    throw;
  }

  adaptiveCovariance = (TMatrixDSym*)tmp_throwMatrix->Clone();
  par_means.resize(tmp_meansvec->GetNrows());
  for(int iMean = 0; iMean < tmp_meansvec->GetNrows(); iMean++){
      par_means[iMean] = (*tmp_meansvec)(iMean);
      updateThrowMatrix(adaptiveCovariance);
  }
  delete tmp_throwMatrix;
  delete tmp_meansvec;
  throwMatrixFile->Close();
  delete throwMatrixFile; // Just in case
}


void covarianceBase::useSeparateThrowMatrix(){
    initialiseNewAdaptiveChain(); 
    updateThrowMatrix(covMatrix);
}

//HW: Truly adaptive MCMC!
void covarianceBase::updateAdaptiveCovariance(){
  // https://projecteuclid.org/journals/bernoulli/volume-7/issue-2/An-adaptive-Metropolis-algorithm/bj/1080222083.full
  // Updates adaptive matrix
  // First we update the total means

#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for(int iRow = 0; iRow < _fNumPar; ++iRow){
    par_means_prev[iRow] = par_means[iRow];
    par_means[iRow] = (_fCurrVal[iRow]+par_means[iRow]*total_steps)/(total_steps+1);
  }

  //Now we update the covariances using cov(x,y)=E(xy)-E(x)E(y)
#ifdef MULTITHREAD
#pragma omp parallel for
#endif
  for(int iRow = 0; iRow < _fNumPar; ++iRow){
    for(int iCol = 0; iCol <= iRow; ++iCol){

      double cov_val = (*adaptiveCovariance)(iRow, iCol)*_fNumPar/5.6644;
      cov_val += par_means_prev[iRow]*par_means_prev[iCol]; //First we remove the current means
      cov_val = (cov_val*total_steps+_fCurrVal[iRow]*_fCurrVal[iCol])/(total_steps+1); //Now get mean(iRow*iCol)
      cov_val -= par_means[iCol]*par_means[iRow];
      cov_val*=5.6644/_fNumPar;
      (*adaptiveCovariance)(iRow, iCol) = cov_val;
      (*adaptiveCovariance)(iCol, iRow) = cov_val;
    }
  }
  //This is likely going to be the slow bit!
  total_steps += 1;
  if(total_steps == lower_adapt)
  {
    resetIndivStepScale();
  }

  if(total_steps >= lower_adapt) {
    updateThrowMatrix(adaptiveCovariance); //Now we update and continue!
  }
  // if(total_steps%1000==0 && total_steps>1000){
  //   suboptimality_vals.push_back(calculatesuboptimality(covSqrt, adaptiveCovariance));
  // }
}

void covarianceBase::initialiseNewAdaptiveChain(){
  // If we don't have a covariance matrix to start from for adaptive tune we need to make one!
  adaptiveCovariance = new TMatrixDSym(_fNumPar);
  par_means.resize(_fNumPar);
  par_means_prev.resize(_fNumPar);
  for(int i = 0; i < _fNumPar; ++i){
    par_means[i] = 0.;
    par_means_prev[i] = 0.;
    for(int j = 0; j <= i; ++j){
      (*adaptiveCovariance)(i,j) = 0.0; // Just make an empty matrix
      (*adaptiveCovariance)(j,i) = 0.0; // Just make an empty matrix
    }
  }
}

void covarianceBase::saveAdaptiveToFile(TString outFileName, TString systematicName){
  TFile* outFile = new TFile(outFileName, "UPDATE");
  if(outFile->IsZombie()){
    MACH3LOG_ERROR("Couldn't find {}", outFileName);
    throw;
  }
  TVectorD* outMeanVec = new TVectorD((int)par_means.size());
  for(int i = 0; i < (int)par_means.size(); i++){
    (*outMeanVec)(i)=par_means[i];
  }
  outFile->cd();
  adaptiveCovariance->Write(systematicName+"_postfit_matrix");
  outMeanVec->Write(systematicName+"_mean_vec");
  outFile->Close();
  delete outFile;
}

//HW: Finds closest possible positive definite matrix in Frobenius Norm ||.||_frob
// Where ||X||_frob=sqrt[sum_ij(x_ij^2)] (basically just turns an n,n matrix into vector in n^2 space
// then does Euclidean norm)
void covarianceBase::makeClosestPosDef(TMatrixDSym *cov)
{
  // Want to get cov' = (cov_sym+cov_polar)/2
  // cov_sym=(cov+cov^T)/2
  // cov_polar-> SVD cov to cov=USV^T then cov_polar=VSV^T

  //Get frob norm of cov
//  Double_t cov_norm=cov->E2Norm();
  
  TMatrixDSym* cov_trans = cov;
  cov_trans->T();
  TMatrixDSym cov_sym = 0.5*(*cov+*cov_trans); //If cov is symmetric does nothing, otherwise just ensures symmetry

  //Do SVD to get polar form
  TDecompSVD cov_sym_svd=TDecompSVD(cov_sym);
  if(!cov_sym_svd.Decompose()){
    MACH3LOG_ERROR("Cannot do SVD on input matrix!");
    throw;
  }
  
  TMatrixD cov_sym_v = (TMatrixD)cov_sym_svd.GetV();
  TMatrixD cov_sym_vt = cov_sym_v;
  cov_sym_vt.T();
  //SVD returns as vector (grrr) so need to get into matrix form for multiplying!
  TVectorD cov_sym_sigvect = (TVectorD)cov_sym_svd.GetSig();
  const Int_t nCols = cov_sym_v.GetNcols(); //square so only need rows hence lack of cols
  TMatrixDSym cov_sym_sig(nCols);
  TMatrixDDiag cov_sym_sig_diag(cov_sym_sig);
  cov_sym_sig_diag=cov_sym_sigvect;

  
  //Can finally get H=VSV
  TMatrixDSym cov_sym_polar = cov_sym_sig.SimilarityT(cov_sym_vt);//V*S*V^T (this took forver to find!)
  
  //Now we can construct closest approximater Ahat=0.5*(B+H)
  TMatrixDSym cov_closest_approx  = 0.5*(cov_sym+cov_sym_polar);//Not fully sure why this is even needed since symmetric B -> U=V
  //Get norm of transformed
//  Double_t approx_norm=cov_closest_approx.E2Norm();

  //std::cout<<"Initial Norm : "<<cov_norm<<" | Norm after transformation : "<<approx_norm<<" | Ratio : "<<cov_norm/approx_norm<<std::endl;
  
  *cov = cov_closest_approx;
  //Now can just add a makeposdef!
  MakePosDef(cov);
}
// ********************************************
std::vector<double> covarianceBase::getNominalArray() {
// ********************************************

  std::vector<double> nominal(_fNumPar);
  for (int i = 0; i < _fNumPar; ++i)
  {
    nominal[i] = _fPreFitValue[i];
  }
 return nominal;
}

// ********************************************
// KS: Convert covariance matrix to correlation matrix and return TH2D which can be used for fancy plotting
TH2D* covarianceBase::GetCorrelationMatrix() {
// ********************************************

  TH2D* hMatrix = new TH2D(getName(), getName(), _fNumPar, 0.0, _fNumPar, _fNumPar, 0.0, _fNumPar);

  for(int i = 0; i < _fNumPar; i++)
  {
    hMatrix->SetBinContent(i+1, i+1, 1.);
    hMatrix->GetXaxis()->SetBinLabel(i+1, GetParFancyName(i).c_str());
    hMatrix->GetYaxis()->SetBinLabel(i+1, GetParFancyName(i).c_str());
  }

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int i = 0; i < _fNumPar; i++)
  {
    for(int j = 0; j <= i; j++)
    {
      const double Corr = (*covMatrix)(i,j) / ( getDiagonalError(i) * getDiagonalError(j));
      hMatrix->SetBinContent(i+1, j+1, Corr);
      hMatrix->SetBinContent(j+1, i+1, Corr);
    }
  }
  return hMatrix;
}

// ********************************************
// KS: After step scale, prefit etc. value were modified save this modified config.
void covarianceBase::SaveUpdatedMatrixConfig() {
// ********************************************

  if (!_fYAMLDoc)
  {
    MACH3LOG_CRITICAL("Yaml node hasn't been initialised for matrix {}, something is not right", matrixName);
    MACH3LOG_CRITICAL("I am not throwing error but should be investigated");
    return;
  }

  YAML::Node copyNode = _fYAMLDoc;
  int i = 0;

  for (YAML::Node param : copyNode["Systematics"])
  {
    //KS: Feel free to update it, if you need updated prefit value etc
    param["Systematic"]["StepScale"]["MCMC"] = std::round(_fIndivStepScale[i] * 100.0) / 100.0; // Round to 2 decimal places

    i++;
  }
  // Save the modified node to a file
  std::ofstream fout("Modified_Matrix.yaml");
  fout << copyNode;
  fout.close();
}


#ifdef DEBUG_PCA
//KS: Let's dump all useful matrices to properly validate PCA
void covarianceBase::DebugPCA(const double sum, TMatrixD temp, TMatrixDSym submat)
{
  (void)submat;//This is used if DEBUG_PCA==2, this hack is to avoid compiler warnings
  TFile *PCA_Debug = new TFile("Debug_PCA.root", "RECREATE");
  PCA_Debug->cd();

  bool PlotText = true;
  //KS: If we have more than 200 plot becomes unreadable :(
  if(_fNumPar > 200) PlotText = false;

  TH1D* heigen_values = new TH1D("eigen_values", "Eigen Values", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  TH1D* heigen_cumulative = new TH1D("heigen_cumulative", "heigen_cumulative", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  TH1D* heigen_frac = new TH1D("heigen_fractional", "heigen_fractional", (int)eigen_values.GetNrows(), 0.0, (int)eigen_values.GetNrows());
  heigen_values->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative = 0;
  for(int i = 0; i < eigen_values.GetNrows(); i++)
  {
    heigen_values->SetBinContent(i+1, (eigen_values)(i));
    heigen_cumulative->SetBinContent(i+1, (eigen_values)(i)/sum + Cumulative);
    heigen_frac->SetBinContent(i+1, (eigen_values)(i)/sum);
    Cumulative += (eigen_values)(i)/sum;
  }
  heigen_values->Write("heigen_values");
  eigen_values.Write("eigen_values");
  heigen_cumulative->Write("heigen_values_cumulative");
  heigen_frac->Write("heigen_values_frac");

  TH2D* heigen_vectors = new TH2D(eigen_vectors);
  heigen_vectors->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors->Write("heigen_vectors");
  eigen_vectors.Write("eigen_vectors");

  TH2D* SubsetPCA = new TH2D(temp);
  SubsetPCA->GetXaxis()->SetTitle("Parameter in Normal Base");
  SubsetPCA->GetYaxis()->SetTitle("Parameter in Decomposed Base");

  SubsetPCA->Write("hSubsetPCA");
  temp.Write("SubsetPCA");
  TH2D* hTransferMat = new TH2D(TransferMat);
  hTransferMat->GetXaxis()->SetTitle("Parameter in Normal Base");
  hTransferMat->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  TH2D* hTransferMatT = new TH2D(TransferMatT);

  hTransferMatT->GetXaxis()->SetTitle("Parameter in Decomposed Base");
  hTransferMatT->GetYaxis()->SetTitle("Parameter in Normal Base");

  hTransferMat->Write("hTransferMat");
  TransferMat.Write("TransferMat");
  hTransferMatT->Write("hTransferMatT");
  TransferMatT.Write("TransferMatT");

  TCanvas *c1 = new TCanvas("c1"," ", 0, 0, 1024, 1024);
  c1->SetBottomMargin(0.1);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.05);
  c1->SetLeftMargin(0.12);
  c1->SetGrid();

  gStyle->SetPaintTextFormat("4.1f");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  // Make pretty correlation colors (red to blue)
  const int NRGBs = 5;
  TColor::InitializeColors();
  Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.25, 1.00, 1.00, 0.50 };
  Double_t green[NRGBs] = { 0.00, 0.25, 1.00, 0.25, 0.00 };
  Double_t blue[NRGBs]  = { 0.50, 1.00, 1.00, 0.25, 0.00 };
  TColor::CreateGradientColorTable(5, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  double maxz = 0;
  double minz = 0;

  c1->Print("Debug_PCA.pdf[");
  TLine *EigenLine = new TLine(nKeptPCApars, 0, nKeptPCApars, heigen_cumulative->GetMaximum());
  EigenLine->SetLineColor(kPink);
  EigenLine->SetLineWidth(2);
  EigenLine->SetLineStyle(kSolid);

  TText* text = new TText(0.5, 0.5, Form("Threshold = %g", eigen_threshold));
  text->SetTextFont (43);
  text->SetTextSize (40);

  heigen_values->SetLineColor(kRed);
  heigen_values->SetLineWidth(2);
  heigen_cumulative->SetLineColor(kGreen);
  heigen_cumulative->SetLineWidth(2);
  heigen_frac->SetLineColor(kBlue);
  heigen_frac->SetLineWidth(2);

  c1->SetLogy();
  heigen_values->SetMaximum(heigen_cumulative->GetMaximum()+heigen_cumulative->GetMaximum()*0.4);
  heigen_values->Draw();
  heigen_frac->Draw("SAME");
  heigen_cumulative->Draw("SAME");
  EigenLine->Draw("Same");
  text->DrawTextNDC(0.42, 0.84,Form("Threshold = %g", eigen_threshold));

  TLegend *leg = new TLegend(0.2, 0.2, 0.6, 0.5);
  leg->SetTextSize(0.04);
  leg->AddEntry(heigen_values, "Absolute", "l");
  leg->AddEntry(heigen_frac, "Fractional", "l");
  leg->AddEntry(heigen_cumulative, "Cumulative", "l");

  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->Draw("Same");

  c1->Print("Debug_PCA.pdf");
  c1->SetRightMargin(0.15);
  c1->SetLogy(0);
  delete EigenLine;
  delete leg;
  delete text;
  delete heigen_values;
  delete heigen_frac;
  delete heigen_cumulative;

  heigen_vectors->SetMarkerSize(0.2);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");

  TLine *Eigen_Line = new TLine(0, nKeptPCApars, LastPCAdpar-FirstPCAdpar, nKeptPCApars);
  Eigen_Line->SetLineColor(kGreen);
  Eigen_Line->SetLineWidth(2);
  Eigen_Line->SetLineStyle(kDotted);
  Eigen_Line->Draw("SAME");
  c1->Print("Debug_PCA.pdf");
  delete Eigen_Line;

  SubsetPCA->SetMarkerSize(0.2);
  minz = SubsetPCA->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else SubsetPCA->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) SubsetPCA->Draw("COLZ TEXT");
  else SubsetPCA->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete SubsetPCA;

  hTransferMat->SetMarkerSize(0.15);
  minz = hTransferMat->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMat->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMat->Draw("COLZ TEXT");
  else hTransferMat->Draw("COLZ");
  c1->Print("Debug_PCA.pdf");
  delete hTransferMat;

  hTransferMatT->SetMarkerSize(0.15);
  minz = hTransferMatT->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else hTransferMatT->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) hTransferMatT->Draw("COLZ TEXT");
  else hTransferMatT->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete hTransferMatT;


  //KS: Crosscheck against Eigen library
  #if DEBUG_PCA == 2
  Eigen::MatrixXd Submat_Eigen(submat.GetNrows(), submat.GetNcols());

  #ifdef MULTITHREAD
  #pragma omp parallel for
  #endif
  for(int i = 0; i < submat.GetNrows(); i++)
  {
    for(int j = 0; j < submat.GetNcols(); j++)
    {
      Submat_Eigen(i,j) = (submat)(i,j);
    }
  }
  Eigen::EigenSolver<Eigen::MatrixXd> EigenSolver;
  EigenSolver.compute(Submat_Eigen);
  Eigen::VectorXd eigen_val = EigenSolver.eigenvalues().real();
  Eigen::MatrixXd eigen_vect = EigenSolver.eigenvectors().real();
  std::vector<std::tuple<double, Eigen::VectorXd>> eigen_vectors_and_values;
  double Sum_Eigen = 0;
  for(int i = 0; i < eigen_val.size(); i++)
  {
    std::tuple<double, Eigen::VectorXd> vec_and_val(eigen_val[i], eigen_vect.row(i));
    eigen_vectors_and_values.push_back(vec_and_val);
    Sum_Eigen += eigen_val[i];
  }
  std::sort(eigen_vectors_and_values.begin(), eigen_vectors_and_values.end(),
            [&](const std::tuple<double, Eigen::VectorXd>& a, const std::tuple<double, Eigen::VectorXd>& b) -> bool
            { return std::get<0>(a) > std::get<0>(b); } );
  int index = 0;
  for(auto const vect : eigen_vectors_and_values)
  {
    eigen_val(index) = std::get<0>(vect);
    eigen_vect.row(index) = std::get<1>(vect);
    index++;
  }
  TH1D* heigen_values_Eigen = new TH1D("eig_values", "Eigen Values", eigen_val.size(), 0.0, eigen_val.size());
  TH1D* heigen_cumulative_Eigen = new TH1D("eig_cumulative", "heigen_cumulative", eigen_val.size(), 0.0, eigen_val.size());
  TH1D* heigen_frac_Eigen = new TH1D("eig_fractional", "heigen_fractional", eigen_val.size(), 0.0, eigen_val.size());
  heigen_values_Eigen->GetXaxis()->SetTitle("Eigen Vector");
  heigen_values_Eigen->GetYaxis()->SetTitle("Eigen Value");

  double Cumulative_Eigen = 0;
  for(int i = 0; i < eigen_val.size(); i++)
  {
    heigen_values_Eigen->SetBinContent(i+1, eigen_val(i));
    heigen_cumulative_Eigen->SetBinContent(i+1, eigen_val(i)/sum + Cumulative_Eigen);
    heigen_frac_Eigen->SetBinContent(i+1, eigen_val(i)/sum);
    Cumulative_Eigen += eigen_val(i)/sum;
  }
  heigen_values_Eigen->SetLineColor(kRed);
  heigen_values_Eigen->SetLineWidth(2);
  heigen_cumulative_Eigen->SetLineColor(kGreen);
  heigen_cumulative_Eigen->SetLineWidth(2);
  heigen_frac_Eigen->SetLineColor(kBlue);
  heigen_frac_Eigen->SetLineWidth(2);

  c1->SetLogy();
  heigen_values_Eigen->SetMaximum(heigen_cumulative_Eigen->GetMaximum()+heigen_cumulative_Eigen->GetMaximum()*0.4);
  heigen_values_Eigen->Draw();
  heigen_cumulative_Eigen->Draw("SAME");
  heigen_frac_Eigen->Draw("SAME");

  TLegend *leg_Eigen = new TLegend(0.2, 0.2, 0.6, 0.5);
  leg_Eigen->SetTextSize(0.04);
  leg_Eigen->AddEntry(heigen_values_Eigen, "Absolute", "l");
  leg_Eigen->AddEntry(heigen_frac_Eigen, "Fractional", "l");
  leg_Eigen->AddEntry(heigen_cumulative_Eigen, "Cumulative", "l");

  leg_Eigen->SetLineColor(0);
  leg_Eigen->SetLineStyle(0);
  leg_Eigen->SetFillColor(0);
  leg_Eigen->SetFillStyle(0);
  leg_Eigen->Draw("Same");

  c1->Print( "Debug_PCA.pdf");
  c1->SetLogy(0);
  delete heigen_values_Eigen;
  delete heigen_cumulative_Eigen;
  delete heigen_frac_Eigen;
  delete leg_Eigen;

  TH2D* heigen_vectors_Eigen = new TH2D("Eigen_Vectors", "Eigen_Vectors", eigen_val.size(), 0.0, eigen_val.size(), eigen_val.size(), 0.0, eigen_val.size());

  for(int i = 0; i < eigen_val.size(); i++)
  {
    for(int j = 0; j < eigen_val.size(); j++)
    {
      //KS: +1 because there is offset in histogram relative to TMatrix
      heigen_vectors_Eigen->SetBinContent(i+1,j+1, eigen_vect(i,j));
    }
  }
  heigen_vectors_Eigen->GetXaxis()->SetTitle("Parameter in Normal Base");
  heigen_vectors_Eigen->GetYaxis()->SetTitle("Parameter in Decomposed Base");
  heigen_vectors_Eigen->SetMarkerSize(0.15);
  minz = heigen_vectors_Eigen->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors_Eigen->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors_Eigen->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));

  if(PlotText) heigen_vectors_Eigen->Draw("COLZ TEXT");
  else heigen_vectors_Eigen->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");

  heigen_vectors->SetTitle("ROOT minus Eigen");
  heigen_vectors->Add(heigen_vectors_Eigen, -1.);
  minz = heigen_vectors->GetMinimum();
  if (fabs(0-maxz)>fabs(0-minz)) heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-maxz),0+fabs(0-maxz));
  else heigen_vectors->GetZaxis()->SetRangeUser(0-fabs(0-minz),0+fabs(0-minz));
  if(PlotText) heigen_vectors->Draw("COLZ TEXT");
  else heigen_vectors->Draw("COLZ");
  c1->Print( "Debug_PCA.pdf");
  delete heigen_vectors_Eigen;

  #endif
  delete heigen_vectors;

  c1->Print( "Debug_PCA.pdf]");
  delete c1;
  PCA_Debug->Close();
  delete PCA_Debug;
}
#endif
