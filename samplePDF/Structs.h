#ifndef _Structs_h_
#define _Structs_h_

// Run low or high memory versions of structs
// N.B. for 64 bit systems sizeof(float) == sizeof(double) so not a huge effect
#define __LOW__MEMORY_STRUCTS__

#ifdef __LOW_MEMORY_STRUCTS__
#define __float__ float
#define __int__ short int
#else
#define __float__ double
#define __int__ int
#endif

// Include some healthy defines for constructors
#define __BAD_DOUBLE__ -999.99
#define __BAD_INT__ -999

#define __LARGE_WEIGHT__ 100

#define __TH2PolyOverflowBins__ 9

#define __DEFAULT_RETURN_VAL__ -999999.123456

// C++ includes
#include <sstream>
#include <fstream> 
#include <iostream>
#include <vector>
#include <iomanip>

// ROOT include
#include "set"
#include "TSpline.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TObjString.h"
#include "unordered_map"
#include "TH2Poly.h"
#include "list"

// *******************
// Template to make vector out of an array of any length
template< typename T, size_t N >
std::vector<T> MakeVector( const T (&data)[N] ) {
// *******************
  return std::vector<T>(data, data+N);
}

// *******************
// Normalisations for cross-section parameters
class XsecNorms2 {
  // *******************
  public:
    // Bins for normalisation parameter
    TAxis *ebins;
    // Name of parameters
    std::string name;
    // Mode which parameter applies to
    int mode;
    // PDG which parameter applies to
    std::vector<int> pdgs;
    // Targets which parameter applies to
    std::vector<int> targets;
    // Parameter number of this normalisation in current NIWG parameterisation
    int startbin;
};

// *******************
// Normalisations for cross-section parameters
class XsecNorms3 {
  // *******************
  public:
    // Bins for normalisation parameter
    TAxis *ebins;
    // Name of parameters
    std::string name;
    // Mode which parameter applies to
    std::vector<int> modes;
    // Horn currents which parameter applies to
    std::vector<int> horncurrents;
    // PDG which parameter applies to
    std::vector<int> pdgs;
    // Preosc PDG which parameter applies to
    std::vector<int> preoscpdgs;
    // Targets which parameter applies to
    std::vector<int> targets;
    //Does this parameter have kinematic bounds
    bool hasKinBounds;
    //Etrue bounds
    double etru_bnd_low;
    double etru_bnd_high;
    //Etrue bounds
    double q2_true_bnd_low;
    double q2_true_bnd_high;
    // Parameter number of this normalisation in current NIWG parameterisation
    int index;
};

// *******************
// ETA - new version of this for refactor
// Normalisations for cross-section parameters
// Carrier for whether you want to apply a systematic to an event or not
class XsecNorms4 {
  // *******************
  public:
    // Bins for normalisation parameter
    TAxis *ebins;
    // Name of parameters
    std::string name;
    // Mode which parameter applies to
    std::vector<int> modes;
    // Horn currents which parameter applies to
    std::vector<int> horncurrents;
    // PDG which parameter applies to
    std::vector<int> pdgs;
    // Preosc PDG which parameter applies to
    std::vector<int> preoscpdgs;
    // Targets which parameter applies to
    std::vector<int> targets;
    //Does this parameter have kinematic bounds
    bool hasKinBounds;
	//Generic vector contain enum relating to a kinematic variable
    //and lower and upper bounds. This can then be passed to IsEventSelected
	std::vector< std::vector<double> > Selection;
	
	//Generic vector containing the string of kinematic type
	//This then needs to be converted to a kinematic type enum 
	//within a samplePDF daughter class
	//The bounds for each kinematic variable are given in Selection
	std::vector< std::string > KinematicVarStr;
	
    // Parameter number of this normalisation in current NIWG parameterisation
    int index;
};

// *******************
// Add a struct to hold info about the splinified xsec parameters
// Used in ND280 code to accelerate CPU TSpline3 evaluation
// Adaptable for SK
struct FastSplineInfo {
// *******************
  // Number of points in spline
  int nPts;

  // Array of the knots positions
  double *xPts;

  // Array of what segment of spline we're currently interested in
  // Gets updated once per MCMC iteration
  int CurrSegment;

  //ETA trying to change things so we don't read in all flat splines
  int flat;
};

// *******************
// The kinematic types to be enjoyed by all
enum KinematicTypes {
  // ******************
  kLeptonMomentum        = 0,
  kLeptonCosTheta        = 1,
  kNeutrinoEnergy        = 2,
  kQ2                    = 3,
  kErecQE                = 4,
  kQ2QE                  = 5,
  kLeptonTheta           = 6,
  kq0                    = 7,
  kq3                    = 8,
  kProtonMomentumReco    = 9,
  kProtonCosThetaReco    = 10,
  kPionMomentum          = 11,
  kPionCosTheta          = 12,
  kPionMomentumFS        = 13,
  kProtonMomentumFS      = 14,
  kLeptonTrueTheta       = 15,
  kW                     = 16,
  kW2                    = 17,
  kM3Mode                = 18,
  kOscChannel            = 19,
  kErec                  = 20,
  //kPDFBinning has to be the final one
  kPDFBinning            = 21,
  kNKinematicParams   
};


// *******************
// Get the name of a kinematic type
// Useful for setting x and y axis of TH2D samples and so on
inline std::string Kinematic_ToString(KinematicTypes type) {
  // *******************
  std::string ReturnString = "";
  switch(type) {
    case kLeptonMomentum:
      ReturnString = "p_{#mu} (MeV)";
      break;
    case kLeptonCosTheta:
      ReturnString = "cos #theta_{#mu}";
      break;
    case kNeutrinoEnergy:
      ReturnString = "E_{#nu}^{true} (GeV)";
      break;
    case kQ2:
      ReturnString = "Q^{2}_{true} (GeV^{2})";
      break;
    case kErecQE:
      ReturnString = "E_{#nu}^{QE} (GeV)";
      break;
    case kQ2QE:
      ReturnString = "Q^{2}_{QE} (GeV^{2})";
      break;
    case kLeptonTheta:
      ReturnString = "#theta_{lep} (degrees)";
      break;
    case kq0:
      ReturnString = "|q_{0}| (Mev)";
      break;
    case kq3:
      ReturnString = "|q_{3}| (Gev)";
      break;
    case kProtonMomentumReco:
      ReturnString = "p_{proton} (MeV)";
      break;
    case kProtonCosThetaReco:
      ReturnString = "cos #theta_{proton}";
      break;
    case kPionMomentum:
      ReturnString = "High. interaction p_{#pi^{#pm}}(MeV)";
      break;
    case kPionCosTheta:
      ReturnString = "cos #theta_{#pi}";
      break;
    case kPionMomentumFS:
      ReturnString = "High. FS p_{#pi^{#pm}}(MeV)";
      break;
    case kProtonMomentumFS:
      ReturnString = "True proton momentum in Final State (MeV)";
      break;
    case kLeptonTrueTheta:
      ReturnString = "True #theta_{lep} (degrees)";
      break;
    case kW2:
      ReturnString = "W^{2} (GeV^{2})";
      break;
    case kW:
      ReturnString = "W (GeV)";
      break;
    case kM3Mode:
      ReturnString = "Mode";
      break;
    case kOscChannel:
      ReturnString = "Oscillation Channel";
      break;
    case kPDFBinning:
      std::cout << "You shouldn't be here because PDF binning is different for each sample." << std::endl;
      std::cout << "Given kPDFBinning. Exitting.." << std::endl;
      throw;
    default:
      std::cerr << "Error, did not find ToString conversion for " << type << std::endl;
      break;
  }
  return ReturnString;
};

inline double StringToKinematicVar(std::string kinematic_str){

  if(strcmp(kinematic_str.c_str(), "Etrue") == 0){
	return kNeutrinoEnergy;
  }
  else if(strcmp(kinematic_str.c_str(), "Q2true") == 0){
	return kQ2;
  }	
  else{
	std::cerr << "[ERROR] " << __FILE__ << ":" << __LINE__ << " Unknown string " << kinematic_str << " passed from xsec cov" << std::endl; 
	std::cerr << "Have you added this kinematic parameter into StringToKinematicVar?" << std::endl;
	throw;
  }

  return 0;
}



// Useful for printing out ot pngs or writing out to histograms
inline std::string Kinematic_ToShortString(KinematicTypes type) {
  // *******************
  std::string ReturnString = "";
  switch(type) {
    case kLeptonMomentum:
      ReturnString = "LepMom";
      break;
    case kLeptonCosTheta:
      ReturnString = "CosThetaLep";
      break;
    case kNeutrinoEnergy:
      ReturnString = "TrueEnu";
      break;
    case kQ2:
      ReturnString = "TrueQ2";
      break;
    case kErecQE:
      ReturnString = "TrueEnuQE";
      break;
    case kQ2QE:
      ReturnString = "TrueQ2QE";
      break;
    case kLeptonTheta:
      ReturnString = "ThetaLep";
      break;
    case kq0:
      ReturnString = "q0";
      break;
    case kq3:
      ReturnString = "q3";
      break;
    case kProtonMomentumReco:
      ReturnString = "ProtonMom";
      break; 
    case kProtonCosThetaReco:
      ReturnString = "CosThetaProton";
      break; 
    case kPionMomentum:
      ReturnString = "PiMomInt";
      break;
    case kPionCosTheta:
      ReturnString = "CosThetaPi";
      break;
    case kPionMomentumFS:
      ReturnString = "LeadingPiMomFS";
      break;
    case kProtonMomentumFS:
      ReturnString = "ProtonMomFS";
      break;
    case kLeptonTrueTheta:
      ReturnString = "TrueThetaLep";
      break;
    case kW2:
      ReturnString = "TrueW2";
      break;
    case kW:
      ReturnString = "TrueW";
      break;
    case kM3Mode:
      ReturnString = "Mode";
      break;
    case kOscChannel:
      ReturnString = "OscChannel";
      break;
    case kErec:
      ReturnString = "Erec";
      break;
    case kPDFBinning:
      std::cout << "You shouldn't be here because PDF binning is different for each sample." << std::endl;
      std::cout << "Given kPDFBinning. Exitting.." << std::endl;
      throw;
    default:
      std::cerr << "Error, did not find ToShortString conversion for " << type << std::endl;
      break;
  }
  return ReturnString;
};

// ********************************************
// Generic xsec class
// Can use TF1 or TSpline3 or TSpline5 here, tjoho
template <class T>
class XSecStruct{
  // ********************************************
  public:
    // The light constructor
    XSecStruct(__int__ NumberOfSplines) {
      nParams = NumberOfSplines;
      Func.reserve(nParams);
      for (int i = 0; i < nParams; ++i) {
        Func[i] = NULL;
      }
    }

    // The empty constructor
    XSecStruct() {
      nParams = 0;
      Func = NULL;
    };

    // The light destructor
    ~XSecStruct() {
      for (int i = 0; i < nParams; ++i) {
        if (Func[i]) delete Func[i];
      }
    }

    // Get number of splines
    inline __int__ GetNumberOfParams() { return nParams; }

    // The Printer
    inline void Print() {
      std::cout << "    Splines: " << std::endl;
      for (int i = 0; i < nParams; ++i) {
        if (!Func[i]) continue;
        std::cout << "    " << std::left << std::setw(25) << Func[i]->GetName() << std::endl;
      }
    }

    // Set the number of splines for this event
    inline void SetSplineNumber(__int__ NumberOfSplines) {
      nParams = NumberOfSplines;
      Func = new T[nParams];
    }

    // Get the function for the nth spline
    inline T GetFunc(__int__ nSpline) { return Func[nSpline]; }
    // Set the function for the nth spline
    inline void SetFunc(__int__ nSpline, T Function) { Func[nSpline] = Function; }
    // Eval the current variation
    inline double Eval(__int__ nSpline, __float__ variation) {
      // Some will be NULL, check this
      if (Func[nSpline]) {
        return Func[nSpline]->Eval(variation);
      } else {
        return 1.0;
      }
    }
  private:
    // Number of parameters
    __int__ nParams;
    // The function
    T* Func;
};

// ************************
// A reduced TF1 class only
// Only saves parameters for each TF1 and how many parameters each parameter set has
class TF1_red {
// ************************
  public:
    // Empty constructor
    TF1_red() {
      length = 0;
      Par = NULL;
    }

    // Empty destructor
    ~TF1_red() {
      if (Par != NULL) {
        delete[] Par;
        Par = NULL;
      }
    }

    // The useful constructor with deep copy
    TF1_red(__int__ nSize, __float__* Array, __int__ Parameter) {
      length = nSize;
      for (int i = 0; i < length; ++i) {
        Par[i] = Array[i];
      }
      ParamNo = Parameter;
    }

    // The TF1 constructor with deep copy
    TF1_red(TF1* &Function, int Param = -1) {
      Par = NULL;
      SetFunc(Function, Param);
    }

    // Get the number
    inline std::string GetName() {
      std::stringstream ss;
      ss << ParamNo;
      return ss.str();
    }

    // Set the function
    inline void SetFunc(TF1* &Func, int Param = -1) {
      length = Func->GetNpar();
      if (Par != NULL) delete[] Par;
      Par = new __float__[length];
      for (int i = 0; i < length; ++i) {
        Par[i] = Func->GetParameter(i);
      }
      ParamNo = Param;
      delete Func;
      Func = NULL;
    }

    // Evaluate a variation
    inline double Eval(__float__ var) {
      // If we have 5 parameters we're using a fifth order polynomial
      if (length == 5) {
        return 1+Par[0]*var+Par[1]*var*var+Par[2]*var*var*var+Par[3]*var*var*var*var+Par[4]*var*var*var*var*var;
      // If we have 2 parameters we're using two linear equations
      } else if (length == 2) {
        return (var<=0)*(1+Par[0]*var)+(var>0)*(1+Par[1]*var);
      } else {
        std::cerr << "*** Error in reduced TF1 class!" << std::endl;
        std::cerr << "    Class only knows about 5th order polynomial and two superposed linear function" << std::endl;
        std::cerr << "    You have tried something else than this, which remains unimplemented" << std::endl;
        std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
        throw;
      }
    }

    // Set a parameter to a value
    inline void SetParameter(__int__ Parameter, __float__ Value) {
      Par[Parameter] = Value;
    }

    // Get a parameter value
    double GetParameter(__int__ Parameter) {
      if (Parameter > length) {
        std::cerr << "Error: you requested parameter number " << Parameter << " but length is " << length << " parameters" << std::endl;
        throw;
        return -999.999;
      }
      return Par[Parameter];
    }

    // Set the size
    inline void SetSize(__int__ nSpline) {
      length = nSpline;
      Par = new __float__[length];
    }
    // Get the size
    inline int GetSize() { return length; }
    inline void Print() {
      std::cout << "Printing TF1_red: " << std::endl;
      std::cout << "  ParamNo = " << ParamNo << std::endl;
      std::cout << "  Length  = " << length << std::endl;
      std::cout << "  a       = " << Par[0] << std::endl;
      std::cout << "  b       = " << Par[1] << std::endl;
      if (length == 5) {
        std::cout << "  c       = " << Par[2] << std::endl;
        std::cout << "  d       = " << Par[3] << std::endl;
        std::cout << "  e       = " << Par[4] << std::endl;
      }
    }

  private:
    // The parameters
    __float__* Par;
    __int__ length;
    // Save the parameter number this spline applies to
    __int__ ParamNo;
};

// ************************
// Reduced TSpline3 class
class TSpline3_red {
// ************************
  public:
    // Empty constructor
    TSpline3_red() {
      nPoints = 0;
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    TSpline3_red(TSpline3* &spline, int Param = -1) {
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
      SetFunc<TSpline3>(spline, Param);
    }

    // Copy constructor, should also accept and copy derived classes, i.e.  akima and monotone splinese
    TSpline3_red(TSpline3_red* &spline, int Param = -1) {
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
      SetFunc<TSpline3_red>(spline, Param);
    }

    // constructor taking parameters
    TSpline3_red(__float__ *X, __float__ *Y, __int__ N, __float__ **P, __int__ parNo){
      nPoints = N;
      ParamNo = parNo;
      // std::cout<<"nPoints: "<<nPoints<<std::endl;
      // Save the parameters for each knot
      Par = new __float__*[nPoints];
      // Save the positions of the knots
      XPos = new __float__[nPoints];
      // Save the y response at each knot
      YResp = new __float__[nPoints];
      for(int j=0; j<N; ++j){
        Par[j] = new __float__[3];
        Par[j][0] = P[j][0];
        Par[j][1] = P[j][1];
        Par[j][2] = P[j][2];
        XPos[j] = X[j];
        YResp[j] = Y[j];

        if((Par[j][0] == -999) | (Par[j][1] ==-999) | (Par[j][2] ==-999) | (XPos[j] ==-999) | (YResp[j] ==-999)){
          std::cerr<<"******************* Bad parameter values when construction TSpline3_red *********************" <<std::endl;
          std::cerr<<"passed val (i, x, y, b, c, d): "<<j<<", "<<X[j]<<", "<<Y[j]<<", "<<P[j][0]<<", "<<P[j][1]<<", "<<P[j][2]<<std::endl;
          std::cerr<<"set values (i, x, y, b, c, d): "<<j<<", "<<XPos[j]<<", "<<YResp[j]<<", "<<Par[j][0]<<", "<<Par[j][1]<<", "<<Par[j][2]<<std::endl;
          std::cerr<<"*********************************************************************************************" <<std::endl;
        }
      }
    }

    // Set a function
    template <typename T>
    inline void SetFunc(T* &spline, int Param = -1) {
      nPoints = spline->GetNp();
      ParamNo = Param;
      if (Par != NULL) {
        for (int i = 0; i < nPoints; ++i) {
          delete[] Par[i];
          Par[i] = NULL;
        }
        delete[] Par;
        Par = NULL;
      }
      if (XPos != NULL) delete[] XPos;
      if (YResp != NULL) delete[] YResp;
      // Save the parameters for each knot
      Par = new __float__*[nPoints];
      // Save the positions of the knots
      XPos = new __float__[nPoints];
      // Save the y response at each knot
      YResp = new __float__[nPoints];
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new __float__[3];
        double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
        spline->GetCoeff(i, x, y, b, c, d);
        XPos[i]   = x;
        YResp[i]  = y;
        Par[i][0] = b;
        Par[i][1] = c;
        Par[i][2] = d;
      }
      delete spline;
      spline = NULL;
    }

    // Empty destructor
    ~TSpline3_red() {
      for (int i = 0; i < nPoints; ++i) {
        if (Par[i] != NULL) {
          delete[] Par[i];
        }
      }
      delete[] Par;
      delete[] XPos;
      delete[] YResp;
      Par = NULL;
      XPos = YResp = NULL;
    }

    // Find the segment relevant to this variation in x
    // See root/hist/hist/src/TSpline3::FindX(double) or samplePDFND....::FindSplineSegment
    inline int FindX(double x) {
      // The segment we're interested in (klow in ROOT code)
      int segment = 0;
      int kHigh = nPoints-1;
      // If the variation is below the lowest saved spline point
      if (x <= XPos[0]){
        segment = 0;
        // If the variation is above the highest saved spline point
      } else if (x >= XPos[nPoints-1]) {
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
          if (x > XPos[kHalf]) {
            segment = kHalf;
            // Else move kHigh down
          } else {
            kHigh = kHalf;
          }
        } // End the while: we've now done our binary search
      } // End the else: we've now found our point
      if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;
      return segment;
    }

    // Evaluate the weight from a variation
    inline double Eval(double var) {
      // Get the segment for this variation
      int segment = FindX(var);
      // The get the coefficients for this variation
      double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
      GetCoeff(segment, x, y, b, c, d);
      double dx = var - x;
      // Evaluate the third order polynomial
      double weight = y+dx*(b+dx*(c+d*dx));
      return weight;
    }

    // Get the number of points
    inline int GetNp() { return nPoints; }
    // Get the ith knot's x and y position
    inline void GetKnot(int i, double &xtmp, double &ytmp) {
      xtmp = XPos[i];
      ytmp = YResp[i];
    }

    // Get the coefficient of a given segment
    inline void GetCoeff(int segment, double &x, double &y, double &b, double &c, double &d) {
      b = Par[segment][0];
      c = Par[segment][1];
      d = Par[segment][2];
      x = XPos[segment];
      y = YResp[segment];
    }

    // Get the number
    inline std::string GetName() {
      std::stringstream ss;
      ss << ParamNo;
      return ss.str();
    }

    // Make a TSpline3 from the reduced splines
    inline TSpline3* ConstructTSpline3() {
      TSpline3 *spline = new TSpline3(GetName().c_str(), XPos, YResp, nPoints);
      return spline;
    }

    // Make a TSpline3_red
    inline TSpline3_red* ConstructTSpline3_red() {
      TSpline3_red *spline = new TSpline3_red(XPos, YResp, nPoints, Par, ParamNo);
      return spline;
    }

  protected: //changed to protected from private so can be accessed by derived classes
    // Number of points/knot in TSpline3
    __int__ nPoints;
    // Always uses a third order polynomial, so hard-code the number of coefficients in implementation
    __float__ **Par;
    // Positions of each x for each knot
    __float__ *XPos;
    // y-value for each knot
    __float__ *YResp;
    // Parameter number (which parameter is this spline for)
    __int__ ParamNo;
};


// ************************
// Akima Spline class
class Akima_Spline: public TSpline3_red {
// ************************
// closely follows TSpline3_red class to fit in easily with existing machinery
// Akima spline is similar to regular cubic spline but is allowed to be discontinuous in 2nd derivative and coefficients in any segment
// only depend on th 2 nearest points on either side

  public:
    // Empty constructor
    Akima_Spline()
    :TSpline3_red()
    {
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    Akima_Spline(TSpline3* &spline, int Param = -1)
    { // need to override this so that daughter class SetFunc Gets called instead of TSpline3_red's version
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
      SetFunc(spline, Param);
    }

    // Empty destructor
    ~Akima_Spline() {
      // this should call base class destructor automatically
    }

    // Set a function
    inline void SetFunc(TSpline3* &spline, int Param = -1) {
      nPoints = spline->GetNp();
      ParamNo = Param;
      if (Par != NULL) {
        for (int i = 0; i < nPoints; ++i) {
          delete[] Par[i];
          Par[i] = NULL;
        }
        delete[] Par;
        Par = NULL;
      }
      if (XPos != NULL) delete[] XPos;
      if (YResp != NULL) delete[] YResp;
      // Save the parameters for each knot
      Par = new __float__*[nPoints];
      // Save the positions of the knots
      XPos = new __float__[nPoints];
      // Save the y response at each knot
      YResp = new __float__[nPoints];

      // get the knot values for the spline
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new __float__[3];

        double x = -999.99, y = -999.99;
        spline->GetKnot(i, x, y);

        XPos[i]   = x;
        YResp[i]  = y;
      }

      __float__ mvals[nPoints + 2];
      __float__ svals[nPoints];

      for (int i = -2; i <= nPoints; ++i) {
        // if segment is first or last or 2nd to first or last, needs to be dealt with slightly differently;
        // need to estimate the values for additinal points which would lie outside of the spline
        if(i ==-2){
          mvals[i+2] = 3.0 * (YResp[1] - YResp[0]) / (XPos[1] - XPos[0]) - 2.0*(YResp[2] - YResp[1]) / (XPos[2] - XPos[1]);
        }
        else if(i==-1){
          mvals[i+2] = 2.0 * (YResp[1] - YResp[0]) / (XPos[1] - XPos[0]) - (YResp[2] - YResp[1]) / (XPos[2] - XPos[1]);
        }
        else if(i==nPoints){
          mvals[i+2] = 3.0 * (YResp[nPoints-1] - YResp[nPoints-2]) / (XPos[nPoints-1] - XPos[nPoints-2]) - 2.0*(YResp[nPoints-2] - YResp[nPoints-3]) / (XPos[nPoints-2] - XPos[nPoints-3]);
        }
        else if(i == nPoints - 1){
          mvals[i+2] = 2.0 * (YResp[nPoints-1] - YResp[nPoints-2]) / (XPos[nPoints-1] - XPos[nPoints-2]) - (YResp[nPoints-2] - YResp[nPoints-3]) / (XPos[nPoints-2] - XPos[nPoints-3]);
        }
        //standard internal segment
        else{
          mvals[i+2] = (YResp[i+1] - YResp[i])/ (XPos[i+1] - XPos[i]);
        }
      }

      for(int i =2; i<=nPoints+2; i++){
        if (abs(mvals[i+1] - mvals[i]) + abs(mvals[i-1] - mvals[i-2]) != 0.0){
          svals[i-2] = (abs(mvals[i+1] - mvals[i]) * mvals[i-1] + abs(mvals[i-1] - mvals[i-2]) *mvals[i]) / (abs(mvals[i+1] - mvals[i]) + abs(mvals[i-1] - mvals[i-2]));
          }
        else{svals[i-2] = mvals[i];}
      }

      // calculate the coefficients for the spline
      for(int i = 0; i <nPoints; i++){
        __float__ b, c, d = -999.999;

        b = svals[i];
        c = (3.0* (YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i]) -2.0 *svals[i] - svals[i +1]) /(XPos[i+1] - XPos[i]);
        d = ((svals[i + 1] +svals[i]) - 2.0*(YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i])) / ((XPos[i+1] - XPos[i]) * (XPos[i+1] - XPos[i]));

        Par[i][0] = b;
        Par[i][1] = c;
        Par[i][2] = d;
      }
      delete spline;
      spline = NULL;
    }
    // finished calculating coeffs
};




// ************************
// Monotone Spline class
class Monotone_Spline: public TSpline3_red {
// ************************
// closely follows TSpline3_red class to fit in easily with existing machinery
// Akima spline is similar to regular cubic spline but enforce the condition that the interpolated value at any point
// must be between its two nearest knots, DOES NOT make the entire spline monotone, only the segments

  public:
    // Empty constructor
    Monotone_Spline()
    :TSpline3_red()
    {
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    Monotone_Spline(TSpline3* &spline, int Param = -1)
    { // need to override this so that daughter class SetFunc Gets called instead of TSpline3_red's version
      Par = NULL;
      XPos = NULL;
      YResp = NULL;
      SetFunc(spline, Param);
    }

    // Empty destructor
    ~Monotone_Spline() {
      // this should call base class destructor automatically
    }

    // Set a function
    inline void SetFunc(TSpline3* &spline, int Param = -1) {
      nPoints = spline->GetNp();
      ParamNo = Param;
      if (Par != NULL) {
        for (int i = 0; i < nPoints; ++i) {
          delete[] Par[i];
          Par[i] = NULL;
        }
        delete[] Par;
        Par = NULL;
      }
      if (XPos != NULL) delete[] XPos;
      if (YResp != NULL) delete[] YResp;
      // Save the parameters for each knot
      Par = new __float__*[nPoints];
      // Save the positions of the knots
      XPos = new __float__[nPoints];
      // Save the y response at each knot
      YResp = new __float__[nPoints];
      // values of the secants at each point (for calculating monotone spline)
      Secants = new __float__[nPoints -1];
      // values of the tangens at each point (for calculating monotone spline)
      Tangents = new __float__[nPoints];

      // get the knot values for the spline
      for (int i = 0; i < nPoints; ++i) {
        // 3 is the size of the TSpline3 coefficients
        Par[i] = new __float__[3];

        double x = -999.99, y = -999.99;
        spline->GetKnot(i, x, y);

        XPos[i]   = x;
        YResp[i]  = y;

        Tangents[i] = 0.0;
      }

      // deal with the case of two points (just do linear interpolation between them)
      if (nPoints ==2){
          Par[0][0] = (YResp[1] - YResp[0]) / ((XPos[1] - XPos[0]) * (XPos[1] - XPos[0]));
          Par[0][1] = 0.0;
          Par[0][2] = 0.0;
          // extra "virtual" segment at end to make Par array shape fit with knot arrays shapes
          Par[1][1] = 0.0;
          Par[1][2] = 0.0;

          return;
      } // if nPoints !=2 do full monotonic spline treatment:

      // first pass over knots to calculate the secants
      for (int i = 0; i < nPoints-1; ++i) {
        Secants[i] = (YResp[i+1] - YResp[i]) / (XPos[i+1] - XPos[i]);
        //std::cout<<"secant "<<i<<": "<<Secants[i]<<std::endl;
      }

      Tangents[0] = Secants[0];
      Tangents[nPoints-1] = Secants[nPoints -2];

      __float__ alpha;
      __float__ beta;

      // second pass over knots to calculate tangents
      for (int i = 1; i < nPoints-1; ++i) {
        if ((Secants[i-1] >= 0.0 && Secants[i] >= 0.0) | (Secants[i-1] < 0.0 && Secants[i] < 0.0)){ //check for same sign
          Tangents[i] = (Secants[i-1] + Secants[i]) /2.0;
        }
      }

      // third pass over knots to rescale tangents
      for (int i = 1; i < nPoints-1; ++i) {
        if (Secants[i] == 0.0){
          Tangents[i] = 0.0;
          Tangents[i+1] = 0.0;
        }

        else{
          alpha = Tangents[i]  / Secants[i];
          beta = Tangents[i+1] / Secants[i];

          if (alpha <0.0){
            Tangents[i] = 0.0;
          }
          if (beta < 0.0){
            Tangents[i+1] = 0.0;
          }

          if (alpha * alpha + beta * beta >9.0){
            __float__ tau = 3.0 / sqrt(alpha * alpha + beta * beta);
            Tangents[i]   = tau * alpha * Secants[i];
            Tangents[i+1] = tau * beta  * Secants[i];
          }
        }
        //std::cout<<"alpha, beta : "<<alpha<<", "<<beta<<std::endl;
        //std::cout<<"tangent "<<i<<": "<<Tangents[i]<<std::endl;

      } // finished rescaling tangents

      //std::cout<<"tangent "<<nPoints-1<<": "<<Tangents[nPoints-1]<<std::endl;

      // fourth pass over knots to calculate the coefficients for the spline
      __float__ dx;
      for(int i = 0; i <nPoints-1; i++){
        double b, c, d = -999.999;
        dx = XPos[i+1] - XPos[i];

        b = Tangents[i] * dx;
        c = 3.0* (YResp[i+1] - YResp[i]) -2.0 *dx * Tangents[i] - dx * Tangents[i +1];
        d = 2.0* (YResp[i] - YResp[i+1]) + dx * (Tangents[i] + Tangents[i+1]);

        Par[i][0] = b /  dx;
        Par[i][1] = c / (dx * dx);
        Par[i][2] = d / (dx * dx * dx);

        if((Par[i][0] == -999) | (Par[i][1] == -999) | (Par[i][2] ==-999) | (Par[i][0] == -999.999) | (Par[i][1] == -999.999) | (Par[i][2] ==-999.999)){
            std::cout<<"bad spline parameters for segment "<<i<<", will cause problems with GPU: (b, c, d) = "<<Par[i][0]<<", "<<Par[i][1]<<", "<<Par[i][2]<<std::endl;
        }
        //std::cout<<"b : "<<b<<std::endl;
        //std::cout<<"dx: "<<dx<<", x_0: "<<XPos[i]<<", x_1: "<<XPos[i+1]<<std::endl;
        //std::cout<<"    "<<" , y_0: "<<YResp[i]<<", y_1: "<<YResp[i+1]<<std::endl;
      }

      // include params for final "segment" outside of the spline so that par array fits with x and y arrays,
      // should never actually get used but if not set then the GPU code gets very angry
      Par[nPoints-1][0] = 0.0;
      Par[nPoints-1][1] = 0.0;
      Par[nPoints-1][2] = 0.0;

      delete spline;
      spline = NULL;
    }
    // finished calculating coeffs

    protected: //changed to protected from private so can be accessed by derived classes
    // values of the secants at each point (for calculating monotone spline)
    __float__ *Secants;
    // values of the tangents at each point (for calculating monotone spline)
    __float__ *Tangents;
};


// ************************
// Truncated spline class
class Truncated_Spline: public TSpline3_red {
// ************************
// cubic spline which is flat (returns y_first or y_last) if x outside of knot range
  public:
    // Empty constructor
    Truncated_Spline()
    :TSpline3_red()
    {
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    Truncated_Spline(TSpline3* &spline, int Param = -1)
    :TSpline3_red(spline, Param)
    {
    }

    // Empty destructor
    ~Truncated_Spline()
    {
    }

    // Find the segment relevant to this variation in x
    // See root/hist/hist/src/TSpline3::FindX(double) or samplePDFND....::FindSplineSegment
    inline int FindX(double x) {
      // The segment we're interested in (klow in ROOT code)
      int segment = 0;
      int kHigh = nPoints-1;
      // If the variation is below the lowest saved spline point
      if (x <= XPos[0]){
        segment = -1;
        // If the variation is above the highest saved spline point
      } else if (x >= XPos[nPoints-1]) {
        segment = -2;
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
          if (x > XPos[kHalf]) {
            segment = kHalf;
            // Else move kHigh down
          } else {
            kHigh = kHalf;
          }
        } // End the while: we've now done our binary search
      } // End the else: we've now found our point
      if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;
      return segment;
    }

    // Evaluate the weight from a variation
    inline double Eval(double var) {
      // Get the segment for this variation
      int segment = FindX(var);
      // The get the coefficients for this variation
      double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;

      if(segment >=0){
        GetCoeff(segment, x, y, b, c, d);
      }

      // if var is outside of the defined range, set the coefficients to 0 so that Eval just returns the value at the end point of the spline
      else if(segment == -1){
        GetKnot(0, x, y);
        b = 0.0;
        c = 0.0;
        d = 0.0;
      }
      else if(segment == -2){
        GetKnot(nPoints-1, x, y);
        b = 0.0;
        c = 0.0;
        d = 0.0;
      }

      double dx = var - x;
      // Evaluate the third order polynomial
      double weight = y+dx*(b+dx*(c+d*dx));
      return weight;
    }
};




// ************************
// Truncated Akima Spline class
class Truncated_Akima_Spline :public Akima_Spline {
// ************************
// Spline with Akima Spline coefficients which is flat outside of the defined knots

  public:

    // Empty constructor
    Truncated_Akima_Spline()
    :Akima_Spline()
    {
    }

    // The constructor that takes a TSpline3 pointer and copies in to memory
    Truncated_Akima_Spline(TSpline3* &spline, int Param = -1)
    :Akima_Spline(spline, Param)
    {
    }

    // Empty destructor
    ~Truncated_Akima_Spline()
    {
    }


    // See root/hist/hist/src/TSpline3::FindX(double) or samplePDFND....::FindSplineSegment
    inline int FindX(double x) {
      // The segment we're interested in (klow in ROOT code)
      int segment = 0;
      int kHigh = nPoints-1;
      // If the variation is below the lowest saved spline point
      if (x <= XPos[0]){
        segment = -1;
        // If the variation is above the highest saved spline point
      } else if (x >= XPos[nPoints-1]) {
        segment = -2;
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
          if (x > XPos[kHalf]) {
            segment = kHalf;
            // Else move kHigh down
          } else {
            kHigh = kHalf;
          }
        } // End the while: we've now done our binary search
      } // End the else: we've now found our point
      if (segment >= nPoints-1 && nPoints > 1) segment = nPoints-2;
      return segment;
    }


    // Evaluate the weight from a variation
    inline double Eval(double var) {
      // Get the segment for this variation
      int segment = FindX(var);
      // The get the coefficients for this variation
      double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;

      if(segment >=0){
        GetCoeff(segment, x, y, b, c, d);
      }

      // if var is outside of the defined range, set the coefficients to 0 so that Eval just returns the value at the end point of the spline
      else if(segment == -1){
        GetKnot(0, x, y);
        b = 0.0;
        c = 0.0;
        d = 0.0;
      }
      else if(segment == -2){
        GetKnot(nPoints-1, x, y);
        b = 0.0;
        c = 0.0;
        d = 0.0;
      }

      double dx = var - x;
      // Evaluate the third order polynomial
      double weight = y+dx*(b+dx*(c+d*dx));
      return weight;
    }
};


// ***************************
// A handy namespace for variables extraction
namespace MaCh3Utils {
  // ***************************

  // ***************************
  // Return mass for given PDG
  double GetMassFromPDG(int PDG);
  // ***************************

  TString SKSampleName_toLatexString(TString String);

  // Neutrino direction
  extern const double ND280NuDir[3];

  extern const double SKNuDir[3];

  extern std::unordered_map<int,int> KnownDetIDsMap;
  extern int nKnownDetIDs;

} // end MaCh3Utils namespace

// *****************
// Enum to track the target material
enum TargetMat {
  // *****************
  kTarget_H  = 1,
  kTarget_C  = 12,
  kTarget_N  = 14,
  kTarget_O  = 16,
  kTarget_Al = 27,
  kTarget_Ti = 48,
  kTarget_Fe = 56,
  kTarget_Pb = 207
};

// *****************
// Converted the Target Mat to a string
inline std::string TargetMat_ToString(TargetMat i) {
  // *****************
  std::string name;

  switch(i) {
    case kTarget_H:
      name = "Hydrogen";
      break;
    case kTarget_C:
      name = "Carbon";
      break;
    case kTarget_N:
      name = "Nitrogen";
      break;
    case kTarget_O:
      name = "Oxygen";
      break;
    case kTarget_Al:
      name = "Aluminium";
      break;
    case kTarget_Ti:
      name = "Titanium";
      break;
    case kTarget_Fe:
      name = "Iron";
      break;
    case kTarget_Pb:
      name = "Lead";
      break;
    default:
      name = "TargetMat_Undefined";
      break;
  }

  return name;
}

// *****************
// Enum to track the incoming neutrino species
enum NuPDG {
  // *****************
  kNue = 12,
  kNumu = 14,
  kNutau = 16,
  kNutau_bar = -16,
  kNumu_bar = -14,
  kNue_bar = -12
};

// Make an enum of the test statistic that we're using
enum TestStatistic {
  kPoisson,
  kBarlowBeeston,
  kIceCube
};

// **************************************************
// Convert a LLH type to a string
inline std::string TestStatistic_ToString(TestStatistic i) {
// **************************************************
    std::string name = "";

    switch(i) {
        case kPoisson:
        name = "Poisson";
        break;
        case kBarlowBeeston:
        name = "BarlowBeeston";
        break;
        case kIceCube:
        name = "IceCube";
        break;
        default:
            std::cerr << "UNKNOWN LIKELHOOD SPECIFIED TO ND280!" << std::endl;
            std::cerr << "You gave test-statistic " << i << std::endl;
            std::cerr << __FILE__ << ":" << __LINE__ << std::endl;
            throw;
    }
    return name;
}

// *******************
// Class containig information about interaction modes
class MaCh3_Modes {
  // *******************
  public:
    //KS: Getter for Modes, each experiemnt has to define it
    virtual int GetNModes();
    virtual std::string Mode_ToString(int i);
    virtual std::string Mode_ToFancyString(int i);
};



// Helper function for calculating unbinned Integral of TH2Poly i.e including overflow
double OverflowIntegral(TH2Poly*);

// Helper function for calculating binned Integral of TH2Poly i.e not including overflow
double NoOverflowIntegral(TH2Poly*);

// Poly Projectors
TH1D* PolyProjectionX(TObject* poly, std::string TempName, std::vector<double> xbins, bool computeErrors = false);
TH1D* PolyProjectionY(TObject* poly, std::string TempName, std::vector<double> ybins, bool computeErrors = false);

// Helper to Normalise histograms
TH2Poly* NormalisePoly(TH2Poly* Histogram);

// Helper to scale th2poly analogous to th2d scale with option "width"
TH2Poly* PolyScaleWidth(TH2Poly *Histogram, double scale);
// Helper to calc integral of th2poly analogous to th2d integra; with option "width"
double PolyIntegralWidth(TH2Poly *Histogram);

inline double returnCherenkovThresholdMomentum(int PDG);
inline double CalculateQ2(double PLep, double PUpd, double EnuTrue, double InitialQ2 = 0.0);
inline double CalculateEnu(double PLep, double cosTheta, double EB, bool neutrino);

// Helper to check if files exist or not
inline std::string file_exists(std::string filename) {
  std::ifstream infile(filename.c_str());
  if (!infile.good()) {
    std::cerr << "*** ERROR ***" << std::endl;
    std::cerr << "File " << filename << " does not exist" << std::endl;
    std::cerr << "Please try again" << std::endl;
    std::cerr << "*************" << std::endl;
    throw;
  }

  return filename;
}

inline int getValue();
inline int parseLine(char* line);

#endif
