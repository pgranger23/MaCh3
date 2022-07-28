#include "NDMCStruct.h"

// ********************
// Constructor
NDEVENT::NDEVENT() {
// ********************
    mom = __BAD_DOUBLE__;
    mom_adj = __BAD_DOUBLE__;
    theta = __BAD_DOUBLE__;
    event_sample = __BAD_INT__;
    det_bin  = __BAD_INT__;
    hist_bin = __BAD_INT__;
    flux_w = __BAD_DOUBLE__;

    nxsec_norm_pointers = __BAD_INT__;
    mode = __BAD_INT__;
    weight = __BAD_DOUBLE__;
}

// ********************
// Empty destructor
NDEVENT::~NDEVENT() {
// ********************

}

// ********************
// Print a given event number
void NDEVENT::Print() {
// ********************
  std::cout << "*** Printing NDEVENT info:" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "pLepton: " << mom << " MeV" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "pLep adj: " << mom_adj << " MeV" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "CosThLepton: " << theta << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Sample: " << event_sample << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Detector bin: " << det_bin << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Histogram bin: " << hist_bin << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Constant Weight: " << flux_w << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Mode: " << mode << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "One time Xsec weight to apply: " << weight << std::endl;
}


// ********************
// Constructor
NDEVENT_AUXILIARY::NDEVENT_AUXILIARY() {
// ********************
    isRHC = false;
    target = __BAD_INT__;
    species = __BAD_INT__;
    Enu = __BAD_DOUBLE__;
    Q2 = __BAD_DOUBLE__;
}

// ********************
// Empty destructor
NDEVENT_AUXILIARY::~NDEVENT_AUXILIARY() {
// ********************

}

// ********************
// Print a given event number
void NDEVENT_AUXILIARY::Print() {
// ********************
  std::cout << "*** Printing NDEVENT_AUXILIARY info:" << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Is RHC: " << isRHC << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Target: " << target << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "species: " << species << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Enu: " << Enu << std::endl;
  std::cout << "    " << std::left << std::setw(20) << "Q2: " << Q2 << std::endl;
}
