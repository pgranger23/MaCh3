#include "Structs.h"

// ********************************************
// Struct for ND280 MC Event
class NDEVENT {
  // ********************************************
  public:
    // Virtual functions
    NDEVENT();
    ~NDEVENT();
    void Print();

    // Momentum of event
    __float__ mom;
    // Eb adjusted Momentum of event
    __float__ mom_adj;
    // Flux weight of event
    __float__ theta;
    // Psyche event sample does event belong to
    __int__ event_sample;
    // Detector bin does event belong to
    __int__ det_bin;
    // Fitting bin does an event belong to
    __int__ hist_bin;
    // Flux weight of event
    __float__ flux_w;

    // number of pointers to norm bin
    __int__  nxsec_norm_pointers;
    // pointer to norm bin
    std::vector < const double* > xsec_norm_pointers;

    //Xsec based variables
    // Mode
    __int__ mode;
    // One time weight to apply to event, e.g. NC1gamma
    __float__ weight;
};


// ********************************************
// Struct for ND280 MC Event wiht AUXILIARY variables, which we only need before starting fit
class NDEVENT_AUXILIARY {
  // ********************************************
  public:
    // Virtual functions
    NDEVENT_AUXILIARY();
    ~NDEVENT_AUXILIARY();
    void Print();

    // is RHC or FHC
    bool isRHC;
    // Which species
    __int__ species;
    // Which target the event is on
    __int__ target;
    // Enu of event
    __float__ Enu;
    // Q2 of event
    __float__ Q2;

    // norm bin used by normalisation dials
    std::vector < __int__ >  xsec_norms_bins;
};

