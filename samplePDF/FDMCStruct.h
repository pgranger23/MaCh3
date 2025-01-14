//constructors are same for all three so put in here
struct fdmc_base {
  int nutype; // 2 = numu/signue | -2 = numub | 1 = nue | -1 = nueb           
  int oscnutype;    
  int nupdg;
  bool signal; // true if signue                                              
  int nEvents; // how many MC events are there                              
  std::string flavourName;
  int *isbound;
  int **Target; // target the interaction was on

  int SampleDetID;

  //THe x_var and y_vars that you're binning in
  double** x_var;
  double** y_var;
  double **rw_etru;

  // spline bins                                                               
  unsigned int *xvar_s_bin;
  unsigned int *yvar_s_bin;
  unsigned int *enu_s_bin;

  // xsec bins  
  std::list< int > *xsec_norms_bins;

  //DB Speedup bits
  double Unity;

  int* nxsec_norm_pointers;
  const double*** xsec_norm_pointers;

  int* nxsec_spline_pointers;
  const double*** xsec_spline_pointers;

  const double** skdet_pointer;
  const double* EScale_pointer;

  int* ntotal_weight_pointers;
  double*** total_weight_pointers;
  double* total_w;

  int* XBin;
  int* YBin;
  int* NomXBin;
  int* NomYBin;

  bool *isNC;

  // histo pdf bins
  double *rw_lower_xbinedge; // lower to check if Eb has moved the erec bin
  double *rw_lower_lower_xbinedge; // lower to check if Eb has moved the erec bin
  double *rw_upper_xbinedge; // upper to check if Eb has moved the erec bin
  double *rw_upper_upper_xbinedge; // upper to check if Eb has moved the erec bin

  int **mode;

  double *osc_w; // oscillation weight                                        
  double *flux_w; // not the same as beam systematics weight!                 
  double *xsec_w;

  splineFDBase *splineFile; 
};
