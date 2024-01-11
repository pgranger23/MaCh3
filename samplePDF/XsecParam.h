#include <string>
#include "Selection.h"

class XsecParam
{
protected:
    std::string _fName;
    double _fGenerated;
    double _fPreFitValue;
    double _fError;
    double _fLowBound;
    double _fUpBound;
    double _fIndivStepScale;
    bool _fFlatPrior;
    int _fDetID;
    std::vector<int> _fModes;
    std::string _fParamType;
    std::string _fDetString;
    std::vector<Selection> _fSelections;

public:
    XsecParam(/* args */);
    ~XsecParam();

    void SetName(std::string name){_fName = name;};
    void SetPreFitValue(double val){_fPreFitValue = val;};
    void SetError(double val){_fError = val;};
    void SetGenerated(double val){_fGenerated = val;};
    void SetBounds(double lb, double ub){_fLowBound = lb; _fUpBound = ub;};
    void SetIndivStepScale(double val){_fIndivStepScale = val;};
    void SetDetID(int id){_fDetID = id;};
    void SetFlatPrior(bool flat){_fFlatPrior = flat;};
    void SetModes(std::vector<int> modes){_fModes = modes;};
    void SetSelections(std::vector<Selection> selections){_fSelections = selections;};

};

class XsecNorm : XsecParam
{
private:
    /* data */
public:
    XsecNorm(/* args */);
    ~XsecNorm();
};

class XsecSpline : XsecParam
{
protected:
    std::string _fFDSplineName;
    // std::string _fNDSplineName;
public:
    XsecSpline(/* args */);
    ~XsecSpline();

    void SetFDSplineName(std::string name){_fFDSplineName = name;};
    // void SetNDSplineName(std::string name){_fNDSplineName = name;};
};

XsecSpline::XsecSpline(/* args */)
{
}

XsecSpline::~XsecSpline()
{
}


XsecNorm::XsecNorm(/* args */)
{
}

XsecNorm::~XsecNorm()
{
}


XsecParam::XsecParam(/* args */)
{
}

XsecParam::~XsecParam()
{
}
