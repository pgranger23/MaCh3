#include "Selection.h"
#include <iostream>
#include <algorithm>
#include <cmath>

const std::map<std::string, Selection::Selector> Selection::str_to_selector = {
        {"BELOW", kBelow}, //Including value
        {"BETWEEN", kBetween}, //Excluding right of range
        {"ABOVE", kAbove}, //Excluding value
        {"AMONG", kAmong}
    };

Selection::Selection(std::string cutVar, std::vector<double> cutValues, std::string selection)
: _cutVar(cutVar), _cutValues(cutValues), _selector_str(selection)
{
    _selector = getSelector(_selector_str);

    //Checking that the size of the cutValues matches the selector used.
    //Safeguards user against some wrong uses.
    if(_selector == kBelow || _selector == kAbove || _selector){
        if(_cutValues.size() != 1){
            std::cerr << "Selector " << _selector_str << " requires 1 cut value. You provided " << _cutValues.size() <<  ". Aborting!" << std::endl;
            abort();
        }
    }
    if(_selector == kBetween){
        if (_cutValues.size() != 2){
            std::cerr << "Selector " << _selector_str << " requires 2 cut values. You provided " << _cutValues.size() <<  ". Aborting!" << std::endl;
            abort();
        }
    }
    if(_selector == kAmong){
        if (_cutValues.size() == 0){
            std::cerr << "Selector " << _selector_str << " requires values. You provided none. Aborting!" << std::endl;
            abort();
        }
    }
}

Selection::~Selection()
{
}

Selection::Selector Selection::getSelector(std::string str){
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);

    if(str_to_selector.count(str) != 0){
        return str_to_selector.at(str);
    }

    std::cerr << "The provided selection rule: '" << str <<"' does not match any of the available selectors: [ ";
    for(const auto &[key, value]: str_to_selector){
        std::cerr << key << " ";
    }
    std::cerr << "]. Aborting." << std::endl;
    abort();
}

std::string Selection::getVar() const{
    return _cutVar;
}
bool Selection::Passes(double value) const{
    if(_selector == kAmong){ //We assume that the selection is made amongs integers as double equality is not so obvious
        int intVal = std::round(value);
        if(std::find_if(_cutValues.begin(), _cutValues.end(), [intVal](double const &v){return std::round(v) == intVal;}) != _cutValues.end()){
            return true;
        }
        return false;
    }
    else if(_selector == kBelow)
    {
        return value <= _cutValues[0];
    }
    else if(_selector == kAbove)
    {
        return value > _cutValues[0];
    }
    else if(_selector == kBetween){
        return value > _cutValues[0] && value <= _cutValues[1];
    }
    else{
        std::cerr << "Looks like someone did not fully implement an additional selector type... Aborting." << std::endl;
        abort();
    }
    
}