#include <vector>
#include <map>
#include <string>

class Selection
{
private:
    //Making these private in the class so that noone tries to do hacky things...
    enum Selector{
        kBelow,
        kBetween,
        kAbove,
        kAmong
    };
    static const std::map<std::string, Selector> str_to_selector;

    int _cutVar;
    std::vector<double> _cutValues;
    std::string _selector_str;
    Selector _selector;

    Selector getSelector(std::string str);

public:
    Selection(int cutVar, std::vector<double> cutValues, std::string selection);
    ~Selection();
    int getVar();
    bool Passes(double value);
};
