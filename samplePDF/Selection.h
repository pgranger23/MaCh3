#include <vector>
#include <map>

class Selection
{
private:
    int _cutVar;
    std::vector<double> _cutValues;
    std::string _selector_str;
    Selector _selector;

    enum Selector{
        kBelow,
        kBetween,
        kAbove,
        kAmong
    };

    static const std::map<std::string, Selector> str_to_selector = {
        {"BELOW", kBelow}, //Including value
        {"BETWEEN", kBetween}, //Excluding right of range
        {"ABOVE", kAbove}, //Excluding value
        {"AMONG", kAmong}
    };

    Selector getSelector(std::string str);

public:
    Selection(int cutVar, std::vector<double> cutValues, std::string selection);
    ~Selection();
    int getVar();
    bool Passes(double value);
};
