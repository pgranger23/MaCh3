#include "styleManager.h"

namespace MaCh3Plotting{
    StyleManager::StyleManager(std::string styleTomlName){
        style_toml = toml_h::parse_card(styleTomlName);
    }

    std::string StyleManager::prettifyParamName(std::string origName) const {
        std::string prettyName = origName;

        auto const &PrettyNames = toml_h::find(style_toml, "PrettyNames");

        if(toml_h::contains(PrettyNames, origName)){
            prettyName = toml_h::find<std::string>(PrettyNames, origName);
        }
        
        return prettyName;
    }

    void StyleManager::setPalette(int rootPlotStyle) const{
        // set the colour palette to one of the root palettes
        gStyle->SetPalette(rootPlotStyle);
    }

    void StyleManager::setPalette(std::string configStyleName) const{
        // set the colour palette to one of the palettes defined in PlottingConfig.toml

        // get the definition of the provided style from the toml
        auto const &card_toml = toml_h::parse_card(std::string(std::getenv("MACH3")) + "/plotting/StyleConfig.toml");
        auto const &palettes  = toml_h::find(card_toml, "ColorPallettes");
        std::vector <std::vector <double>> paletteDef = toml_h::find <std::vector <std::vector <double>>> (palettes, configStyleName);

        const Int_t NCont = (Int_t)(paletteDef[0][0]);

        std::vector <double> stopVec   = paletteDef[1];
        std::vector <double> redsVec   = paletteDef[2];
        std::vector <double> greensVec = paletteDef[3];
        std::vector <double> bluesVec  = paletteDef[4];
        
        // get the number of colour stops and check all vectors are same size
        const Int_t NRGBs = stopVec.size();
        if((Int_t)redsVec.size()   != NRGBs ||
           (Int_t)greensVec.size() != NRGBs ||
           (Int_t)bluesVec.size()  != NRGBs)
        {
            std::cerr << "ERROR: invalid colour palettet defined in config toml: " << configStyleName << std::endl;
            std::cerr << "       RGB arrays dont all have the same size, please fix that" << std::endl;
            throw;
        }

        // root only likes arrays so convert to those
        Double_t stops[NRGBs];
        Double_t red[NRGBs];
        Double_t green[NRGBs];
        Double_t blue[NRGBs];
        for(int i = 0; i < NRGBs; i++){
            stops[i] = stopVec[i];
            red[i]   = redsVec[i];
            green[i] = greensVec[i];
            blue[i]  = bluesVec[i];
        }
        
        // now actually set the palette
        TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
        gStyle->SetNumberContours(NCont);
    }

    void StyleManager::setTH1Style(TH1* hist, std::string styleName) const{
        // get the definition of the provided style from the toml
        auto const &card_toml = toml_h::parse_card(std::string(std::getenv("MACH3")) + "/plotting/StyleConfig.toml");
        auto const &TH1Styles = toml_h::find(card_toml, "TH1Styles");
        auto const &styleDef  = toml_h::find(TH1Styles, styleName);

        if(toml_h::contains(styleDef, "MarkerColor")){
            hist->SetMarkerColor(toml_h::find<int>(styleDef, "MarkerColor"));
        }
        if(toml_h::contains(styleDef, "MarkerStyle")){
            hist->SetMarkerStyle(toml_h::find<int>(styleDef, "MarkerStyle"));
        }
        if(toml_h::contains(styleDef, "FillColor")){
            hist->SetFillColor(toml_h::find<int>(styleDef, "FillColor"));
        }
        if(toml_h::contains(styleDef, "FillStyle")){
            hist->SetFillStyle(toml_h::find<int>(styleDef, "FillStyle"));
        }
        if(toml_h::contains(styleDef, "LineColor")){
            hist->SetLineColor(toml_h::find<int>(styleDef, "LineColor"));
        }
        if(toml_h::contains(styleDef, "LineStyle")){
            hist->SetLineStyle(toml_h::find<int>(styleDef, "LineStyle"));
        }
    }

    std::string StyleManager::prettifySampleName(std::string fullName) const{
        std::string newName = fullName;

        std::string numu = "#nu_{#mu}";
        std::string numuBar = "#bar{#nu}_{#mu}";
        std::string pi = "#pi";
        std::string gamma = "#gamma";
        std::string nuBar = "#bar{#nu}";
        std::string boldP = "#bf{P}";

        if(newName.find("anti-numu") <= newName.length())newName.replace(newName.find("anti-numu"), 9, numuBar);
        if(newName.find("numu")      <= newName.length())newName.replace(newName.find("numu"), 4, numu);
        if(newName.find("NuMu")      <= newName.length())newName.replace(newName.find("NuMu"), 4, numu);
        if(newName.find("photon")    <= newName.length())newName.replace(newName.find("photon"), 6, gamma);
        if(newName.find("AntiNu")    <= newName.length())newName.replace(newName.find("AntiNu"), 6, nuBar);
        if(newName.find("protons")   <= newName.length())newName.replace(newName.find("protons"), 7, boldP);

        if(newName.find("Mode")      <= newName.length())newName.erase(newName.find("Mode"), 4);

        return newName;
    }


}