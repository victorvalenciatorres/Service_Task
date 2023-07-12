// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.


///
/// @author  Victor Valencia 


#include "boost/program_options.hpp"
#include "MCHMappingInterface/CathodeSegmentation.h"
#include "MCHMappingInterface/Segmentation.h"
#include "MCHMappingSegContour/CathodeSegmentationContours.h"
#include "MCHMappingSegContour/CathodeSegmentationSVGWriter.h"
#include "MCHGeometryTransformer/Transformations.h"
#include "TGeoManager.h"
#include "MCHContour/SVGWriter.h"
#include "MCHConstants/DetectionElements.h"
#include "MCHGlobalMapping/ChannelCode.h"
#include "MCHGlobalMapping/DsIndex.h"
#include <fstream>
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TH1F.h"
#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "DetectorsBase/GeometryManager.h"
#include "MCHContour/Polygon.h"


using namespace o2::mch::mapping;

namespace po = boost::program_options;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// Function with 156 translation offsets to not overlap the deIds of each chamber
std::pair<double, double> getTranslationOffset(int deId) {

    std::map<int, std::pair<double, double>> translationOffsets = {

        //1st type of Chamber Shape (chambers 1,2,3,4)
        {100, {5.0, 0.0}}, {101, {-5.0, 0.0}}, {102, {-5.0, 10.0}}, {103, {5.0, 10.0}},
        {200, {5.0, 0.0}}, {201, {-5.0, 0.0}}, {202, {-5.0, 10.0}}, {203, {5.0, 10.0}},
        {300, {5.0, 0.0}}, {301, {-5.0, 0.0}}, {302, {-5.0, 10.0}}, {303, {5.0, 10.0}},
        {400, {5.0, 0.0}}, {401, {-5.0, 0.0}}, {402, {-5.0, 10.0}}, {403, {5.0, 10.0}},
        //2nd type of Chamber Shape (chambers 5,6)
        {500, {5.0, 5.0}}, {501, {5.0, -10.0}}, {502, {5.0, -20.0}}, {503, {5.0, -30.0}},
        {504, {5.0, -40.0}}, {505, {0, -40.0}}, {506, {0, -30.0}}, {507, {0, -20.0}}, {508, {0, -10.0}},
        {509, {0.0, 5.0}}, {510, {0.0, 20.0}}, {511, {0.0, 30.0}}, {512, {0.0, 40.0}},
        {513, {0.0, 50.0}}, {514, {5.0, 50.0}}, {515, {5.0, 40.0}}, {516, {5.0, 30.0}}, {517, {5.0, 20.0}},
        {600, {5.0, 5.0}}, {601, {5.0, -10.0}}, {602, {5.0, -20.0}}, {603, {5.0, -30.0}},
        {604, {5.0, -40.0}}, {605, {0, -40.0}}, {606, {0, -30.0}}, {607, {0, -20.0}}, {608, {0, -10.0}},
        {609, {0.0, 5.0}}, {610, {0.0, 20.0}}, {611, {0.0, 30.0}}, {612, {0.0, 40.0}},
        {613, {0.0, 50.0}}, {614, {5.0, 50.0}}, {615, {5.0, 40.0}}, {616, {5.0, 30.0}}, {617, {5.0, 20.0}},

        //3rd type of Chamber Shape (chambers 7,8,9,10)
        {700, {3.0, 5.0}}, {701, {3.0, -10.0}}, {702, {3.0, -20.0}}, {703, {3.0, -30.0}}, {704, {3.0, -50.0}}, {705, {3, -60.0}}, {706, {3, -80.0}}, {707, {0, -80.0}}, {708, {0, -60.0}}, {709, {0.0, -50}}, {710, {0.0, -30.0}}, {711, {0.0, -20.0}}, {712, {0.0, -10.0}}, {713, {0.0, 5.0}}, {714, {0.0, 20.0}}, {715, {0.0, 30.0}}, {716, {0.0, 42.0}}, {717, {0.0, 64.0}}, {718, {0.0, 79.0}}, {719, {0.0, 101.0}}, {720, {3.0, 101.0}}, {721, {3.0, 79}}, {722, {3.0, 64.0}}, {723, {3.0, 42.0}}, {724, {3.0, 30}}, {725, {3.0, 20}},
        {800, {3.0, 5.0}}, {801, {3.0, -10.0}}, {802, {3.0, -20.0}}, {803, {3.0, -30.0}}, {804, {3.0, -50.0}}, {805, {3, -60.0}}, {806, {3, -80.0}}, {807, {0, -80.0}}, {808, {0, -60.0}}, {809, {0.0, -50}}, {810, {0.0, -30.0}}, {811, {0.0, -20.0}}, {812, {0.0, -10.0}}, {813, {0.0, 5.0}}, {814, {0.0, 20.0}}, {815, {0.0, 30.0}}, {816, {0.0, 42.0}}, {817, {0.0, 64.0}}, {818, {0.0, 79.0}}, {819, {0.0, 101.0}}, {820, {3.0, 101.0}}, {821, {3.0, 79}}, {822, {3.0, 64.0}}, {823, {3.0, 42.0}}, {824, {3.0, 30}}, {825, {3.0, 20}},  // deId 825 has an extra reflexion, this is why the value is adjusted like this..
        {900, {0.0, 5.0}}, {901, {0.0, -10.0}}, {902, {0.0, -20.0}}, {903, {0.0, -30.0}}, {904, {0.0, -50.0}}, {905, {0, -60.0}}, {906, {0, -80.0}}, {907, {0, -80.0}}, {908, {0, -60.0}}, {909, {0.0, -50}}, {910, {0.0, -30.0}}, {911, {0.0, -20.0}}, {912, {0.0, -10.0}}, {913, {0.0, 5.0}}, {914, {0.0, 20.0}}, {915, {0.0, 30.0}}, {916, {0.0, 50.0}}, {917, {0.0, 70.0}}, {918, {0.0, 90.0}}, {919, {0.0, 110.0}}, {920, {0.0, 110.0}}, {921, {0.0, 90}}, {922, {0.0, 70}}, {923, {0.0, 50}}, {924, {0.0, 30}}, {925, {0.0, 20}},
        {1000, {0.0, 5.0}}, {1001, {0.0, -10.0}}, {1002, {0.0, -20.0}}, {1003, {0.0, -30.0}}, {1004, {0.0, -50.0}}, {1005, {0, -60.0}}, {1006, {0, -80.0}}, {1007, {0, -80.0}}, {1008, {0, -60.0}}, {1009, {0.0, -50}}, {1010, {0.0, -30.0}}, {1011, {0.0, -20.0}}, {1012, {0.0, -10.0}}, {1013, {0.0, 5.0}}, {1014, {0.0, 20.0}}, {1015, {0.0, 30.0}}, {1016, {0.0, 50.0}}, {1017, {0.0, 70.0}}, {1018, {0.0, 90.0}}, {1019, {0.0, 110.0}}, {1020, {0.0, 110.0}}, {1021, {0.0, 90}}, {1022, {0.0, 70}}, {1023, {0.0, 50}}, {1024, {0.0, 30}}, {1025, {0.0, 20}}
       
    };
    
        return translationOffsets[deId];
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// Function to transform local contours to global contours coordinates 
std::vector<o2::mch::contour::Contour<double>> transformLocalToGlobal(int deId, bool bending, const o2::mch::geo::TransformationCreator& transformation) {
    CathodeSegmentation cSeg{deId, bending}; 
    std::vector<o2::mch::contour::Contour<double>> dualSampaContoursIn = getDualSampaContours(cSeg);
    std::vector<o2::mch::contour::Contour<double>> dualSampaContoursOut; 
    auto transformDeId = transformation(deId);

    for(int i = 0; i < dualSampaContoursIn.size(); i++) {
        auto dsContourIn = dualSampaContoursIn[i];
        o2::mch::contour::Contour<double> dsContourOut;

        for (auto p = 0; p < dsContourIn.size(); p++) {
            auto polygIn = dsContourIn[p];
            std::vector<o2::mch::contour::Vertex<double>> verticesOut;

            for (auto v = 0; v < polygIn.size(); v++) {
                auto vertexIn = polygIn[v];
                o2::math_utils::Point3D<float> lpos(vertexIn.x, vertexIn.y, 0.0); 
                o2::math_utils::Point3D<float> gpos;        

                // Transformations: Rotations + Translations (local --> global) 
                transformDeId.LocalToMaster(lpos, gpos);
                o2::mch::contour::Vertex<double> vertexOut;

                auto offset = getTranslationOffset(deId);
                   
            
                    vertexOut.y = -gpos.Y();       // reflection anti-clockwise (minus value)
                    vertexOut.y += offset.second;
                    vertexOut.x = gpos.X() + offset.first; 
             
                verticesOut.push_back(vertexOut);
        
            }

            o2::mch::contour::Polygon<double> polygOut(verticesOut.begin(), verticesOut.end()); 
            dsContourOut.addPolygon(polygOut);
        }

        dualSampaContoursOut.push_back(dsContourOut);
    }


    return dualSampaContoursOut;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// GET ALL deID of a given Chamber
std::vector<int> getAllDeIds(int nChamber) {
    std::vector<int> deIds;
    int i_1 = *std::move(o2::mch::constants::deId2DeIndex(nChamber*100));       //Using move to convert optional integer into normal integer
    int i_2 = *std::move(o2::mch::constants::deId2DeIndex((nChamber+1)*100));
    int Chamberlength = i_2 -i_1;

    if(nChamber<10){
        for (auto i = 0; i < Chamberlength; i++) {                              //Getting deIds for chambers from 1 - 9
            deIds.push_back(o2::mch::constants::deIdsForAllMCH[i+i_1]);
        }
    }
    else{
        for (int i = 0; i < 26; i++) {
            deIds.push_back(o2::mch::constants::deIdsForAllMCH[i+i_1]);         //Getting deIds for chamber 10
        }
    }
    return deIds;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// GET ALL DualSampa of a given deId
std::vector<int> getDualSampas(int deId) {
    std::vector<int> dualSampas;
    const o2::mch::mapping::Segmentation& seg = o2::mch::mapping::segmentation(deId);
    seg.forEachDualSampa([&dualSampas](int ds) { dualSampas.push_back(ds); });
    return dualSampas;  
}


// GET ALL Bending or Non-Bending DualSampa of a given deId:
std::vector<int> getDualSampasBorNB(int deId, bool isBending) {
    std::vector<int> dualSampas;
    const o2::mch::mapping::Segmentation& seg = o2::mch::mapping::segmentation(deId);
    seg.forEachDualSampa([&dualSampas, isBending, &seg](int dualSampaId) {
        bool foundDualSampa = false;
        seg.forEachPad([&foundDualSampa, isBending, dualSampaId, &seg](int dePadIndex) {
            if (seg.isBendingPad(dePadIndex) == isBending && seg.padDualSampaId(dePadIndex) == dualSampaId) {
                foundDualSampa = true;
                return;  
            }
        });
        if (foundDualSampa) {
            dualSampas.push_back(dualSampaId);
        }
    });
    
  
    return dualSampas; 
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// Get All Dual Sampas of n Chambers:
std::vector<int> getAllDualSampas(int numChambers) {
    std::vector<int> allDualSampas;
    int count = 0; 

    for (int nChamber = 1; nChamber <= numChambers; nChamber++) {
        std::vector<int> deIds = getAllDeIds(nChamber);
        std::vector<int> dualSampas;
        for (auto deId : deIds) {
            std::vector<int> dsIds = getDualSampas(deId);
            for (auto dsId : dsIds) {
                dsId += count; 
                dualSampas.push_back(dsId);
            }
        }
        count = dualSampas.back() + 1; 
        allDualSampas.insert(allDualSampas.end(), dualSampas.begin(), dualSampas.end());
    }

    std::cout << std::endl;

    return allDualSampas;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

//Convert DsIndex (Global) to DsId (Local)
uint16_t convertDsIndextoDsId(o2::mch::DsIndex dsIndex)
{
    o2::mch::raw::DsDetId dsDetId = o2::mch::getDsDetId(dsIndex);
    uint16_t dsId = dsDetId.dsId();
    return dsId;
}


//Convert from DsId & deId (Local) to DsIndex (Global)
uint16_t getDsIndexFromDsIdAndDeId(uint16_t dsId, uint16_t deId)
{
  o2::mch::raw::DsDetId dsDetId(deId, dsId); // Create a DsDetId object with the given deId and dsId
  return o2::mch::getDsIndex(dsDetId); // Get the corresponding dsIndex from the getDsIndex function
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

//Getting ClustersPerDualSampa TH1F Histogram stored in root file
TH1F* getrootHistogram1() {

    TFile* file = TFile::Open("/Users/valencia/test/ClustersMCH_LHC22t.root");
    TH1F* ClustersperDualSampa = (TH1F*)file->Get("ClustersPerDualSampa");
    return ClustersperDualSampa;
}

TH1F* getrootHistogram2() {

    TFile* file = TFile::Open("/Users/valencia/Desktop/emilie/clusters1.root");
    TH1F* ClustersperDualSampa = (TH1F*)file->Get("Clusters");
    return ClustersperDualSampa;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

//Storing clusters and dsindex from TH1F histogram
std::pair<std::vector<int>, std::vector<uint16_t>> processClustersperDualSampa(const TH1F* ClustersperDualSampa) {

    std::vector<int> nClusters;
    std::vector<uint16_t> dsindex;
    std::vector<uint16_t> dsId;

    for (int i = 1; i <= ClustersperDualSampa->GetNbinsX(); i++) {
        int clusters = ClustersperDualSampa->GetBinContent(i);
        nClusters.push_back(clusters);
        dsindex.push_back(i - 1);
        dsId.push_back(convertDsIndextoDsId(i - 1));
        o2::mch::raw::DsDetId dsDetId = o2::mch::getDsDetId(i - 1);
        uint16_t currentDsId = dsDetId.dsId();
        uint16_t deId = dsDetId.deId();
    }

    return std::make_pair(nClusters, dsindex);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

//Getting transformations from the aligned geometry 
o2::mch::geo::TransformationCreator loadGeometry(const std::string& name) {
    
    o2::base::GeometryManager::loadGeometry(name.c_str());

    auto transformation = o2::mch::geo::transformationFromTGeoManager(*gGeoManager);

    return transformation;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

 //Creation of Color Gradiant Palette
std::vector<std::string> colorGradiant()
{ 
   
   std::vector<std::string> hexcolors;
   gStyle->SetPalette(kSunset);
   TColor::InvertPalette();
   TArrayI colors = TColor::GetPalette();
   for (int i = 0; i < colors.GetSize(); i++) {
    Int_t col = colors.At(i);
    TColor *tcol = gROOT->GetColor(col);
    hexcolors.push_back(tcol->AsHexString());
   }

    TColor::InvertPalette();     

    return hexcolors;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

std::vector<int> gradient( int n, int m ) {
  std::div_t q { 0, 0 };
  std::vector<int> grad(m);
  for( int i=1 ; i<m ; ++i ) {
    q = std::div( n + q.rem, m-1 );
    grad[i] = grad[i-1] + q.quot;
  }
  return grad;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Add Rentangle with 255 colors + Numbers (from 0 to 1)
void addRectangleContour(int nChamber, o2::mch::contour::Contour<double>& contour, o2::mch::contour::SVGWriter& w) {
    double rectWidth;  
    double rectHeight;  
    double rectX;  
    double rectY;  

    double step = 1.0 / 10; 
    double textX; 
    double textY; 
    double startY; 
    double spacing; 

    std::string title = ""; 

    //size of rectangle depending on chamber
    if (nChamber == 1 || nChamber == 2) {
        rectWidth = 10;
        rectHeight = 200;
        rectX = 110;
        rectY = -96;

        textX = rectX + rectWidth + 5;
        startY = 82 + rectY + rectHeight / 2;
        spacing = 1 + rectHeight / 11;
        w.text(title, textX - 5, startY -190); 
    } else if (nChamber == 3 || nChamber == 4) {
        rectWidth = 20;
        rectHeight = 200;
        rectX = 130;
        rectY = -94;

        textX = rectX + rectWidth + 5;
        startY = 82 + rectY + rectHeight / 2;
        spacing = 1 + rectHeight / 11;
        w.text(title, textX - 7, startY -190); 
    } else if (nChamber == 5 || nChamber == 6) { 
        rectWidth = 20;
        rectHeight = 225;
        rectX = 175;
        rectY = -20;

        textX = rectX + rectWidth + 5;
        startY = 90 + rectY + rectHeight / 2;
        spacing = 0.7 + rectHeight / 11;
        w.text(title, textX - 5, startY -210); 
    } else if (nChamber == 7 || nChamber == 8) {
        rectWidth = 30;
        rectHeight = 325;
        rectX = 270;
        rectY = -20;

        textX = rectX + rectWidth + 5;
        startY = 130 + rectY + rectHeight / 2;
        spacing = 1.5 + rectHeight / 11;
        w.text(title, textX - 5, startY -300); 
    } else if (nChamber == 9 || nChamber == 10) {
        rectWidth = 35;
        rectHeight = 325;
        rectX = 270;
        rectY = -20;

        textX = rectX + rectWidth + 5;
        startY = 129 + rectY + rectHeight / 2;
        spacing = 1.5 + rectHeight / 11;
        w.text(title, textX - 5, startY -300); 
    }

    double stepX = rectWidth / 255.0;
    double stepY = rectHeight / 255.0;

    for (int i = 0; i < 255; i++) {
        double startY = rectY + i * stepY;
        double endY = rectY + (i + 1) * stepY;

        double startX = rectX;
        double endX = rectX + rectWidth;

        contour.addPolygon({
            {startX, startY},
            {endX, startY},
            {endX, endY},
            {startX, endY}
        });
    }

    //Produces Numbers from 0 to 1
    double n;
    for (int i = 10; i >= 0; i--) {
        n = i * step; 
        double textY = startY + (1 - i) * spacing; 

        int integerPart = static_cast<int>(n); 
        int decimalPart = static_cast<int>((n - integerPart) * 10); 

        std::string numberString = std::to_string(integerPart) + "." + std::to_string(decimalPart);
        w.text(numberString, textX, textY);
    }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

//Calculation of the maximum ratio of all Clusters/DSArea
double calculateMaxRatio(int nChamber, bool bending, const TH1F* ClustersperDualSampa, o2::mch::geo::TransformationCreator transformation) {

    auto nClusters_dsindex = processClustersperDualSampa(ClustersperDualSampa);
    std::vector<int> nClusters = nClusters_dsindex.first;
    std::vector<uint16_t> dsindex = nClusters_dsindex.second;

    // Getting All DeIds for all Chambers
    auto deIds = getAllDeIds(nChamber);

    std::vector<std::pair<double, int>> areas; // Vector with Areas of DS Contours paired with dsIndex
    double maxRatio = 0.0;
    double ratio = 0.0;
    int dsIndex;

    // Contours of all deId transformed + all dsContourOut
    for (auto deId : deIds) {
        auto dualSampaContoursOut = transformLocalToGlobal(deId, bending, transformation);
        // Get the map index to dsId
        auto dsIds = getDualSampasBorNB(deId, bending); 
        for (auto i = 0; i < dualSampaContoursOut.size(); i++) {
            auto& contour = dualSampaContoursOut[i]; // Get the current contour
            double area = 0.0; // Area of a DS Contour
            for (const auto& poly : contour.getPolygons()) {
                area += poly.signedArea();
            }
            // Get the local dsId
            auto dsId = dsIds[i];
            // Convert local dsId to global dsIndex (for a given deId)
            dsIndex = getDsIndexFromDsIdAndDeId(dsId, deId);
            // Create a pair of area and dsIndex and push it to the areas vector
            areas.push_back(std::make_pair(area, dsIndex));
           
            ratio = nClusters[dsIndex] / std::abs(areas.back().first);
                  
            maxRatio = std::max(maxRatio, ratio);

        }
    }
    

    return maxRatio;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// Creating Chambers in SVG format
void svgChamber(o2::mch::contour::SVGWriter& w, int nChamber, bool bending, const TH1F* ClustersperDualSampa, o2::mch::geo::TransformationCreator transformation, double maxRatio) {


    int nclustermax = ClustersperDualSampa->GetMaximum();

    //Load clusters from TH1F Histogram
    auto nClusters_dsindex = processClustersperDualSampa(ClustersperDualSampa);
    std::vector<int> nClusters = nClusters_dsindex.first;
    std::vector<uint16_t> dsindex = nClusters_dsindex.second;

    //Colors Vector in HEX RBG format
    std::vector<std::string> colors = colorGradiant();

    //Getting all deIds for all Chambers
    auto deIds = getAllDeIds(nChamber);

    double epsilon=1.e-6; //Small shift

    std::vector<std::pair<double, int>> areas; // Vector with Areas of DS Contours paired with dsIndex

     double maxratio = calculateMaxRatio( nChamber, bending, ClustersperDualSampa, transformation);  //Maximun Ratio

    // Contours of all deId transformated  + SVGWRITER of all dsContourOut bien
    for (auto deId : deIds) {
        auto dualSampaContoursOut = transformLocalToGlobal(deId, bending, transformation);
        w.svgGroupStart("dualsampas");
                std::string  str = ".dualsampas { fill:";
                str += colors[0]; // Assign the first color in the vector
                str += "; stroke-width: 0.25px; stroke: #333333;}";
                w.addStyle(str);
        // Get the map index to dsId
        auto dsIds = getDualSampasBorNB(deId, bending); 
        for (auto i = 0; i < dualSampaContoursOut.size(); i++) {

            auto& contour = dualSampaContoursOut[i]; // Get the current contour

            double area = 0.0; // Area of a DS Contour
            for (const auto& poly : contour.getPolygons()) {
                area += poly.signedArea();
            }

            auto dsId = dsIds[i];
            int dsIndex = getDsIndexFromDsIdAndDeId(dsId, deId);
            areas.push_back(std::make_pair(area, dsIndex));

            int colorId = int((nClusters[dsIndex] / std::abs(areas.back().first)) / (maxratio + epsilon) * colors.size());
            
            if (nClusters[dsIndex] == 0) {
                w.contour(dualSampaContoursOut[i], "#FFFFFF");      // White color for 0 Cluster  <--->  "#00FF00" for Bright green color
            } else {
                w.contour(dualSampaContoursOut[i], colors[colorId]);
            }

        }

        w.svgGroupEnd();
    }

   
    // Add rectangle for color scale + text:
    o2::mch::contour::Contour<double> rectangleContour;
    addRectangleContour(nChamber, rectangleContour, w);

    for(auto i =0; i<rectangleContour.size();i++){

        w.polygon(rectangleContour[i],colors[rectangleContour.size()-1-i]);  //rectangle white to black

    }
    

  
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

int main(int argc, char* argv[])
{

  std::string prefix;
  std::vector<int> detElemIds; 
  std::vector<int> dualSampas; 

  using Point = std::pair<double, double>;
  std::vector<std::string> pointStrings;
  std::vector<Point> points;
  po::variables_map vm;
  po::options_description generic("Generic options");

  generic.add_options()("help", "produce help message")("hidepads", "hide pad outlines")(
    "hidedualsampas", "hide dualsampa outlines")("hidedes", "hide detection element outline")(
    "hidepadchannels", "hide pad channel numbering")("de", po::value<std::vector<int>>(&detElemIds),
                                                     "which detection element to consider")(
    "prefix", po::value<std::string>(&prefix)->default_value("seg"), "prefix used for outfile filename(s)")(
    "point", po::value<std::vector<std::string>>(&pointStrings), "points to show")("all", "use all detection elements");

  po::options_description cmdline;
  cmdline.add(generic);

  po::store(po::command_line_parser(argc, argv).options(cmdline).run(), vm);
  po::notify(vm);

    std::cout << "deId:" << detElemIds[0] << "\n";

  if (vm.count("help")) {
    std::cout << generic << "\n";
    return 2;
  }

  if (vm.count("de") && vm.count("all")) {
    std::cout << "--all and --de options are mutually exclusive. --all will be used\n";
    detElemIds.clear();
  }

  if (vm.count("all")) {
    o2::mch::mapping::forOneDetectionElementOfEachSegmentationType(
      [&detElemIds](int detElemId) { detElemIds.push_back(detElemId); });
  }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

// Get All dualSampas for 10 chambers
getAllDualSampas(10);

// Define the bounding boxes for the 10 images:
std::vector<o2::mch::contour::BBox<double>> bboxes = {
    {-150, -150, 150, 150},
    {-150, -150, 150, 150},
    {-200, -200, 200, 200},
    {-200, -200, 200, 200},
    {-240, -240, 240, 240},
    {-240, -240, 240, 240},
    {-350, -350, 350, 350},
    {-450, -450, 450, 450},
    {-450, -450, 450, 450},
    {-450, -450, 450, 450}
};

//Load Aligned Geometry from file.root
std::string name = "o2sim_geometry-aligned.root";
    
// Create with SVGWriter the 10 Chambers with their respective bounding boxes
for (auto isBendingPlane : {true, false}) {
    for (int i = 0; i < 10; i++) {
        std::ofstream outv("CHAMBERS-" + std::to_string(i+1) + "-" +   //output file
                        (isBendingPlane ? "B" : "NB") + ".html");

        // Creating bboxes 
        o2::mch::contour::SVGWriter wSegLeft(bboxes[i]);
        o2::mch::contour::SVGWriter wSegRight(bboxes[i]);
       
        // Creating Left and Right Chambers  
        double maxRatioLeft = calculateMaxRatio(i+1, isBendingPlane, getrootHistogram1(), loadGeometry(name));
        double maxRatioRight = calculateMaxRatio(i+1, isBendingPlane, getrootHistogram2(), loadGeometry(name));
        svgChamber(wSegLeft, i+1, isBendingPlane, getrootHistogram1(), loadGeometry( name), maxRatioLeft);
        svgChamber(wSegRight, i+1, isBendingPlane, getrootHistogram2(), loadGeometry(name), maxRatioRight);
        

        // Write in HTML left and right chambers (using <div> tag)
        outv << "<div style='display:flex;justify-content:center'>" << std::endl;
        wSegLeft.writeHTML(outv);
        outv << "<div style='margin-left:20px;'></div>" << std::endl;
        wSegRight.writeHTML(outv);
        outv << "</div>" << std::endl;

      
    } 
}
  double realNmax = 254.4;
  int roundNmax = int(realNmax) +1;

  for( int i : gradient(roundNmax,11) ){

    std::cout << i << ' ';
     
  }
  
  std::cout << '\n';

  return 0;
}
