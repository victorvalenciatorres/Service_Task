# Service_Task



Map Muon Chambers:





--------------------------------------------------------------------------------------

Input Files used: 

    - ClustersMCH.root
    - ClustersMCH_LHC22t.root
    - o2sim_geometry-aligned.json

--------------------------------------------------------------------------------------


MODIFICATIONS:

    Modification in src File:

        - map_mch.cxx

        PATH: alice/O2/Detectors/MUON/MCH/GlobalMapping/src/map_mch.cxx


    Modification in Include File:

        - svgWriter.h

        PATH:  alice/O2/Detectors/MUON/MCH/Contour/include/MCHContour/SVGWriter.h


    Modification in CMakeList:

        - CMakeLists.txt 

        PATH: alice/O2/Detectors/MUON/MCH/GlobalMapping/CMakeLists.txt


--------------------------------------------------------------------------------------


COMPILACION (commands):

    cd alice ==>    alienv enter O2/latest ninja/latest

    cd sw/BUILD/O2-latest/O2 ==>   cmake --build . 

EXECUTION (command):

    stage/bin/o2-mch-map_mch --hidepadchannels --hidepads --de 100 


--------------------------------------------------------------------------------------

OPEN OUTPUT FILES: 

    open CHAMBERS-1-NB.html CHAMBERS-2-NB.html CHAMBERS-3-NB.html CHAMBERS-4-NB.html CHAMBERS-5-NB.html CHAMBERS-6-NB.html CHAMBERS-7-NB.html CHAMBERS-8-NB.html CHAMBERS-9-NB.html CHAMBERS-10-NB.html

    open  CHAMBERS-1-B.html CHAMBERS-2-B.html CHAMBERS-3-B.html  CHAMBERS-4-B.html CHAMBERS-5-B.html CHAMBERS-6-B.html CHAMBERS-7-B.html CHAMBERS-8-B.html CHAMBERS-9-B.html  CHAMBERS-10-B.html

--------------------------------------------------------------------------------------
