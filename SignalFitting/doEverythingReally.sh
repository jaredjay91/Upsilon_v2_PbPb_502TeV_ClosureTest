#!/bin/bash

./doEverything.sh 0.05
./doEverything.sh 0.1
./doEverything.sh 0.2
./doEverything.sh 0.5

collId=kAADATA
ptLow=0.0
ptHigh=30.0
yLow=0.0
yHigh=2.4
cLow=0
cHigh=100
muPtCut=3.5

#<<COMMENT
root -b -q -l "PlotMCmassAndEpAngles.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, $cLow, $cHigh, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, 0, 3, $yLow, $yHigh, $cLow, $cHigh, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, 3, 6, $yLow, $yHigh, $cLow, $cHigh, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, 6, 30, $yLow, $yHigh, $cLow, $cHigh, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 0, 10, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 10, 30, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 30, 50, $muPtCut)"
root -b -q -l "PlotMCmassAndEpAngles.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 50, 100, $muPtCut)"

root -b -q -l "drawResCor.C()"
#COMMENT
