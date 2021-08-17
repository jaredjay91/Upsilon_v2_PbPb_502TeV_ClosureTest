#!/bin/bash

whichv2=$1
date=2021_08_16

collId=kAADATA
ptLow=0.0
ptHigh=30.0
yLow=0.0
yHigh=2.4
cLow=0
cHigh=100
muPtCut=3.5

echo $date

#root -b -q SkimMCTree_recenter_Jared.C
#root -b -q SkimMCTree_Jared.C
#.root -b -q SkimMCTree_flatten_Jared.C

mkdir ExtractedYields
mkdir Plots

#mv ExtractedYields_2020_06_17_3_v2_$whichv2 ExtractedYields
#mv Plots_2020_06_17_3_v2_$whichv2 Plots

#<<COMMENT
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, 0, 3, $yLow, $yHigh, $cLow, $cHigh, $muPtCut, $whichv2)"
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, 3, 6, $yLow, $yHigh, $cLow, $cHigh, $muPtCut, $whichv2)"
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, 6, 30, $yLow, $yHigh, $cLow, $cHigh, $muPtCut, $whichv2)"
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 0, 10, $muPtCut, $whichv2)"
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 10, 30, $muPtCut, $whichv2)"
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 30, 50, $muPtCut, $whichv2)"
root -b -q -l "PlotMCmassAnddphiEp2.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 50, 100, $muPtCut, $whichv2)"


root -b -q -l "GetYieldsvsPhi.C($collId, 0, 3, $yLow, $yHigh, $cLow, $cHigh, $muPtCut, $whichv2)"
root -b -q -l "GetYieldsvsPhi.C($collId, 3, 6, $yLow, $yHigh, $cLow, $cHigh, $muPtCut, $whichv2)"
root -b -q -l "GetYieldsvsPhi.C($collId, 6, 30, $yLow, $yHigh, $cLow, $cHigh, $muPtCut, $whichv2)"
root -b -q -l "GetYieldsvsPhi.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 0, 10, $muPtCut, $whichv2)"
root -b -q -l "GetYieldsvsPhi.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 10, 30, $muPtCut, $whichv2)"
root -b -q -l "GetYieldsvsPhi.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 30, 50, $muPtCut, $whichv2)"
root -b -q -l "GetYieldsvsPhi.C($collId, $ptLow, $ptHigh, $yLow, $yHigh, 50, 100, $muPtCut, $whichv2)"
#COMMENT

root -b -q -l Get_v2_vs_var.C
root -b -q -l "draw_v2_pt.C(1,$whichv2)"
root -b -q -l "draw_v2_centrality.C(1,$whichv2)"

mv ExtractedYields ExtractedYields_${date}_v2_$whichv2
mv Plots Plots_${date}_v2_$whichv2

