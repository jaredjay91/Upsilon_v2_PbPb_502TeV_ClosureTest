#!/bin/bash

nevt=-1
todaysDate=20200626
flattenBinByBin=kTRUE

#root -l -q "SkimMCTree_raw_GetAverageQ.C($nevt)"
#root -l -q "SkimMCTree_recentering_GetAverageEP.C($nevt)"
root -l -q "SkimMCTree_flatten_weight_GetResCor.C($nevt,$todaysDate,$flattenBinByBin,0.5)"
root -l -q "SkimMCTree_flatten_weight_GetResCor.C($nevt,$todaysDate,$flattenBinByBin,0.2)"
root -l -q "SkimMCTree_flatten_weight_GetResCor.C($nevt,$todaysDate,$flattenBinByBin,0.1)"
root -l -q "SkimMCTree_flatten_weight_GetResCor.C($nevt,$todaysDate,$flattenBinByBin,0.05)"
