#!/bin/bash

nevt=1000000
todaysDate=20200612
flattenBinByBin=kTRUE

root -l -q "SkimMCTree_raw_GetAverageQ.C($nevt)"
root -l -q "SkimMCTree_recentering_GetAverageEP.C($nevt)"
root -l -q "SkimMCTree_flatten_weight_GetResCor.C($nevt,$todaysDate,$flattenBinByBin,0.5)"
