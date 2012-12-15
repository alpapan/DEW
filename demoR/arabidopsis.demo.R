source('/home/pap056/workspace/dew/R/edgeR_funcs.R')
targetsFile='arabidopsis.demo.files'
datafile='tmp'; dispersion='auto'; FDR=0.05;kclusters = 10;minCPM = 2;minReps = 3
edgeR_DE_analysis_explore(targetsFile,datafile,dispersion,FDR,kclusters,minCPM,minReps)

