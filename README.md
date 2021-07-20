# Gyrfalcon Mark-Recapture-Recovery modelling

Analyses of capture-recapture-recovery/resighting data from Iceland by F. Barraquand and O.K. Nielsen. We estimate survival rates for adults and juvenile, as well as effects of prey abundance and weather on juvenile survival. Final estimation is performed in Stan within the `MRR` folder, where we fit the Mark-Recapture-Recovery models. These combine recaptures from live and dead animals. These have been modified from Kéry & Schaub's BPA Chapter 9 (bear in mind, we do use different models though!). The Stan code has itself been adapted from the [Stan translation of BPA codes by Hiroki Itô](https://github.com/stan-dev/example-models/tree/master/BPA). RMarkdown files at the base of the repo produce the capture histories from the raw data. 
