# Pyrogenic organic matter dark remineralization raw data repository

Collaborators and contributors: Zachary A. Quinlan, Chris B. Wall, Loreto Paulino Jr., Riley Johnson, Justin Han, Craig E. Nelson, Claire Moreland-Ochoa

## Disclaimer:
This is where we will store all of the raw data that can be used for publications on the dark remineralization experiment. If you would like to use this, you can add this reposoitory as a submodule to your own repository. 

## How to use:
```
git submodule add <url> <path>
```

If the submodule is out of date within your repository and needs to be updated you can cd into the submodule location and then run:
```
git checkout <branch>
git pull
git cd <higher level repo>
git add -A 
git commit -m <commit message>
git push
```

### Funding sources:
This research was sponsored by the US National Science Foundation Division of Ocean Sciences program grant #2516814

## Overview:
### AutoBOD:
- Oxygen measurements

### FCM:
- Flow cytometry data

### Metabolomics:
Results of LC MS/MS
- backgroundRemoval.R: Code to take the raw data and remove the backgound features defined from the methods blanks
- cleanedMSData.csv: The cleaned data exported from the backgroundRemoval.R
- libraryHits.tsv: library match data from GNPS2
- molecularNetworking.tsv: molecular subnetwork information for each features
- msMetadata.csv: metadata sheet to make sense of the filnames in the raw data
- pyroDOM_GNPS_quant.csv: quantification of ion feature intensity for the entire run
- qcPlots/: subdirectory that holds two QC plots
