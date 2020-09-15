[![CC BY 4.0][cc-by-shield]][cc-by]

# mphil-thesis

Source for the PDF is in the LaTeX directory. Full text of the thesis is in `Edridge DSouza MPhil Thesis.pdf`.

Within the `Computation` directory, each subdirectory represents a different portion of the project. All files are provided except for the files in the `data` directories, as well as any `.RDS` data files. However, following the `00_setup.sh` script in any of the project directories will create the data directory necessary to run the scripts.

- `2019_11_04_scRNAseq analysis` contains the embryonic dataset from Karaiskos et al. (referred to as `zinzen` in the code). 

- `2020_04_01_Avalos_scRNA` contains the larval brain dataset from Avalos et al. (referred to as `avalos`). 

- `2020_04_29_Allen_VentralNerveCord` contains the adult VNC dataset from Allen et al. (referred to as `allen`). 

- `2020_04_29_ClusterMap` contains the ClusterMap attempt on the local Ubuntu machine (the remainder can be found in `Rustbucket Analyses`). 

- `2020_07_03_PutItTogether` contains an exhaustive set of Gene Ontology plots for all the SoxB-positive clusters in each of the datasets. 

- `Cytoscape` contains figures as well as the network file for the Cytoscape analysis.

- `Rustbucket Analyses` contains the "regular" and "by the book" analyses for the adult VNC dataset, which was performed on the Rustbucket server. It also contains the other half of the ClusterMap attempted analysis.