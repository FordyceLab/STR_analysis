# STR_MITOMI_analysis

Accompanying code for "Short tandem repeats recruit transcription factors to tune eukaryotic gene expression." Accompanying data for the study can be found at https://osf.io/gbxhz.

Citation: COMING SOON

### AffinityDistillation
AffinityDistillation code can be found at `TBD`. MAX ChIP-seq Pho4 ChIP-nexus training data can be found at `TBD`.

### Affinity predictions
Here you can find a Jupyter notebook with code to predict binding with PSAMs, a sliding Z-score model, and a partition function-based model. Sample data to predict MAX and Pho4 binding to the sequences used in this work are provided in the repo.

Additional data to predict binding for your TF of interest can be found:
- PSAMs: https://jaspar.genereg.net/
- uPBM 8-mer Z-scores: https://osf.io/gbxhz under `uPBM_data` folder
- uPBM 8-mer intensity scores: http://cisbp.ccbr.utoronto.ca/entireDownload.php (under `PBM Intensities`. Corresponding information under `TF Information`) or http://the_brain.bwh.harvard.edu/uniprobe/.

### Enhancer enrichment
Here you can find two Jupyter notebooks to calculate enrichment of repeats in STRs and enrichment for broadly active/most active enhancers. Due to GitHub repo storage limits, the necessary data files are located in our OSF repository: https://osf.io/gbxhz under `enhancer_annotations` and `TRF_output`. STRs were located in the human genome using https://github.com/Benson-Genomics-Lab/TRF. Because the code to calculate enrichment for the genomic data and shuffled take several hours to run, we also provide the results of our calculations that were plotted in Figure 7. The code to reproduce the enrichment calculations is provided at the bottom of the notebook if you wish to rerun it.

### Gillespie simulation
Here you can find code to run Gillespie simulations for TF search in the yeast nucleus. We also provide bash scripts for running these simulations on a SLURM-based HPC. (Analysis and code by J. Schaepe)

### Kinetic model fitting
Here you can find code for fitting microscopic rate parameters in a 4-state kinetic model of TF-DNA binding. To repeat the analysis, run `model_kinetic_repeats/model_kinetic_repeats_211014_reallynewnorm_Pho4.m` and `model_kinetic_repeats/model_kinetic_repeats_211014_reallynewnorm_Max.m` (Analysis and code by E. Marklund)

### Sample MITOMI notebooks
Here you can find sample Jupyter notebooks used to analyze MITOMI binding data, k-MITOMI kinetic data, or STAMMP binding data. We also include some sample data for each type of analysis.

### Statistical mechanics simulations
Here you can find code to repeat simulations interrogating the role of entropy and enthalpy in binding STRs.

### uPBM analysis
Here you can find code comparing repeat differences between paralogous TFs (`paralog_analysis.ipynb`) and comparing TF binding motifs to repeat preferences (`motif_vs_repeat.ipynb`).