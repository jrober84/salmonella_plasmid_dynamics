# Salmonella plasmid dynamics framework
Analysis and visualization tools for plasmid manuscript
The scripts used to process the MOB-suite and abricate results into visualizations are provided here in addition to HTML versions of the plots and the processed result files.

#Prepare input files for use in visualization tools
Due to size limitations in github, the input files can be downloaded from zenodo here: []

python perform_analysis.py --abricate ./data/abricate-ncbi-amr-genes.txt.gz --contigs /data/salmonella-contig-reports.txt.gz --mobtyper ./data/2020-11-Salmonella-mobtyper_results.txt.gz --metadata ./data/metadata.txt.gz --outdir ./results

#Create visualizations
python create_plots.py


