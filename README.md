# btseq
Bisulfite targeted sequencing analysis

This program is designed to be run on a SLURM cluster.
It is designed to process raw nextseq results and analyze methylation of sampels provided in a sample sheet.
The results will be printed to the user's current directory.
In the directory should be a tab-delimited sample sheet (named sample_sheet.txt) describing the experiment (see example).
The pipeline can take as input a BCL directory. Alternatively, if fastq.gz files are already present in the working directory, the current directory may be given as input (ex. btseq .). fastq.gz files should have the form A001_S1_R1_001.fastq.gz with the number following A00 and S signifying the order that the relevant barcode appears in the sample sheet.
The pipeline is modular and can be run in seperate steps using the relevant tools.
Note: paths in scripts may need to be adjusted.
Usage: btseq [path_to_nextseq_BCL_folder]

Steps of pipeline:
1. setup - prepare relevant input files for downstream analysis - fasta files for each sample-target pair, bc_ss.csv, BSTarget_input.txt.
2. get_reads - extract fastq.gz files from BCL directory (will fail but not disrupt pipeline if no BCL input given).
3. prepare_trim - prepare run files in trim_galore directory for trim_galore.
4. prepare_align - prepare run files for bismark alignment in run directory.
5. trim_galore - trim_galore run seperately in each fastq file.
6. alignment - dependent on completion of trim_galore. output summary are printed to Run directory. Results are printed to seperate directory for each sample (ex. './Sample_1/Output'). Output includes custom .hist and and .summary files which concisely summarize methylation results.
7. summary - creates Summary_results.txt file in main directory which summarizes methylation patterns for all samples.
(8. redo_missing - not part of built-in pipeline. May be run if parts of run failed (ex. due to faulty cluster nodes). Command is 'btseq_redo_missing'.)
(9. cleanup - not part of built-in pipeline but should be run after completion in order to remove unnecessary files. Command is 'btseq_cleanup'.)
