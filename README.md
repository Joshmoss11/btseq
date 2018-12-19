# btseq
Bisulfite targeted sequencing analysis

This program is designed to be run on a SLURM cluster.
It is designed to process raw nextseq results and analyze methylation of sampels provided in a sample sheet.
The results will be printed to the user's current directory.
In the directory should be a tab-delimited sample sheet (named sample_sheet.txt) describing the experiment (see example).
Note: paths in scripts may need to be adjusted.
Usage: btseq [path_to_nextseq_BCL_folder]
