# CHYRON

The code repository associated with the paper: *insert link here*

~~~~~~~~~~~~~~~~~~ * NGS Sequencing Pipeline * ~~~~~~~~~~~~~~~~~~~~

We used a simple wrapper to feed the necesary information into pipeline.sh:

	$1 = the path to the forward read fastq file
	$2 = the path to the reverse read fastq file
	$3 = the path to the barcode file with the following format:
		"Sample name"	"forward barcode"	"reverse barcode"	"20bp around cutsite"	"Experiment"
	$4 = the path to the reference sequence, with the format:
		"Experiment name"	"reference sequence"
	$5 = The name of the experiment, for file naming purposes


You'll need to download:
1) PEAR, to pair Paired End Illumina reads: https://sco.h-its.org/exelixis/web/software/pear/

And then compile:

2) Mapp, an optimal sequence alignment used in Perli et al, 2016. The C++ code is available at: http://www.rle.mit.edu/sbg/resources/stgRNA/

I) Perli, S. D., Cui, C. H. & Lu, T. K. Continuous genetic recording with self-targeting CRISPR-Cas in human cells. Science 353, aag0511 (2016).

II) Zhang, J., Kobert, K., Flouri, T. & Stamatakis, A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30, 614–620 (2014).

We used Enthought Python 7.3.2 and Python 2.7.15

This pipeline was written by Elmira Forouzmand and Mason Schechter

~~~~~~~~~~~~~~~~~ * Average Prefix Reconstruction Method * ~~~~~~~~~~~~~~~~~~

You will need ins_parsing_no_umi.py, lib.py, and VectorSpaceWellMatrix.py

1) Gather the .newinfo.withlengths.sorted.txt output files from NGS pipeline
2) In a separate text file (generally, we call it "exp"_well_data.txt), put each output file from step 1 into its own line
3) Run ins_parsing_no_umi.py with the well_data.txt file as the the only command line argument. The output will be a pickle with a .pkl ending. To run:
```
    $ python ins_parsing_no_umi.py myExp_well_data.txt ```
    
4) To run the well pairing script, VectorSpaceWellMatrix.py, change the SYSTEM parameter "PICKLE_DATA_FILE" to the name of the output .pkl from step 3 and adjust your "LEN_CUTOFF" (insertion length)
and "COUNT_CUTOFF" (insertion count) to produce your desired lineage reconstruction. To run:
```
    $ python VectorSpaceWellMatrix.py ```
    
This reconstruction algorithm was written by Mason Schechter and Bijan Agahi

