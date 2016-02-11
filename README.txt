

This is a script to validate that ROI genelist filtering is done correctly for studies in NextCode.

Dependencies:
Tested with ython 2.7.10
Python packages pandas and argparse are needed
 
Run on example data in roi_validation_demo directory:
python roiResultsValidation.py -roivariant roi_validation_demo/roi_test.gor -exomevariant roi_validation_demo/exome_test.gor -prefix BMF_demo -genelist roi_validation_demo/BMF_v1.genelist.txt

The input are:
- Two gor tables exported from nextcode containing the variants in the roi study and the exome study, e.g. roi_test.gor and exome_test.gor. 
- Genelist used to product the roi data from the exome data, e.g. BMF_v1.genelist.txt
- Prefix for the output files

The two gore files are produced by: 
1) opening a study in nextcode
2) generate advanced report
3) open the all variants->all by variant tab
4) Click button "open query in sequence miner" 
4) Selecting columns (You can also select all columns in all categories, but it might be clearer to interpret if only select columns appear)
Basic: CHROM, POS, Reference, Call
VEP: _max_consequence
GENE: symbol
5) Click disk icon to save gor file to the user_data folder on the nextcode clinical cloud directory (not a nested folder)
6) Right click or click on the created file in user_data and select copy
7) Navigate to local directory and paste the copied file

The three outputs are:
- BMF_demo_exome_subtract_roi_variants.txt: Table of variants found in the exome but not in the roi
- BMF_demo_exome_subtract_roi_variants_in_genelist.txt: Subset of above that are in the genelist
- BMF_demo_roi_variants_not_in_genelist.txt: Subset roi variants that are not in the genelist

No variants should be found in the later two categories. If they are, it should be investigated.
