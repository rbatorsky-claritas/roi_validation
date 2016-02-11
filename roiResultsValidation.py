# Script for ROI genelist filtering in nextcode
# Test whether roi selects only genes of interest
# Input is one roi vcf and one exome vcf and one genelist

# example data in roi_validation_demo directory:
# python roiResultsValidation.py -roivariant roi_validation_demo/roi_test.gor -exomevariant roi_validation_demo/exome_test.gor -prefix demo -genelist roi_validation_demo/BMF_v1.genelist.txt

import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-roivariant", dest="roivariant", help="roi variant file exported from nextcode", required=True)
    parser.add_argument("-exomevariant", dest="exomevariant", help="exome variant file exported from nextcode", required=True)
    parser.add_argument("-genelist", dest="genelist", help="roi genelist file to be validated", required=True)
    parser.add_argument("-prefix", dest="prefix", help="prefix for output", required=True)

    args = parser.parse_args()
    return args.roivariant, args.exomevariant, args.genelist, args.prefix


def read_genelist(genelistFile):
    genelistTable = pd.read_csv(genelistFile, sep='\t')
    return set(genelistTable['Approved_Symbol'])


def find_variants_outside_genelist(geneSymbolList,variantFile):
    variantTable = pd.read_csv(variantFile, sep='\t')

    variantsInList = variantTable.loc[(variantTable['GENE_symbol'].isin(geneSymbolList)) &
                         (~variantTable['VEP_max_consequence'].isin(excludedVEPCategories))]

    variantsToCheck = variantTable.loc[(~variantTable['GENE_symbol'].isin(geneSymbolList)) &
                     (~variantTable['VEP_max_consequence'].isin(excludedVEPCategories))]

    variantsToCheck['inROI'] = 1

    for index, row in variantsToCheck.iterrows():
        if (~((variantsInList['CHROM'] == row['CHROM']) & (variantsInList['POS'] == [row['POS']])).any()):
            variantsToCheck['inROI'].ix[index] = 0

    variantsNotInList = variantsToCheck[variantsToCheck['inROI'] == 0]

    return variantsNotInList

def find_variants_inside_genelist(geneSymbolList,variantFile):
    variantTable = pd.read_csv(variantFile, sep='\t')

    variantsInList = variantTable.loc[(variantTable['GENE_symbol'].isin(geneSymbolList)) &
                         (~variantTable['VEP_max_consequence'].isin(excludedVEPCategories))]
    return variantsInList

## add a check that the two files have the same header, otherwise this wont work
def subtract_roi_variants_from_exome_variants():

    roiTable = pd.read_csv(roiVariantFile, sep='\t')
    exomeTable = pd.read_csv(exomeVariantFile,sep='\t')

    joinTable = pd.merge(roiTable,exomeTable,on=roiTable.columns.tolist(),how='inner')
    joinTable['common'] = 1
    temp_df = pd.merge(exomeTable, joinTable, on=exomeTable.columns.tolist(), how='left')
    exomeSubtractRoi = temp_df[temp_df['common'].isnull()].drop('common', axis=1)

    return exomeSubtractRoi


### main

excludedVEPCategories = ['downstream_gene_variant','upstream_gene_variant','5_prime_UTR_variant','3_prime_UTR_variant','intron_variant']
roiVariantFile, exomeVariantFile, genelistFile, prefix = parse_args()
geneSymbolList = read_genelist(genelistFile)

roiVariantsNotInList = find_variants_outside_genelist(geneSymbolList,roiVariantFile)
roiVariantsNotInList.to_csv(prefix + '_roi_variants_not_in_genelist.txt',sep='\t')

exomeSubtractRoi = subtract_roi_variants_from_exome_variants()

filename = 'exome_subtract_roi_variants.txt'
exomeSubtractRoi.to_csv(prefix + '_' + filename,sep='\t')

ExomeSubtractRoiInList = find_variants_inside_genelist(geneSymbolList, prefix + '_' + filename)
ExomeSubtractRoiInList.to_csv(prefix + '_exome_subtract_roi_variants_in_genelist.txt',sep='\t')
