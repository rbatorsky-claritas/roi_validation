# Test whether roi selects only genes of interest
# Input is one roi vcf and one exome vcf and one genelist


import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-roivariant", dest="roivariant", help="roi variant file", required=True)
    parser.add_argument("-exomevariant", dest="exomevariant", help="exome variant file", required=True)
    parser.add_argument("-genelist", dest="genelist", help="roi genelist file", required=True)
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
def subtract_roi_variants_from_exome():

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

exomeSubtractRoi = subtract_roi_variants_from_exome()

filename = 'exome_subtract_roi_variants.txt'
exomeSubtractRoi.to_csv(prefix + '_' + filename,sep='\t')

ExomeSubtractRoiInList = find_variants_inside_genelist(geneSymbolList, prefix + '_' + filename)
ExomeSubtractRoiInList.to_csv(prefix + '_exome_subtract_roi_variants_in_genelist.txt',sep='\t')









    #
    # chromosomes = variantsToCheck['CHROM'].tolist()
    # positions = variantsToCheck['POS'].tolist()
    # reference = variantsToCheck['Reference'].tolist()
    # variant = variantsToCheck['Call'].tolist()
    # geneSymbol = variantsToCheck['GENE_symbol'].tolist()
    #
    # variantsOutROI = {}
    #
    # outChrom = []
    # outPos = []
    # outVariant = []
    #
    # for ind in xrange(0, len(chromosomes)):
    #     if (~((variantsInROI['CHROM'] == chromosomes[ind]) & (variantsInROI['POS'] == [positions[ind]])).any()):
    #         outChrom.append(chromosomes[ind])
    #         outPos.append(positions[ind])
    #
    # variantOutROI = dict(zip(outChrom,outPos))
    # print variantOutROI
    # return variantOutROI

    #
    # print outvar
    # print "invar = " + str(len(invar))
    # print "outvar = " + str(len(outvar))
    # return outvar





#
#
# invar = roiTable.loc[(roiTable['GENE_symbol'].isin(geneSymbolList)) &
#                      (~roiTable['VEP_max_consequence'].isin(['downstream_gene_variant','upstream_gene_variant']))]
# outvar = roiTable.loc[(~roiTable['GENE_symbol'].isin(geneSymbolList)) &
#                      (~roiTable['VEP_max_consequence'].isin(['downstream_gene_variant','upstream_gene_variant']))]

#
# def read_genelist(genelistFile):
#     geneSymbolList = []
#     with open(genelistFile, 'r') as csvfile:
#         csvreader = csv.reader(csvfile, delimiter="\t")
#         for row in csvreader:
#             geneSymbolList.append(row[2])
#     return set(geneSymbolList)




# def find_exome_variants_inside_genelist():
#     df1 = pd.read_csv(roiVariantFile, sep='\t')
#     df2 = pd.read_csv(exomeVariantFile, sep='\t')
#     common_cols = roiTable.tolist()
#     df12 = pd.merge(df1,df2,on=common_cols,how='inner')
#     df2 = df2[~df2['A'].isin(df12['A'])]
#
#     print subtractionResult
#     return subtractionResult


# find_exome_variants_inside_genelist()
#
# common_cols = df1.columns.tolist()                         #generate list of column names
# df12 = pd.merge(df1, df2, on=common_cols, how='inner')     #extract common rows with merge
# df2 = df2[~df2['A'].isin(df12['A'])]
#
# http://stackoverflow.com/questions/29464234/compare-python-pandas-dataframes-for-matching-rows
