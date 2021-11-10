#!/usr/bin/python

import os, sys, re, glob, json
import pandas as pd
from pandas.tseries import frequencies

flist = glob.glob(
    "/storm_KDNA/Analysis/KDNA/Tools/pharmcat/test/pharm_vcf/*_v1.0.0.pharmcat.report.json"
)
countall = []
count_genotype = 0
frequency = 0
table = []
cnt_table = []
count_sample = {}
DIC = {}
for infile in flist:
    fpin = json.load(open(infile))
    SAMPLE = infile.split("/")[-1].split("_")[0]

    for line in fpin["genotypes"]:
        Gene = str(line["gene"])
        Genotype = str("|".join(line["calls"]))
        Function = str("|".join(line["functions"]))
        Phenotype = str("|".join(line["phenotype"]))
        missingVariants = str(line["missingVariants"])
        if SAMPLE not in DIC:
            DIC[SAMPLE] = {}
        if Gene not in DIC:
            DIC[SAMPLE][Gene] = {}
        DIC[SAMPLE][Gene]["Genotype"] = Genotype
        DIC[SAMPLE][Gene]["Function"] = Function
        DIC[SAMPLE][Gene]["Phenotype"] = Phenotype
        DIC[SAMPLE][Gene]["missingVariants"] = missingVariants
        # print(DIC)

        # cnt_genotype
        if Genotype != "not called":
            if Gene not in count_sample:
                count_sample[Gene] = {}
            if Genotype not in count_sample[Gene]:
                count_sample[Gene][Genotype] = 1
            else:
                count_sample[Gene][Genotype] += 1


#################
####cnt sample###
#################

allitem = count_sample.items()
gene_cnt = len(allitem)

for k in range(0, gene_cnt):
    aa = allitem[k][0]
    key_list = allitem[k][1].keys()
    key_cnt = len(allitem[k][1].keys())
    value_list = allitem[k][1].values()
    for i in range(0, key_cnt):
        bb = key_list[i]
        cc = value_list[i]
        cnt_table.append([aa, bb, cc])

cnt_total = pd.DataFrame(cnt_table, columns=["Gene", "Genotype", "Frequency"])
cnt_total.to_excel("Genotype_Frequency_final.xlsx")


##################
##total summary###
##################

if not Genotype == "not called":
        table.append(
            [SAMPLE, Gene, Genotype, Function, Phenotype, missingVariants]
        )
        df_total = pd.DataFrame(
            table,
            columns=[
                "SampleID",
                "Gene",
                "Genotype",
                "GeneFrequency",
                "Phenotype",
                "MissingVariants",
            ],
        )
df_total.to_excel("Total_summary_nocall.xlsx")
