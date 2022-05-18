#!/usr/bin/env python

#import necessary packages
import pandas as pd
import sys;

#read the files
colnames = ['CHROM','POS','ID','AF','REF','ALT','GENE','BIOTYPE','IMPACT','HGVSP', 'GT', 'AD', 'DP','GQ','PL']
variants = pd.read_table(sys.argv[1],sep='\t', dtype = 'str',names = colnames, header = None)
variants["CHROM"].replace("chr", "",regex = True,  inplace=True)

gnomad = pd.read_table(sys.argv[2],sep='\t', dtype = 'str')


# combine the files
merged = pd.merge(variants , gnomad,on = ["CHROM","POS","REF","ALT"], indicator = True)
merged.to_csv(sys.argv[1] + '_left_join.txt', sep="\t",index = False)
