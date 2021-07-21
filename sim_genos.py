#!/usr/bin/python

#this script will take in some arguments, 2 chromosome lengths rounded to the nearest Mbp and location of desired QTL on those chromosomes to the nearest Mbp, and output a VCF file with genotypes for a herbicide resistant and susceptible parents as well as for 200 segregating F2s
#the first 50 will be homozygous for the high effect resistant QTL (13 homozygous for low-effect, 25 heterozygous, 12 homozygous sensitive), next 100 heterozygous for high effect (25 homozygous low effect, 50 heterozygous, 25 homozygous sensitive), last 50 homozygous sensitive for high-effect (same distribution for low effect)
#usage /path/to/simulate_qtl_genotypes.py chr_length1 chr_length2 QTL1_location QTL2_location path/to/out.vcf

from sys import argv, version_info
from os.path import realpath, splitext
import os
import random

chr6_len, chr1_len, qtl1_loc, qtl2_loc, out_vcf = int(argv[1]), int(argv[2]), int(argv[3]), int(argv[4]), realpath(argv[5])

#create file handle for output file
out_fh=open(out_vcf, 'w+')

#create a dictionary to store each individual as a key and an array of genotypes 
chr6={}
chr1={}

#load the dictionaries with keys ,R parent, S parent, and 1 through 200 for the F2s, and initialize arrays for each entry
chr6['R_parent']=[]
chr6['S_parent']=[]
chr1['R_parent']=[]
chr1['S_parent']=[]
for a in range(200):
	chr6[a+1]=[]
	chr1[a+1]=[]

#load the dictionaries with genotypes for the parents
for a in range(chr6_len+1):
	chr6['R_parent'].append('1/1')
	chr6['S_parent'].append('0/0')
for a in range(chr1_len+1):
	chr1['R_parent'].append('1/1')
	chr1['S_parent'].append('0/0')
#load the arrays of each F2 with filler so the arrays are ready to take the genotypes
for a in range(200):
	for b in range(chr6_len+1):
		chr6[a+1].append('a')
	for b in range(chr1_len+1):
		chr1[a+1].append('a')

rand=0
#start at the QTL positions and step out 1 MB at a time, using a random number generator to determine if recombination changes your genotype at the next marker
for a in range(chr6_len+1-qtl1_loc):
	if a==0:
		for b in range(200):
			if b+1 <51:
				chr6[b+1][a+qtl1_loc]=str('1/1')
			elif b+1 <151:
				chr6[b+1][a+qtl1_loc]=str('0/1')
			else:
				chr6[b+1][a+qtl1_loc]=str('0/0')
	else:
		for b in range(200):
			rand=random.randint(1, 1000)
			if chr6[b+1][a+qtl1_loc-1]=='1/1':
				if rand==1:
					chr6[b+1][a+qtl1_loc]=str('0/0')
				elif rand<22:
					chr6[b+1][a+qtl1_loc]=str('0/1')
				else:
					chr6[b+1][a+qtl1_loc]=str('1/1')
			elif chr6[b+1][a+qtl1_loc-1]=='0/1':
				if rand<11:
					chr6[b+1][a+qtl1_loc]=str('0/0')
				elif rand>989:
					chr6[b+1][a+qtl1_loc]=str('1/1')
				else:
					chr6[b+1][a+qtl1_loc]=str('0/1')
			else:
				if rand==1:
					chr6[b+1][a+qtl1_loc]=str('1/1')
				elif rand<22:
					chr6[b+1][a+qtl1_loc]=str('0/1')
				else:
					chr6[b+1][a+qtl1_loc]=str('0/0')

#start at the qtl position and work backwards to the beginning of chr6, additing genotypes along the way
for a in range(qtl1_loc):
	for b in range(200):
		rand=random.randint(1,1000)
		if chr6[b+1][qtl1_loc-a]=='1/1':
			if rand==1:
				chr6[b+1][qtl1_loc-a-1]=str('0/0')
			elif rand<22:
				chr6[b+1][qtl1_loc-a-1]=str('0/1')
			else:
				chr6[b+1][qtl1_loc-a-1]=str('1/1')
                if chr6[b+1][qtl1_loc-a]=='0/1':
                        if rand<11:
                                chr6[b+1][qtl1_loc-a-1]=str('0/0')
                        elif rand>989:
                                chr6[b+1][qtl1_loc-a-1]=str('1/1')
                        else:
                                chr6[b+1][qtl1_loc-a-1]=str('0/1')
                if chr6[b+1][qtl1_loc-a]=='0/0':
                        if rand==1:
                                chr6[b+1][qtl1_loc-a-1]=str('1/1')
                        elif rand<22:
                                chr6[b+1][qtl1_loc-a-1]=str('0/1')
                        else:
                                chr6[b+1][qtl1_loc-a-1]=str('0/0')

#do the same two steps for chr1
for a in range(chr1_len+1-qtl2_loc):
        if a==0:
                for b in range(200):
                        if b+1<13:
                                chr1[b+1][a+qtl2_loc]=str('1/1')
                        elif 12<b+1<39:
                                chr1[b+1][a+qtl2_loc]=str('0/1')
                        elif 38<b+1<51:
                                chr1[b+1][a+qtl2_loc]=str('0/0')
                        if 50<b+1<76:
                                chr1[b+1][a+qtl2_loc]=str('1/1')
                        elif 75<b+1<126:
                                chr1[b+1][a+qtl2_loc]=str('0/1')
                        elif 125<b+1<151:
                                chr1[b+1][a+qtl2_loc]=str('0/0')
                        if 150<b+1<163:
                                chr1[b+1][a+qtl2_loc]=str('1/1')
                        elif 162<b+1<189:
                                chr1[b+1][a+qtl2_loc]=str('0/1')
                        elif 188<b+1<201:
                                chr1[b+1][a+qtl2_loc]=str('0/0')
        else:
                for b in range(200):
                        rand=random.randint(1, 1000)
                        if chr1[b+1][a+qtl2_loc-1]=='1/1':
                                if rand==1:
                                        chr1[b+1][a+qtl2_loc]=str('0/0')
                                elif rand<22:
                                        chr1[b+1][a+qtl2_loc]=str('0/1')
                                else:
                                        chr1[b+1][a+qtl2_loc]=str('1/1')
                        elif chr1[b+1][a+qtl2_loc-1]=='0/1':
                                if rand<11:
                                        chr1[b+1][a+qtl2_loc]=str('0/0')
                                elif rand>989:
                                        chr1[b+1][a+qtl2_loc]=str('1/1')
                                else:
                                        chr1[b+1][a+qtl2_loc]=str('0/1')
                        else:
                                if rand==1:
                                        chr1[b+1][a+qtl2_loc]=str('1/1')
                                elif rand<22:
                                        chr1[b+1][a+qtl2_loc]=str('0/1')
                                else:
                                        chr1[b+1][a+qtl2_loc]=str('0/0')

#start at the qtl position and work backwards to the beginning of chr6, additing genotypes along the way
for a in range(qtl2_loc):
        for b in range(200):
                rand=random.randint(1,1000)
                if chr1[b+1][qtl2_loc-a]=='1/1':
                        if rand==1:
                                chr1[b+1][qtl2_loc-a-1]=str('0/0')
                        elif rand<22:
                                chr1[b+1][qtl2_loc-a-1]=str('0/1')
                        else:
                                chr1[b+1][qtl2_loc-a-1]=str('1/1')
                if chr1[b+1][qtl2_loc-a]=='0/1':
                        if rand<11:
                                chr1[b+1][qtl2_loc-a-1]=str('0/0')
                        elif rand>989:
                                chr1[b+1][qtl2_loc-a-1]=str('1/1')
                        else:
                                chr1[b+1][qtl2_loc-a-1]=str('0/1')
                if chr1[b+1][qtl2_loc-a]=='0/0':
                        if rand==1:
                                chr1[b+1][qtl2_loc-a-1]=str('1/1')
                        elif rand<22:
                                chr1[b+1][qtl2_loc-a-1]=str('0/1')
                        else:
                                chr1[b+1][qtl2_loc-a-1]=str('0/0')

#write some comment lines at the top of the vcf
out_fh.write('##fileformat=VCFv4.1'+'\n'+'##source="GATK haplotype Caller"'+'\n'+'##some more comment lines about how the variant calling was done'+'\n')

#write the column names
out_fh.write('#CHROM'+'\t'+'POS'+'\t'+'ID'+'\t'+'REF'+'\t'+'ALT'+'\t'+'QUAL'+'\t'+'FILTER'+'\t'+'INFO'+'\t'+'FORMAT')

#loop through one of the dictionaries and write the sorted keys as column names
for key in sorted(chr6):
	out_fh.write('\t'+str(key))
out_fh.write('\n')

#loop through the chr1 dictionary and write information about each marker
for marker in range(len(chr1[1])):
	out_fh.write("chr1"+'\t'+str(int(marker)*1000000)+'\t'+"."+'\t'+"A"+'\t'+"T"+'\t'+str(random.randint(4000, 400000)/100)+'\t'+"."+'\t'+"."+'\t'+"GT:AD:DP:GQ:PL")
	for a in sorted(chr1):
		out_fh.write('\t'+str(chr1[a][marker])+":0,0:0:0,0:0")
	out_fh.write('\n')


#loop through the chr6 dictionary and write information about each marker
for marker in range(len(chr6[1])):
        out_fh.write("chr6"+'\t'+str(int(marker)*1000000)+'\t'+"."+'\t'+"A"+'\t'+"T"+'\t'+str(random.randint(4000, 400000)/100)+'\t'+"."+'\t'+"."+'\t'+"GT:AD:DP:GQ:PL")      
        for a in sorted(chr6):
                out_fh.write('\t'+str(chr6[a][marker])+":0,0:0:0,0:0")
	out_fh.write('\n')

