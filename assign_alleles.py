#!/usr/bin/python

#This script will take in a csv file of extracted genotypes from GATK v4.2.0.0 (use extract.gt() function from vcf package in R) and convert variant calls (1/1) into parental allele values (R or S for resistant and sensitive) to be used by r/qtl2 in a genome scan
#The vcf file should have been filtered to eliminate poor quality variant sites
#Output will be a data frame in csv format with allele information for all markers represented in the input vcf file
#usage: /path/to/assign_alleles.py /path/to/in.vcf R_parent_name S_parent_name /path/to/out.csv

from sys import argv, version_info
from os.path import realpath, splitext
import os

in_csv, R_name, S_name, out_csv = realpath(argv[1]), argv[2], argv[3], realpath(argv[4])

line_count=0
parts=[]
R_parent={}
S_parent={}
col_count=0
R_col=0
S_col=0
subparts=[]
in_fh=open(in_csv, 'r')
out_fh=open(out_csv, 'w+')

#loop through the genotypes file
for line in in_fh:
  line=line.strip('\n')
  parts=line.split(',')
  #for the first line, find which column the R and S parent data are housed in
  if line_count==0:
    for sample in parts:
      sample=sample.strip('"')
      sample=sample.strip("'")
      if sample==R_name:
        R_col=col_count
      elif sample==S_name:
        S_col=col_count
      col_count+=1
    line_count+=1
  #for the rest of the lines, store the R and S parent genotypes in an array within a dictionary with line number (marker) as the key and genotype as two numbers in an array
  else:
    col_count=0
    for sample in parts:
      sample=sample.strip('"')
      sample=sample.strip("'")
      if col_count==S_col:
        if '/' in sample:
          S_parent[line_count]=sample.split('/')
        elif '|' in sample:
          S_parent[line_count]=sample.split('|')
      elif col_count==R_col:
        if "/" in sample:
          R_parent[line_count]=sample.split('/')
        elif "|" in sample:
          R_parent[line_count]=sample.split('|')
      col_count+=1
    line_count+=1
in_fh.close
line_count=0

#loop through the file again, this time figuring out which parent donated the allele and converting 1's and 0's into R's or S's
for line in open(in_csv, 'r'):
  line=line.strip('\n')
  if line_count==0:
    out_fh.write(str(line)+'\n')
    line_count+=1
  else:
    if R_parent[line_count]!=S_parent[line_count]:
      col_count=0
      R_alleles=[]
      S_alleles=[]
      parts=line.split(',')
      out_fh.write(str(parts[0]))
        
      for i in R_parent[line_count]:
        if i not in S_parent[line_count]:
          R_alleles.append(i)
      for i in S_parent[line_count]:
        if i not in R_parent[line_count]:
          S_alleles.append(i)
      if R_alleles == []:
        R_alleles.append(R_parent[line_count][1])
      elif S_alleles == []:
        S_alleles.append(S_parent[line_count][1])
                    
      for sample in parts:
        if col_count!=0:
          sample=sample.strip('"')
          sample=sample.strip("'")
          out_fh.write(',')
          if "/" in sample:
            subparts=sample.split('/')
          elif "|" in sample:
            subparts=sample.split('|')
          
          for allele in subparts:
            if allele not in R_alleles and allele not in S_alleles:
              out_fh.write("NA")
            elif allele in R_alleles:
              out_fh.write("R")
            elif allele in S_alleles:
              out_fh.write("S")
        col_count+=1
    out_fh.write('\n')
            
              
            
    line_count+=1
in_fh.close
