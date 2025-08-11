#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  7 09:28:07 2025

@author: mabin
"""

###Pan methylomes framework
import pandas as pd 
import numpy as np
from scipy.stats import pearsonr
import argparse
from rpy2.robjects import r as Rcode
from rpy2.robjects.packages import importr as Rrequire
from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument('--index',required=True,help="Strain sequence index location.")
parser.add_argument('--gene_index',required=True,help="Gene index location.")

parser.add_argument('--methy',required=True,help="Tombo methylation result location.")
parser.add_argument('--seq_dir',required=True,help="Sequences location.")
parser.add_argument('--motif',default='GATC',help="Target motif.")
parser.add_argument('--gene_info',required=True,help="Gene information.")


parser.add_argument('--out',required=True,help="Output location.")
parser.add_argument('--gene_suffix',default='fa',help="Sequences suffix.")


args = parser.parse_args()

def motif_index(seq,motif):
    count = seq.count(motif)
    motif_list = []
    start_location = 0
    for i in range(count):
        index = seq.index(motif,start = start_location)
        motif_list.append(index)
        start_location = index+1
        
    return motif_list


####
NGS_info = list(pd.read_csv(args.index,header = None)[1])

seq_index = list(pd.read_csv(args.gene_index,header = None)[0])

total_methy_info = []
total_seq_name_info = []
total_methy_score = []
total_unmethy_score= []
for h in range(len(seq_index)):
    gene_name = seq_index[h]

    select_gene_location = args.methy+gene_name+'.'+args.gene_suffix
    align_gene_location = args.methy+gene_name+'_mafft.'+args.gene_suffix
    # os.system('mafft --thread 20 %s > %s' %(select_gene_location,align_gene_location))
    
    align_gene_id = []
    align_gene_seq = []
    for seq_record in SeqIO.parse(select_gene_location,'fasta'):
        align_gene_id.append(seq_record.id)

    for seq_record in SeqIO.parse(align_gene_location,'fasta'):
        align_gene_seq.append(seq_record.seq.upper())
        
    ##search total motif location
    ref_motif_location = []
    for i in range(len(align_gene_id)):
        if args.motif == 'CCWGG':
            tmp_motif_location_A = motif_index(align_gene_seq[i],'CCTGG')
            tmp_motif_location_B = motif_index(align_gene_seq[i],'CCAGG')
            tmp_motif_location = np.array(tmp_motif_location_A+tmp_motif_location_B)
        else:
            tmp_motif_location = np.array(motif_index(align_gene_seq[i],args.motif))
        for j in range(len(tmp_motif_location)):
            if tmp_motif_location[j] not in ref_motif_location:
                ref_motif_location.append(tmp_motif_location[j])
                
    tmp_info_matrix = np.zeros((len(ref_motif_location),6))
    tmp_info_matrix[:,0] = np.array(ref_motif_location)
    tmp_seq_name_info = []
    for k in range(len(ref_motif_location)):
        tmp_seq_name_info.append([])
        
    tmp_methy_score = []    
    for k in range(len(ref_motif_location)):
        tmp_methy_score.append([])

    for i in range(len(align_gene_id)):

        seq_name = align_gene_id[i].split('_')[0]
        chrom = align_gene_id[i].split('_')[1]+'_'+ align_gene_id[i].split('_')[2]
        start = int(align_gene_id[i].split('_')[3].split('-')[0])
        end = int(align_gene_id[i].split('_')[3].split('-')[1].split(':')[0])
        strand = align_gene_id[i].split('_')[3].split(':')[1]
        ##calculate overlap
        seq_info_location = args.gene_info+seq_name+'_info.csv'
        seq_info = pd.read_csv(seq_info_location,header = None)

        select_gene_name = chrom+'_'+str(start)+'-'+str(end)+':'+strand
        gene_index = list(seq_info[5]).index(select_gene_name)
        overlap = seq_info[4][gene_index]
        
        if args.motif=='CCWGG':
            tmp_motif_location_A = motif_index(align_gene_seq[i],'CCTGG')
            tmp_motif_location_B = motif_index(align_gene_seq[i],'CCAGG')
            motif_location = np.array(tmp_motif_location_A+tmp_motif_location_B)
        else:
            motif_location = np.array(motif_index(align_gene_seq[i],args.motif))
        
        if len(motif_location) == 0:
            tmp_info_matrix[:,3] = tmp_info_matrix[:,3] +1
            continue
     
        #####
        if args.motif == 'GATC':
            if strand == '+':
                motif_location_filt = motif_location+start
            else:
                motif_location_filt = end-motif_location-2
        elif args.motif == 'CCWGG':
            if strand == '+':
                motif_location_filt = motif_location+start
            else:
                motif_location_filt = end-motif_location-2
        elif args.motif == 'ACTGA':
            if strand == '+':
                motif_location_filt = motif_location+start+3
            else:
                motif_location_filt = end-motif_location-5
        else:
            ######change by user
            if strand == '+':
                motif_location_filt = motif_location+start
            else:
                motif_location_filt = end-motif_location
            ####################
        
    ############

        deepsignal_location = args.methy + seq_index[h] +'combine_freq_cover.tsv'
        deepsignal_result = pd.read_table(deepsignal_location)
    
        deepsignal_select = deepsignal_result[deepsignal_result['Chr'] == chrom]
        deepsignal_select = deepsignal_select[deepsignal_select['Strand'] == strand]

        for k in range(len(motif_location_filt)):
            index = list(ref_motif_location).index(motif_location[k])
            tmp_info_matrix[index,1] = tmp_info_matrix[index,1] +1 
            tmp_info_matrix[index,4] = tmp_info_matrix[index,4] + overlap
            tmp_info_matrix[index,5] = tmp_info_matrix[index,5] + 1
            if len(deepsignal_select[deepsignal_select['Pos_start'] == motif_location_filt[k]]) != 0:
                deepsignal_match = deepsignal_select[deepsignal_select['Pos_start'] == motif_location_filt[k]].reset_index(drop=True)
                if deepsignal_match['Methylation'][0]>=1:
                    tmp_info_matrix[index,2] = tmp_info_matrix[index,2] + 1
                    
                    tmp_seq_name_info[index].append(seq_name)
                    
                    tmp_methy_score[index].append(deepsignal_match['Methylation'][0])
                    
        tmp_info_matrix[:,3] = tmp_info_matrix[:,3] +1
    
    tmp_info_matrix[:,4] = tmp_info_matrix[:,4] / tmp_info_matrix[:,5]

    
    nan_location = np.isnan(tmp_info_matrix)
    tmp_info_matrix[nan_location] = -1
    
    total_methy_info.append([gene_name,tmp_info_matrix,len(align_gene_seq[0])])
    total_seq_name_info.append(tmp_seq_name_info)
    total_methy_score.append(tmp_methy_score)

out_location = args.out + 'pan_methy_sense.tab'
with open(out_location,'w') as f:
      for i in range(len(total_methy_info)):
          gene_name = total_methy_info[i][0]
          gene_length = total_methy_info[i][2]
          for j in range(len(total_methy_info[i][1])):
              if total_methy_info[i][1][j,1] == 0:
                  continue
              f.write(gene_name+'\t'+str(total_methy_info[i][1][j,0])+'\t'+str(total_methy_info[i][1][j,1])+'\t'+str(total_methy_info[i][1][j,2])+'\t'+str(total_methy_info[i][1][j,3])+'\t'+str(total_methy_info[i][1][j,4])+'\t'+str(gene_length)+'\t')
              if len(total_seq_name_info[i][j]) == 0:
                  f.write('NA\n')
              else:            
                  for h in range(len(total_seq_name_info[i][j])):
                      if h == len(total_seq_name_info[i][j]) - 1:
                          f.write(total_seq_name_info[i][j][h]+'\t')
                      else:  
                          f.write(total_seq_name_info[i][j][h]+',')

                  for h in range(len(total_methy_score[i][j])):
                      if h == len(total_methy_score[i][j]) - 1:
                          f.write(str(total_methy_score[i][j][h])+'\n')
                      else:  
                          f.write(str(total_methy_score[i][j][h])+',')
                
                          
total_methy_info = []
total_seq_name_info = []
total_methy_score = []
total_unmethy_score= []
for h in range(len(seq_index)):
    gene_name = seq_index[h]

    select_gene_location = args.methy+gene_name+'.'+args.gene_suffix
    align_gene_location = args.methy+gene_name+'_mafft.'+args.gene_suffix
    # os.system('mafft --thread 20 %s > %s' %(select_gene_location,align_gene_location))
    
    align_gene_id = []
    align_gene_seq = []
    for seq_record in SeqIO.parse(select_gene_location,'fasta'):
        align_gene_id.append(seq_record.id)

    for seq_record in SeqIO.parse(align_gene_location,'fasta'):
        align_gene_seq.append(seq_record.seq.upper())
        
    ##search total motif location
    ref_motif_location = []
    for i in range(len(align_gene_id)):
        if args.motif == 'CCWGG':
            tmp_motif_location_A = motif_index(align_gene_seq[i],'CCTGG')
            tmp_motif_location_B = motif_index(align_gene_seq[i],'CCAGG')
            tmp_motif_location = np.array(tmp_motif_location_A+tmp_motif_location_B)
        else:
            tmp_motif_location = np.array(motif_index(align_gene_seq[i],args.motif))
        for j in range(len(tmp_motif_location)):
            if tmp_motif_location[j] not in ref_motif_location:
                ref_motif_location.append(tmp_motif_location[j])
                
    tmp_info_matrix = np.zeros((len(ref_motif_location),6))
    tmp_info_matrix[:,0] = np.array(ref_motif_location)
    tmp_seq_name_info = []
    for k in range(len(ref_motif_location)):
        tmp_seq_name_info.append([])
        
    tmp_methy_score = []    
    for k in range(len(ref_motif_location)):
        tmp_methy_score.append([])

    for i in range(len(align_gene_id)):

        seq_name = align_gene_id[i].split('_')[0]
        chrom = align_gene_id[i].split('_')[1]+'_'+ align_gene_id[i].split('_')[2]
        start = int(align_gene_id[i].split('_')[3].split('-')[0])
        end = int(align_gene_id[i].split('_')[3].split('-')[1].split(':')[0])
        strand = align_gene_id[i].split('_')[3].split(':')[1]
        ##calculate overlap
        seq_info_location = args.gene_info+seq_name+'_info.csv'
        seq_info = pd.read_csv(seq_info_location,header = None)

        select_gene_name = chrom+'_'+str(start)+'-'+str(end)+':'+strand
        gene_index = list(seq_info[5]).index(select_gene_name)
        overlap = seq_info[4][gene_index]

        if args.motif=='CCWGG':
            tmp_motif_location_A = motif_index(align_gene_seq[i],'CCTGG')
            tmp_motif_location_B = motif_index(align_gene_seq[i],'CCAGG')
            motif_location = np.array(tmp_motif_location_A+tmp_motif_location_B)
        else:
            motif_location = np.array(motif_index(align_gene_seq[i],args.motif))
                
        if len(motif_location) == 0:
            tmp_info_matrix[:,3] = tmp_info_matrix[:,3] +1
            continue
     
        #####
        if args.motif == 'GATC':
            if strand == '+':
                motif_location_filt = motif_location+start+1
            else:
                motif_location_filt = end-motif_location-3
        elif args.motif == 'CCWGG':
            if strand == '+':
                motif_location_filt = motif_location+start+2
            else:
                motif_location_filt = end-motif_location-4
        elif args.motif == 'ACTGA':
            if strand == '+':
                motif_location_filt = motif_location+start
            else:
                motif_location_filt = end-motif_location-2
        else:
            ######change by user
            if strand == '+':
                motif_location_filt = motif_location+start
            else:
                motif_location_filt = end-motif_location
            ####################
    ############

        deepsignal_location = args.methy + seq_index[h] +'combine_freq_cover.tsv'
        deepsignal_result = pd.read_table(deepsignal_location)
    
        deepsignal_select = deepsignal_result[deepsignal_result['Chr'] == chrom]
        deepsignal_select = deepsignal_select[deepsignal_select['Strand'] != strand]

        for k in range(len(motif_location_filt)):
            index = list(ref_motif_location).index(motif_location[k])
            tmp_info_matrix[index,1] = tmp_info_matrix[index,1] +1 
            tmp_info_matrix[index,4] = tmp_info_matrix[index,4] + overlap
            tmp_info_matrix[index,5] = tmp_info_matrix[index,5] + 1
            if len(deepsignal_select[deepsignal_select['Pos_start'] == motif_location_filt[k]]) != 0:
                deepsignal_match = deepsignal_select[deepsignal_select['Pos_start'] == motif_location_filt[k]].reset_index(drop=True)
                if deepsignal_match['Methylation'][0]>=1:
                    tmp_info_matrix[index,2] = tmp_info_matrix[index,2] + 1
                    
                    tmp_seq_name_info[index].append(seq_name)
                    
                    tmp_methy_score[index].append(deepsignal_match['Methylation'][0])
                    
        tmp_info_matrix[:,3] = tmp_info_matrix[:,3] +1
    
    tmp_info_matrix[:,4] = tmp_info_matrix[:,4] / tmp_info_matrix[:,5]

    
    nan_location = np.isnan(tmp_info_matrix)
    tmp_info_matrix[nan_location] = -1
    
    total_methy_info.append([gene_name,tmp_info_matrix,len(align_gene_seq[0])])
    total_seq_name_info.append(tmp_seq_name_info)
    total_methy_score.append(tmp_methy_score)

out_location = args.out + 'pan_methy_antisense.tab'
with open(out_location,'w') as f:
      for i in range(len(total_methy_info)):
          gene_name = total_methy_info[i][0]
          gene_length = total_methy_info[i][2]
          for j in range(len(total_methy_info[i][1])):
              if total_methy_info[i][1][j,1] == 0:
                  continue
              f.write(gene_name+'\t'+str(total_methy_info[i][1][j,0])+'\t'+str(total_methy_info[i][1][j,1])+'\t'+str(total_methy_info[i][1][j,2])+'\t'+str(total_methy_info[i][1][j,3])+'\t'+str(total_methy_info[i][1][j,4])+'\t'+str(gene_length)+'\t')
              if len(total_seq_name_info[i][j]) == 0:
                  f.write('NA\n')
              else:            
                  for h in range(len(total_seq_name_info[i][j])):
                      if h == len(total_seq_name_info[i][j]) - 1:
                          f.write(total_seq_name_info[i][j][h]+'\t')
                      else:  
                          f.write(total_seq_name_info[i][j][h]+',')

                  for h in range(len(total_methy_score[i][j])):
                      if h == len(total_methy_score[i][j]) - 1:
                          f.write(str(total_methy_score[i][j][h])+'\n')
                      else:  
                          f.write(str(total_methy_score[i][j][h])+',')
                          
















































