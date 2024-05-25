from Bio import AlignIO
from collections import defaultdict
from matplotlib import pyplot as plt
import os
import csv
import argparse

def main():
    arguments=argument_parser()
    
    with open(arguments.alignment_path) as fh:
        MSA=AlignIO.read(fh,'clustal')
    
    #create folder named conserved_epitopes_argparse if it doesn't already exist
    
    file_name=os.path.basename(arguments.alignment_path).split('.')[0]
    file_name=os.path.join('conserved_epitopes_argparse',file_name+'_conserved_epitopes.csv')
    
    epitope_counts=get_epitope_counts(MSA)
    epitope_conservations=get_epitope_conservations(epitope_counts,MSA,file_name)
    conserved_epitopes_to_csv(epitope_conservations, file_name)
    
#parse input file path with argparse
def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment_path", help='Path to the alignment file')
    args = parser.parse_args()
    
    return args

#break the multiple sequence alignment into 8, 9, 10, and 11mer epitopes
def get_epitopes(sequence):
    epitopes=[]

    for i in range(8,12):

        for j in range(len(sequence)-i):
            epitope=sequence[j:j+i]
            if 'X' not in epitope:
                epitopes.append(epitope)
    
    return epitopes

#get counts of each epitope 
def get_epitope_counts(MSA):
    epitope_counts=defaultdict(set)

    for sequence in MSA:
        sequence_nopadding=sequence.seq.replace('-','')
        epitopes_list=get_epitopes(sequence_nopadding)
        for epitope in epitopes_list:
            isolate_set=epitope_counts[epitope]
            isolate_set.add(sequence.id)
    return epitope_counts

#calculate conservations of each epitope across the alignment and return histogram
def get_epitope_conservations(epitope_counts,MSA,file_name):

    epitope_conservations={}

    for epitope, isolate_set in epitope_counts.items():
        epitope_cons=len(isolate_set)/len(MSA)
        epitope_conservations[epitope]=epitope_cons

    plt.hist(epitope_conservations.values(),label=file_name)
    plt.legend()
    plt.xlabel('Percent conservation')
    plt.ylabel('Number of epitopes')
    plt.savefig(file_name+'_conservation_plot.jpg',dpi=300)
    
    return epitope_conservations

#for each epitope with conservation >=0.9, write out to csv
def conserved_epitopes_to_csv(epitope_conservations,file_name):
    
    with open(file_name,'wt') as fh:
        writer=csv.writer(fh, delimiter=',')
        writer.writerow(['epitope','conservation'])
        
        for epitope, conservation in epitope_conservations.items():
            if conservation >=0.9:
                writer.writerow([epitope,conservation])


if __name__ == '__main__':
    main()