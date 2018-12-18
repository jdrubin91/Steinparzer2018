__author__ = 'Jonathan D. Rubin'

#==============================================================================
#Imports
#==============================================================================
from Steinparzer2018 import functions

def differential_transcription(comparisons=None, raw_data=None, 
                                bedtools_intersect=None, 
                                bedtools_multicov=None, bedtools_merge=None, 
                                bedtools_sort=None, output_folder=None):

    functions.intersect(comparisons=comparisons, raw_data=raw_data, 
                        bedtools_intersect=bedtools_intersect, 
                        output_folder=output_folder)

    functions.merge(comparisons=comparisons, bedtools_sort=bedtools_sort, 
                    bedtools_merge=bedtools_merge, output_folder=output_folder)

    functions.count(comparisons=comparisons, output_folder=output_folder, 
                    raw_data=raw_data, bedtools_multicov=bedtools_multicov) 

    functions.DESeq2(comparisons=comparisons, raw_data=raw_data, 
                    output_folder=output_folder)

    functions.DESeq_to_bed(comparisons=comparisons, 
                            output_folder=output_folder)

def MDS_analysis():
    return None

def differential_MDS():
    return None

if __name__ == "__main__":
    differential_transcription()