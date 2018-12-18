__author__ = 'Jonathan D. Rubin'

#==============================================================================
#Imports
#==============================================================================
import functions

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

def MDS_analysis(comparisons=None, output_folder=None, pval_cutoff=False, 
                    percentage_rank_cutoff=False, bedtools_intersect=None, 
                    MDS_src=None, fasta_file=None, PSSM_DB=None):

    functions.pull_regions_from_comparisons(comparisons=comparisons, 
                                            padj_cutoff=pval_cutoff, 
                                            percentage_rank_cutoff=percentage_rank_cutoff, 
                                            output_folder=output_folder, 
                                            padj_column=-1,
                                            bedtools_intersect=bedtools_intersect)

    functions.run_MDS_from_comparisons(comparisons=comparisons, 
                                        output_folder=output_folder, 
                                        MDS_src=MDS_src, fasta_file=fasta_file,
                                        PSSM_DB=PSSM_DB)

def differential_MDS(comparisons=None):
    for condition1, condition2 in comparisons:
        MDS1=''
        MDS2=''
        output=''
        X, Y, pvals, genelist = functions.MDS(MDS1=MDS1, MDS2=MDS2)
        functions.write_MDS_output(genelist=genelist, X=X, Y=Y, p_vals=pvals, 
                            output=output)
    return None

if __name__ == "__main__":
    differential_transcription()