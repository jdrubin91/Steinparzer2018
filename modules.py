__author__ = 'Jonathan D. Rubin'

#==============================================================================
#Imports
#==============================================================================
import os
<<<<<<< HEAD
import math
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

=======
>>>>>>> ff222d89782f2fec243caa8b18c03b7823a7a1d7
import functions
import matplotlib.pyplot as plt

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

def differential_MDS(comparisons=None, output_folder=None):
    print("Performing differential MDS analysis...")
    for condition1, condition2 in comparisons:
        print('\t', condition1, 'vs.', condition2)
        MDS1 = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.' + condition1
                            + '.MDS_MDS.csv')

        MDS2 = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.' + condition2
                            + '.MDS_MDS.csv')

        output = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.MDS_diff.txt')

        X, Y, pvals, genelist = functions.MDS(MDS1=MDS1, MDS2=MDS2)
        functions.write_MDS_output(genelist=genelist, X=X, Y=Y, p_vals=pvals, 
                            output=output)

def MA_plot_figure(comparisons=None, output_folder=None, save=False, 
                    label=False):
    F = plt.figure(figsize=(15,15))
<<<<<<< HEAD
    subplot_n = math.ceil(len(comparisons)/2)
=======
    subplot_n = int(len(comparisons)/2)
>>>>>>> ff222d89782f2fec243caa8b18c03b7823a7a1d7
    for i, (condition1, condition2) in enumerate(comparisons):
        subplot = (subplot_n*100) + (subplot_n*10) + i + 1
        ax = F.add_subplot(subplot)
        functions.plot_MA(condition1=condition1, condition2=condition2, 
                            output_folder=output_folder, ax=ax, label=label)
    plt.tight_layout()
    plt.show()
    if save:
        save_path = os.path.join(output_folder, condition1 + '_' + condition2
                                                + '.MA-plot.svg')
        plt.savefig(save_path)

<<<<<<< HEAD
def heatmap_plot_figure(comparisons=None, output_folder=None, save=False, 
                        motif=None):

    F = plt.figure(figsize=(15,15))
    print(motif)
    subplot_n = math.ceil(len(comparisons)/2)
    outer = gridspec.GridSpec(subplot_n, subplot_n, wspace=0.25, hspace=0.25)
    for i, (condition1, condition2) in enumerate(comparisons):
        MDS1 = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.' + condition1
                            + '.MDS_MDS.csv')

        MDS2 = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.' + condition2
                            + '.MDS_MDS.csv')

        inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[i], wspace=0.2, hspace=0.2)

        ax = plt.Subplot(F, inner[0])
        F.add_subplot(ax)
        functions.plot_heatmap(input_file=MDS1, motif=motif, ax=ax)

        ax = plt.Subplot(F, inner[1])
        F.add_subplot(ax)
        functions.plot_heatmap(input_file=MDS2, motif=motif, ax=ax)
    if save:
        save_path = os.path.join(output_folder, condition1 + '_' + condition2
                                                + '.heatmaps.svg')
        plt.savefig(save_path)
    plt.show()
=======

>>>>>>> ff222d89782f2fec243caa8b18c03b7823a7a1d7

if __name__ == "__main__":
    differential_transcription()