__author__ = 'Jonathan D. Rubin'

#==============================================================================
#Imports
#==============================================================================
import os
import sys
import subprocess
import math
from scipy.stats import norm
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#==============================================================================
#Functions
#==============================================================================
def intersect(comparisons=None, raw_data=None, bedtools_intersect=None, 
                output_folder=None):
    print("Intersecting bed file replicates...")
    for sample_list in comparisons:
        for sample in sample_list:
            print("\t", sample)
            replicate_number = len(raw_data[sample])
            if replicate_number == 1:
                intersect_command = ("cat " + raw_data[sample][1][0] + " > " 
                                    + os.path.join(output_folder, 
                                    sample + ".bidir_predictions.intersect.bed"))
            else:
                intersect_command = (bedtools_intersect 
                                    + " -a " + raw_data[sample][1][0] 
                                    + " -b " + raw_data[sample][1][1])
                for i in range(replicate_number-2):
                    intersect_command = (intersect_command + " | " 
                                        + bedtools_intersect + " -a stdin -b " 
                                        + raw_data[sample][1][i])

                intersect_command = (intersect_command + " > " 
                                    + os.path.join(output_folder, 
                                    sample + ".bidir_predictions.intersect.bed"))

            os.system(intersect_command)

def merge(comparisons=None, bedtools_sort=None, bedtools_merge=None, 
            output_folder=None):
    print("Merging bed files across conditions...")
    for condition1, condition2 in comparisons:
        print('\t', condition1, 'vs.', condition2)
        merge_command = ("cat <(cat " + os.path.join(output_folder, 
                            condition1 + ".bidir_predictions.intersect.bed") 
                        + ") <(cat " + os.path.join(output_folder, 
                            condition2 + ".bidir_predictions.intersect.bed") 
                        + ") | " + bedtools_sort + " -i stdin | "
                        + bedtools_merge + " -i stdin > " 
                        + os.path.join(output_folder, 
                            condition1 + "_" + condition2 + ".merged.bed"))
        subprocess.call(['bash', '-c', merge_command])

def count(comparisons=None, output_folder=None, raw_data=None, 
            bedtools_multicov=None):
    print("Counting reads over merged regions using bam files...")
    for condition1, condition2 in comparisons:
        bed_file = os.path.join(output_folder, 
                condition1 + "_" + condition2 + ".merged.bed")
        bam1 = raw_data[condition1][0]
        bam2 = raw_data[condition2][0]
        print('\t', condition1, 'vs.', condition2)
        count_reads(bedfile=bed_file, bam1=bam1, bam2=bam2, 
                    tempdir=output_folder, label1=condition1, 
                    label2=condition2, bedtools_multicov=bedtools_multicov)

def count_reads(bedfile=None, bam1=None, bam2=None, tempdir=None, label1=None, 
                label2=None, bedtools_multicov=None):
    '''Counts reads across regions in a given bed file using bam files inputted
        by a user

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will 
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #This os.system call runs bedtools multicov to count reads in all specified
    #BAMs for given regions in BED
    os.system(bedtools_multicov + " -bams " + " ".join(bam1) + " " 
                + " ".join(bam2) + " -bed " + bedfile + " > " 
                + os.path.join(tempdir, 
                    label1 + "_" + label2 + ".count_file.bed"))

    #This section adds a header to the count_file and reformats it to remove 
    #excess information and add a column with the region for later use
    count_file = os.path.join(tempdir, 
                    label1 + "_" + label2 + ".count_file.header.bed")
    outfile = open(count_file, 'w')
    outfile.write("chrom\tstart\tstop\tregion\t" 
                    + '\t'.join([label1]*len(bam1)) + "\t" 
                    + '\t'.join([label2]*len(bam2)) + "\n")

    with open(os.path.join(tempdir, label1 + "_" + label2 + ".count_file.bed")) as F:
        for line in F:
            line = line.strip('\n').split('\t')
            chrom,start,stop = line[:3]
            counts = line[-(len(bam1)+len(bam2)):]
            outfile.write('\t'.join([chrom,start,stop]) + "\t" 
                            + chrom + ":" + start + "-" + stop + "\t"
                            + '\t'.join(counts) + "\n")

def DESeq2(comparisons=None, raw_data=None, output_folder=None):
    print("Performing DE-Seq2 analysis on comparisons...")
    for condition1, condition2 in comparisons:
        count_file = os.path.join(output_folder, 
                    condition1 + "_" + condition2 + ".count_file.header.bed")
        bam1 = raw_data[condition1][0]
        bam2 = raw_data[condition2][0]
        print('\t', condition1, 'vs.', condition2)
        write_deseq_script(bam1=bam1, bam2=bam2, tempdir=output_folder, 
                            count_file=count_file, label1=condition1, 
                            label2=condition2)
        os.system("R < " + os.path.join(output_folder, "DESeq.R") + " --no-save")

def write_deseq_script(bam1=None, bam2=None, tempdir=None, count_file=None, 
                        label1=None, label2=None):
    '''Writes an R script within the tempdir directory in TFEA output to run 
        either DE-Seq or DE-Seq2 depending on the number of user-inputted 
        replicates.

    Parameters
    ----------
    bedfile : string
        full path to a bed file containing full regions of interest which will
        be counted using bedtools multicov

    bam1 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    bam2 : list or array
        a list of full paths to bam files pertaining to a single condition 
        (i.e. replicates of a single treatment)

    tempdir : string
        full path to temp directory in output directory (created by TFEA)

    label1 : string
        the name of the treatment or condition corresponding to bam1 list

    label2 : string
        the name of the treatment or condition corresponding to bam2 list

    Returns
    -------
    None
    '''
    #If more than 1 replicate, use DE-Seq2
    if (len(bam1) > 1 and len(bam2) > 1):
        outfile = open(os.path.join(tempdir, 'DESeq.R'),'w')
        outfile.write('sink("' + os.path.join(tempdir,'DESeq.Rout') + '")\n')
        outfile.write('library("DESeq2")\n')
        outfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')
        outfile.write('countsTable <- subset(data, select=c('
                +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
                +'))\n')

        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- as.data.frame(c(' 
                        + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + '))\n')

        outfile.write('colnames(conds) <- c("treatment")\n')
        outfile.write('ddsFullCountTable <- DESeqDataSetFromMatrix(\
                                                    countData = countsTable, \
                                                    colData = conds, \
                                                    design = ~ treatment)\n')

        outfile.write('dds <- DESeq(ddsFullCountTable)\n')
        outfile.write('res1 <- results(dds,alpha = 0.05, \
                                        contrast=c("treatment",\
                                                        "'+label2+'",\
                                                        "'+label1+'"))\
                                                        \n')

        outfile.write('resShrink <- lfcShrink(dds, res = res1, \
                                                contrast = c("treatment",\
                                                "'+label2+'",\
                                                "'+label1+'"))\n')

        outfile.write('resShrink$fc <- 2^(resShrink$log2FoldChange)\n')
        outfile.write('res <- resShrink[c(1:3,7,4:6)]\n')
        outfile.write('write.table(res, file = "'
                        + os.path.join(tempdir, 
                            label1 + '_' + label2 + '.DESeq.res.txt') 
                        + '", append = FALSE, sep= "\t" )\n')
        outfile.write('sink()')
    else:
        outfile = open(os.path.join(tempdir, 'DESeq.R'),'w')
        outfile.write('sink("'+os.path.join(tempdir, 'DESeq.Rout') + '")\n')
        outfile.write('library("DESeq")\n')
        outfile.write('data <- read.delim("'+count_file+'", sep="\t", \
                        header=TRUE)\n')

        outfile.write('countsTable <- subset(data, select=c('
            +', '.join([str(i) for i in range(5,5+len(bam1)+len(bam2))])
            +'))\n')

        outfile.write('rownames(countsTable) <- data$region\n')
        outfile.write('conds <- c(' + ', '.join(['"'+label1+'"']*len(bam1)) 
                        + ', ' 
                        + ', '.join(['"'+label2+'"']*len(bam2)) 
                        + ')\n')

        outfile.write('cds <- newCountDataSet( countsTable, conds )\n')
        outfile.write('cds <- estimateSizeFactors( cds )\n')
        outfile.write('sizeFactors(cds)\n')                                                               
        outfile.write('cds <- estimateDispersions( cds ,method="blind", \
                        sharingMode="fit-only")\n')

        outfile.write('res <- nbinomTest( cds, "'+label1+'", "'+label2+'" )\n')
        outfile.write('rownames(res) <- res$id\n')                      
        outfile.write('write.table(res, file = "'
                        + os.path.join(tempdir,
                            label1 + '_' + label2 + '.DESeq.res.txt') 
                        + '", append = FALSE, sep= "\t" )\n')

        outfile.write('sink()')
    outfile.close()

def DESeq_to_bed(comparisons=None, output_folder=None):
    for condition1, condition2 in comparisons:
        deseq_file = os.path.join(output_folder, 
                        condition1 + '_' + condition2 + '.DESeq.res.txt')
        parse_deseq_file(deseq_file=deseq_file, condition1=condition1, 
                        condition2=condition2, output_folder=output_folder)

def parse_deseq_file(deseq_file=None, condition1=None, condition2=None, 
                    output_folder=None):
    output_file = os.path.join(output_folder, 
                    condition1 + '_' + condition2 + '.deseq.bed')
    with open(output_file,'w') as outfile:
        with open(deseq_file) as F:
            header = F.readline().strip('\n').split('\t')
            pval_index = header.index('"pvalue"')+1
            for line in F:
                region = line.strip('\n').split('\t')[0].strip('"')
                chrom = region.split(':')[0]
                start, stop = region.split(':')[1].split('-')
                pvalue = line.strip('\n').split('\t')[pval_index]
                if pvalue == 'NA':
                    pvalue = '1.0'
                outfile.write('\t'.join([chrom, start, stop, pvalue]) + '\n')

def pull_regions_from_comparisons(comparisons=None, padj_cutoff=False, 
                                    percentage_rank_cutoff=False, 
                                    output_folder=None,  padj_column=None,
                                    bedtools_intersect=None):
    print("Retrieving individual regions from merged bed file...")
    for condition_names in comparisons:
        condition1, condition2 = condition_names
        print('\t', condition1, 'vs.', condition2)
        bed_file = os.path.join(output_folder, 
                    condition1 + '_' + condition2 + '.deseq.bed')
        original_beds = (os.path.join(output_folder, 
                                        condition1 
                                        + ".bidir_predictions.intersect.bed"),
                        os.path.join(output_folder, 
                                        condition2 
                                        + ".bidir_predictions.intersect.bed"))
        pull_individual_regions(condition_names=condition_names, 
                                padj_cutoff=padj_cutoff, 
                                percentage_rank_cutoff=percentage_rank_cutoff, 
                                bed_file=bed_file, original_beds=original_beds,
                                padj_column=padj_column,
                                bedtoolsintersect=bedtools_intersect, 
                                output_folder=output_folder) 

def pull_individual_regions(condition_names=None, padj_cutoff=False, 
                            percentage_rank_cutoff=False, bed_file=None, 
                            original_beds=None, padj_column=None,
                            bedtoolsintersect=None, output_folder=None):
    '''Takes an intersected/merged bed file of regions and splits it up into
        regions originating from each specified condition.

        Parameters
        ----------
        condition_names : tuple
            a tuple of exactly two strings corresponding to the condition
            names for the comparison to be analyzed.
        bed_file : string
            full path to the intersected/merged bed file
        original_beds : tuple
            a tuple with exactly two strings that correspond to replicate 
            intersected bed files from the original conditions

    '''
    if padj_cutoff != False:
        outfilename = os.path.join(output_folder, 
                                    '_'.join(condition_names) 
                                    + '.differential_padj-' 
                                    + str(padj_cutoff) + '.bed')
    elif percentage_rank_cutoff != False:
        outfilename = os.path.join(output_folder, 
                                    '_'.join(condition_names) 
                                    + '.differential_percent-' 
                                    + str(percentage_rank_cutoff) + '.bed')
    else:
        sys.exit("One of padj_cutoff or percentage_rank_cutoff must not be False")
    outfile = open(outfilename,'w')
    if padj_cutoff == False:
        temp_file = list()

    with open(bed_file) as F:
        for line in F:
            linelist = line.strip('\n').split('\t')
            padj = float(linelist[padj_column])
            if padj_cutoff != False:
                if padj < padj_cutoff:
                    outfile.write(line)
            else:
                temp_file.append((padj, line))

    if padj_cutoff == False:
        for _, line in sorted(temp_file)[:int(len(temp_file)*percentage_rank_cutoff)]:
            outfile.write(line)

    outfile.close()

    outbeds = list()
    for i in range(len(original_beds)):
        bed = original_beds[i]
        label = condition_names[i]
        outbed = os.path.join(output_folder, 
                                '_'.join(condition_names) 
                                + '.' + label 
                                + '.differential.bed')
        command = (bedtoolsintersect + ' -wb -a ' + outfilename 
                    + ' -b ' + bed + ' > ' + outbed)
        subprocess.call(command, shell=True)
        outbeds.append(outbed)

def run_MDS_from_comparisons(comparisons=None, output_folder=None, 
                                MDS_src=None, fasta_file=None, PSSM_DB=None):
    print("Performing MDS analysis for each condition from given regions...")
    for condition_names in comparisons:
        for label in condition_names:
            print('\t', ' vs. '.join(condition_names), ':', label)
            bed = os.path.join(output_folder, 
                                    '_'.join(condition_names) 
                                    + '.' + label 
                                    + '.differential.bed')
            ID = '_'.join(condition_names) + '.' + label + '.MDS'
            run_MDS(MDS_src=MDS_src, bed=bed, fasta_file=fasta_file, 
                        PSSM_DB=PSSM_DB, output_folder=output_folder, ID=ID)

def run_MDS(MDS_src=None, bed=None, fasta_file=None, PSSM_DB=None, 
            output_folder=None, ID=None, bsn='150', H='1500', pv='0.000001'):
    MDS_command = [MDS_src, 'EVAL', 
                    '-bed', bed, 
                    '-fasta', fasta_file, 
                    '-DB', PSSM_DB, 
                    '-o', output_folder,
                    '-log_out', output_folder,
                    '-ID', ID,
                    '-bsn',bsn,
                    '-H', H,
                    '-pv', pv]
    subprocess.call(MDS_command)

def compute_MDS(histlist):
    '''Takes in a list of ints and returns h/H and H
    '''
    windowsize = 150
    H = sum(histlist)
    middle = int(len(histlist)/2)
    h = sum(histlist[middle-windowsize:middle]) + sum(histlist[middle:middle+windowsize])
    if H != 0:
        return float(h)/float(H), H
    else:
        return 0,0

def MDS(MDS1=None, MDS2=None):
    d = dict()
    with open(MDS1) as F:
        F.readline()
        for line in F:
            line = line.strip().split(',')
            # line = line.strip().split('\t')
            md,H = compute_MDS([int(x) for x in line[1:]])
            d[line[0]] = [md,H]

    with open(MDS2) as F:
        F.readline()
        for line in F:
            line = line.strip().split(',')
            # line = line.strip().split('\t')
            md,H = compute_MDS([int(x) for x in line[1:]])
            d[line[0]].append(md)
            d[line[0]].append(H)
    

    X = [0] * len(d)
    Y = [0] * len(d)
    pvals = [0] * len(d)
    genelist = [''] * len(d)
    diff = list()
    for key in d:
        mdk=d[key][0]
        mdj=d[key][2]
        diff.append(float(mdj)-float(mdk))
    mean = sum(diff)/len(diff)
    i = 0
    for key in d:
        mdk=float(d[key][0])
        mdj=float(d[key][2])
        Nk=float(d[key][1])
        Nj=float(d[key][3])
        if Nj > 10 and Nk > 10:
            p=((mdj*Nj)+(mdk*Nk))/(Nj+Nk)
            SE=(p*(1-p))*((1/Nj)+(1/Nk))
            Y[i] = mdj-mdk-mean
            X[i] = math.log((Nj+Nk)/2.0,10)
            genelist[i] = key.split('.')[0].split('_')[0]+'_'+key.split('.')[-2]
            try:
                z = (mdj-mdk-mean)/math.sqrt(SE)
            except ZeroDivisionError:
                z=0
            cdf=norm.cdf(z,0,1)
            p=min(cdf,1-cdf)*2
            pvals[i] = p
        i += 1
        
    return X, Y, pvals, genelist

def write_MDS_output(genelist=None, X=None, Y=None ,p_vals=None, 
                        output=None):
    sorted_data = [(p,x,y,g) for p,x,y,g in sorted(zip(p_vals,X,Y,genelist)) if len(g) > 1]
    with open(output, 'w') as outfile:
        outfile.write('TF\tLog10Events\tMDS_difference\tp-value\n')
        for i in range(len(sorted_data)):
            outfile.write('\t'.join([sorted_data[i][3],
                                    str(sorted_data[i][1]),
                                    str(sorted_data[i][2]),
                                    str(sorted_data[i][0])]) + '\n')

def parse_MDS_output(input_file=None):
    genelist = list()
    p_vals = list()
    x = list()
    y = list()
    with open(input_file) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            genelist.append(line[0])
            x.append(float(line[1]))
            y.append(float(line[2]))
            p_vals.append(float(line[3]))
    
    return x, y, genelist, p_vals

def print_MDS(condition1=None, condition2=None, output_folder=None, 
                pval_cut=None, n=None, motif=None):
    input_file = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.MDS_diff.txt')
    x, y, genelist, p_vals = parse_MDS_output(input_file=input_file)
    print('TF\tLog10Events\tMDS_difference\tp-value')
    if pval_cut != None:
        for i in range(len(p_vals)):
            p = p_vals[i]
            if p < pval_cut:
                print(genelist[i], x[i], y[i], p_vals[i])
    if n != None:
        for i in range(len(p_vals)):
            if i <= n:
                print(genelist[i], x[i], y[i], p_vals[i])
    
    if motif != None:
        for i in range(len(genelist)):
            if motif in genelist[i]:
                print(genelist[i], x[i], y[i], p_vals[i])

def plot_MA(condition1=None, condition2=None, output_folder=None, ax=None, 
            pval_cut=None, label=False):
    input_file = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.MDS_diff.txt')
    x, y, genelist, p_vals = parse_MDS_output(input_file=input_file)
    upx = [x1 for x1,y1,p in zip(x,y,p_vals) if p < pval_cut and y1 > 0]
    upy = [y1 for x1,y1,p in zip(x,y,p_vals) if p < pval_cut and y1 > 0]
    dnx = [x1 for x1,y1,p in zip(x,y,p_vals) if p < pval_cut and y1 < 0]
    dny = [y1 for x1,y1,p in zip(x,y,p_vals) if p < pval_cut and y1 < 0]
    ax.scatter(x,y, c="k", edgecolor="", s=30)
    ax.scatter(upx, upy, c="r", edgecolor="", s=30, 
                label="Up: p < " + str(pval_cut))
    ax.scatter(dnx, dny, c="g", edgecolor="", s=30, 
                label="Down: p < " + str(pval_cut))
    if label != False:
        lbx = [x1 for x1,y1,g in zip(x,y,genelist) if label in g]
        lby = [y1 for x1,y1,g in zip(x,y,genelist) if label in g]
        ax.scatter(lbx, lby, c="orange", edgecolor="", s=30, 
                    label=label + " motifs")
    ax.set_title(condition1 + " vs. " + condition2)
    ax.set_ylabel('MD Score Difference')
    ax.set_xlabel('Mean Overlap Events (log10)')
    plt.legend(loc='lower right')

def label(condition1=None, condition2=None, output_folder=None, name=None, 
            label=None, xytext=None, ax=None):
    input_file = os.path.join(output_folder, 
                            '_'.join([condition1, condition2])
                            + '.MDS_diff.txt')
    x, y, genelist, p_vals = parse_MDS_output(input_file=input_file)
    for i, txt in enumerate(genelist):
        if name in txt:
            ax.annotate(label, xy=(x[i], y[i]), xytext=xytext, 
                            arrowprops=dict(facecolor='black', shrinkA=0.05, 
                                            arrowstyle='-|>'))


def plot_heatmap(input_file=None, motif=None, ax=None):
    with open(input_file) as file1:
        for line in file1:
            line = line.strip('\n').split(',')
            TF = line[0]
            if TF == motif:
                #Data in the form of a histogram with motif counts per bp
                x = [float(x) for x in line[1:]]

                window_size = len(x)/2

                bins=100
                edges = np.linspace(-window_size, window_size, bins+1)
                
                counts = np.zeros(bins)
                window = int(len(x)/bins)
                for i in range(bins):
                    j = i*window
                    average = sum(x[j:j+window])
                    counts[i] = average

                # counts, edges = np.histogram(x, bins=bins)
                edges        = (edges[1:]+edges[:-1])/2. 
                ax.bar(edges, counts, width=(edges[-1]-edges[0])/bins)
            
                norm    = mpl.colors.Normalize(vmin=min(counts), vmax=max(counts))
                cmap    = cm.Blues
                m       = cm.ScalarMappable(norm=norm, cmap=cmap)
                colors  = [m.to_rgba(c) for c in counts] 

                # print(colors)
                
                ax.bar(edges, np.ones((len(edges),)), color=colors, 
                        width=(edges[-1]-edges[0])/len(edges),
                        edgecolor=colors)
                ax.set_xlim(-1500, 1500)
                ax.set_ylim(0, 1)
                ax.set_title(input_file.split('.')[-3] 
                                + " (n=" + str(sum(x)) + ")")
                ax.tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='on') # labels along the bottom edge are off
                ax.tick_params(
                    axis='y',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    right='off',      # ticks along the bottom edge are off
                    left='off',         # ticks along the top edge are off
                    labelbottom='off',
                    labelleft='off') 
                    


if __name__ == "__main__":
    import matplotlib.gridspec as gridspec
    import numpy as np

    F = plt.figure(figsize=(15,15))
    outer = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.2)

    inner = gridspec.GridSpecFromSubplotSpec(2, 1,
                        subplot_spec=outer[0], wspace=0.1, hspace=0.1)

    ax = plt.Subplot(F, inner[0])
    F.add_subplot(ax)
    x = range(-1500,1500)
                
    bins=100
    # counts, edges = np.histogram(x, bins=bins)
    edges = np.linspace(-1500,1500,101)
    counts = range(100)
    edges        = (edges[1:]+edges[:-1])/2. 
    ax.bar(edges, counts, width=(edges[-1] - edges[0])/bins  )

    norm    = mpl.colors.Normalize(vmin=min(counts), vmax=max(counts))
    cmap    = cm.Blues
    m       = cm.ScalarMappable(norm=norm, cmap=cmap)
    colors  = [m.to_rgba(c) for c in counts] 
    
    ax.bar(edges,np.ones((len(edges),)), color=colors, 
            width=(edges[-1]-edges[0])/len(edges),
            edgecolor=colors )
    ax.set_xlim(-1500, 1500)
    ax.set_ylim(0, 1)

    plt.show()
