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
import matplotlib.pyplot as plt
#==============================================================================
#Functions
#==============================================================================
def intersect_merge_bed(bed1=None, bed2=None, tempdir=None, 
                        bedtools_intersect=None, bedtools_sort=None,
                        bedtools_merge=None):
    '''Takes in two lists of bed files, each containing replicates for one 
        condition. Intersects replicates using bedtools then merges intersected
        regions.

    Parameters
    ----------
    bed1 : list or array
        full paths to bed files for condition 1 (strings)
    
    bed2 : list or array
        full paths to bed files for condition 2 (strings)
        
    tempdir : string
        full path to tempdir directory in output directory (created by TFEA)

    Returns
    -------
    combined_input_merged_bed : string 
        full path to a bed file containing the merged regions inputted by the 
        user 
    '''
    #Define the output file
    combined_input_merged_bed = os.path.join(tempdir, 
                                                "combined_input.merge.bed")

    if len(bed1) > 1:
        #Build command to perform bedtools intersect on condition1 beds
        intersect1 = (bedtools_intersect + " -a " + bed1[0] + " -b " + bed1[1])
        for bedfile in bed1[2:]:
            intersect1 = (intersect1 + " | " + bedtools_intersect 
                        + " -a stdin -b " + bedfile)
    else:
        intersect1 = "cat " + bed1[0]

    if len(bed2) > 1:
        #Build command to perform bedtools intersect on condition2 beds
        intersect2 = (bedtools_intersect + " -a " + bed2[0] + " -b " + bed2[1])
        for bedfile in bed2[2:]:
            intersect2 = (intersect2 + " | " + bedtools_intersect 
                        + " -a stdin -b " + bedfile)
    else:
        intersect2 = "cat " + bed2[0]

    #Build full command which pipes both intersect commands into cat, then 
    # sorts and merges this resulting bed file
    command = ("cat <(" + intersect1 + ") <(" + intersect2 
                + ") | "+bedtools_sort+" -i stdin | " + bedtools_merge 
                + " -i stdin > " + combined_input_merged_bed)
    
    #Need to use subprocess here because this command is bash not shell
    subprocess.call(['bash', '-c', command])

    return combined_input_merged_bed
#==============================================================================

#==============================================================================
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
    outfile = open(count_file,'w')
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

    return count_file
#==============================================================================

#==============================================================================
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
#==============================================================================

#==============================================================================
def pull_individual_regions(condition_names=None, padj_cutoff=False, 
                            percentage_rank_cutoff=False, bed_file=None, 
                            original_beds=None, padj_column=-3,
                            bedtools_intersect=None):
    if padj_cutoff != False:
        outfilename = './differential_padj-' + str(padj_cutoff) + '.bed'
    elif percentage_rank_cutoff != False:
        outfilename = ('./differential_percent-' + str(percentage_rank_cutoff) 
                        + '.bed')
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
        outbed = os.path.join(src_path, label + '_differential.bed')
        command = bedtools_intersect + ' -wb -a ' + outfilename 
                    + ' -b ' + bed + ' > ' + outbed
        subprocess.call(command, shell=True)
        outbeds.append(outbed)
        
    return outbeds
#==============================================================================

#==============================================================================
def parse_deseq_file(deseq_file=None, condition1=None, condition2=None, 
                    output_folder=None):
    output_file = os.path.join(output_folder, condition1 + '_' + condition2 + '.deseq.bed')
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
#==============================================================================

#==============================================================================


#Initialization module from MDS_Differentially_Transcribed
#==============================================================================
#Paths
#==============================================================================


#==============================================================================
#Functions
#==============================================================================
def pull_individual_regions(condition_names=None, padj_cutoff=False, percentage_rank_cutoff=False, 
                            bed_file=None, original_beds=None, padj_column=-1,
                            bedtools_intersect=None):
    if padj_cutoff != False:
        outfilename = './differential_padj-' + str(padj_cutoff) + '.bed'
    elif percentage_rank_cutoff != False:
        outfilename = './differential_percent-' + str(percentage_rank_cutoff) + '.bed'
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
        outbed = os.path.join(src_path, label + '_differential.bed')
        command = bedtoolsintersect + ' -wb -a ' + outfilename + ' -b ' + bed + ' > ' + outbed
        subprocess.call(command, shell=True)
        outbeds.append(outbed)
        
    return outbeds
#====================================================================================================================

#====================================================================================================================
def run_MDS(MDS_src=MDS_src, bed=None, fasta_file=None, PSSM_DB=None, src_path=src_path, ID=None, 
            bsn='150', H='1500', pv='0.000001'):
    MDS_command=[MDS_src, 'EVAL', 
             '-bed', bed, 
             '-fasta', fasta_file, 
             '-DB', PSSM_DB, 
             '-o', src_path,
             '-log_out', src_path,
             '-ID', ID,
             '-bsn',bsn,
             '-H', H,
             '-pv', pv]
    subprocess.call(MDS_command)
#====================================================================================================================

#====================================================================================================================
#Functions associated with MDS calculation and plotting
def plot_MA(X=None, Y=None, pvals=None, genelist=None, annotate=True, TF=None, 
            name1=None, name2=None, figure=None, subplot=None, title=None, pval_cut=0.1):
    ax = figure.add_subplot(subplot)
    ax.scatter(X, Y, c="k", edgecolor="", s=14) 
    for i, txt in enumerate(genelist):
        if TF != None and TF in txt:
            ax.annotate(txt, (X[i], Y[i]))
            
    sig1 = [x for x,p in zip(X,pvals) if p < pval_cut and x > 0]
    sig2 = [y for x,y,p in zip(X,Y,pvals) if p < pval_cut and x > 0]
    ax.scatter(sig1,sig2,c='r',edgecolor="",s=20, label="Up: p < 0.1")
    
    sig3 = [x for x,p in zip(X,pvals) if p < pval_cut and x < 0]
    sig4 = [y for x,y,p in zip(X,Y,pvals) if p < pval_cut and x < 0]
    ax.scatter(sig3,sig4,c='g',edgecolor="",s=20, label="Down: p < 0.1")
    
    statx = [x for x,g in zip(X,genelist) if 'STAT' in g or 'STA5B' in g]
    staty = [y for y,g in zip(Y,genelist) if 'STAT' in g or 'STA5B' in g]
    ax.scatter(statx, staty, facecolor="None", edgecolor="",linewidth='5', s=20, label='STAT motifs')
    
    if annotate:
        for i, txt in enumerate(annotate_list):
                ax.annotate(txt, (sig1[i], sig2[i]))

    if annotate:
        for i, txt in enumerate(annotate_list):
                ax.annotate(txt, (sig3[i], sig4[i]))
                


    ax.set_title(title)
    ax.set_ylabel('MD Score Difference')
    ax.set_xlabel('Mean Overlap Events (log10)')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.legend(loc='lower right')

    # plt.savefig(savedir + 'MA_plot_' + name1 + '-' + name2 + '.png')
    # plt.savefig(savedir + 'MA_plot_intersect_30_30_CA.png')
#====================================================================================================================

#====================================================================================================================
#Takes in a list of ints and returns h/H and H
def compute_MDS(histlist):
    windowsize = 150
    H = sum(histlist)
    middle = int(len(histlist)/2)
    h = sum(histlist[middle-windowsize:middle]) + sum(histlist[middle:middle+windowsize])
    if H != 0:
        return float(h)/float(H), H
    else:
        return 0,0
#====================================================================================================================

#====================================================================================================================
def MDS(MDS1=None,MDS2=None,print_statements=False):
    name1 = MDS1.split('/')[-1].split('_')[0]
    name2 = MDS2.split('/')[-1].split('_')[0]
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
            if p < 0.1:
                if mdj-mdk-mean > 0:
                    if print_statements:
                        print('up', key, math.log((Nj+Nk)/2.0,10), mdj-mdk-mean)
                else:
                    if print_statements:
                        print('down', key, math.log((Nj+Nk)/2.0,10), mdj-mdk-mean)
        i += 1
        
    return X, Y, pvals, genelist
#====================================================================================================================

#====================================================================================================================
def write_MDS_output(genelist=None, X=None,Y=None,p_vals=None,output=None):
    sorted_data = [(p,x,y,g) for p,x,y,g in sorted(zip(p_vals,X,Y,genelist))]
    with open(output, 'w') as outfile:
        outfile.write('TF\tMDS_difference\tLog10Events\tp-value\n')
        for i in range(len(X)):
            outfile.write('\t'.join(sorted_data[3][i],
                                    str(sorted_data[1][i]),
                                    str(sorted_data[2][i]),
                                    str(sorted_data[0][i])) + '\n')
#====================================================================================================================

#====================================================================================================================
def parse_MDS_output(output=None):
    genelist = list()
    p_vals = list()
    x = list()
    y = list()
    with open(output) as F:
        F.readline()
        for line in F:
            line = line.strip('\n').split('\t')
            genelist.append(line[0])
            x.append(float(line[1]))
            y.append(float(line[2]))
            p_vals.append(float(line[3]))
    
    return x,y,genelist,p_vals
#====================================================================================================================

#====================================================================================================================
def differential_MDS(condition_names=None, padj_cutoff=None, percentage_rank_cutoff=None, 
                    bed_file=None, original_beds=None, fasta_file=None, PSSM_DB=None, src_path=src_path,
                    subplot=None, figure=None):
    
    outbeds = pull_individual_regions(condition_names=condition_names, padj_cutoff=padj_cutoff, 
                                  percentage_rank_cutoff=percentage_rank_cutoff, bed_file=bed_file, original_beds=original_beds)
    
    for ID, bed in zip(condition_names,outbeds):
        run_MDS(bed=bed, fasta_file=fasta_file, PSSM_DB=PSSM_DB, src_path=src_path, ID=ID)
    
    MDS1 = os.path.join(src_path, condition_names[0] + '_MDS.csv')
    MDS2 = os.path.join(src_path, condition_names[1] + '_MDS.csv')
    X, Y, pvals, genelist = MDS(MDS1=MDS1,MDS2=MDS2,print_statements=False)
    
    if padj_cutoff != False:
        title = condition_names[0] + " vs. " + condition_names[1] + " p-val < " + str(padj_cutoff)
    elif percentage_rank_cutoff != False:
        title = condition_names[0] + " vs. " + condition_names[1] + " percent = " + str(percentage_rank_cutoff)
    
    output = os.path.join(src_path, condition_names[0] + "_" + condition_names[1] + ".mds.txt")
    write_MDS_output(X=X, Y=Y, p_vals=pvals, genelist=genelist, output=output)
