FuTER (Fusion TE Reporter) was developed to allow users to efficiently identify fusion transcripts arising from Transposable Elements (TEs).
Although it was developed using human samples, it should be easily adapted to any kind of model organism (although additional work by the end user will be required).

To provide an easier start up, a set of the required input files is provided for Gencode version 44 in the HG38_Genome_Data directory. 
More information on how each file was created can be found in the README file found in that directory.  
In order to perform this analysis on a different genome/model organism, these input files must be generated/obtained. 

There are two main steps in the pipeline. 
The first generates the reference libraries that will be used by the program to identify the fusion transcripts. This will generate a new folder within FuTER called TEIF (an ode to the CTAT Viral Integration Finder). It will contain the combined reference genome and TE scaffolds as well as the final TE database that is filtered to remove TEs that are very similar to one another. If users are utilizing HG38, we have made the input data used below available at: 

To run this step, the TE_lib_integration.py must be run as follows:

python TE-lib-integration.py \
	--TE_db_fasta $datadir/HG38_Genome_Data/RM_merged_TEs_gencode.v44_sequences.fa \
	--genome_fasta $datadir/HG38_Genome_Data/gencode.v44.genome.fa \
	--genome_gtf $datadir/HG38_Genome_Data/gencode.v44.annotation.gtf \
	--repeatmasker_gtf $datadir/HG38_Genome_Data/RM_merged_TEs_gencode.v44.gtf \
	--TE_splice_acceptor $datadir/HG38_Genome_Data/RM_merged_TEs_gencode.v44_acceptor_predictions_0.5.gtf \
	--TE_splice_donor $datadir/HG38_Genome_Data/RM_merged_TEs_gencode.v44_donor_predictions_0.5.gtf \



The second performs the actual identification of fusion transcripts. 
To run this step, the script must be run as follows:

python $FUTER_basedir/FuTER_Pipeline.py  \
    --reads_fastq $FUTER_basedir/test_reads/RM_testset_newTEs_AllFusions_and_Transcripts_sim_reads.fasta \
    --genome_fasta $FUTER_basedir/HG38_Genome_Data/gencode.v44.genome.fa \
    --genome_gtf   $FUTER_basedir/HG38_Genome_Data/gencode.v44.annotation.gtf \
    --outdir outdir/data \
    --CPU N

This will create a series of folders within the output directory containing the outputs of intermediate processing steps as well as the final output of the pipeline. Contents are as follows:
1) Phase1_init_alignments - the results of aligning long reads to the user defined TE genome. These are broken down into three groups (Pure TE - reads that aligned 100% to a TE sequence; chimeric - reads that partially aligned to a TE sequence, as well as the human genome; No_TE - reads that aligned less than 5% to a TE sequence. These reads are discarded and do not proceed through the pipeline as they are not potential contributors of major TE chimeric transcripts). 
2) Phase2_LR_fusion_cand - the output of LR reads aligned to the combined TE / Human Genome. This identifies the reads that were identified by Ctat minimap2 as being chimeric reads which in this case aligned to both gene-encoding regions of the human genome as well as TE scaffolds. The "_chims_described" file serves as the reference document for the next Phase of the pipeline as it identiifes and organizes the reads that aligned to the two areas of the genome (gene-encoding regions and TE scaffold). 
3) Phase3_TELR_Fusion - This folder contains the majority of output files from the pipeline. This phase serially filters the identified chimeric reads based on mapping quality, promiscuity of read alignment, etc (most settings can be adjusted by the user). 
4) Final_output - Congrats. This contains the output! 

In order to create visualizations, we also provide two further scripts that can be used to further process the long read sequencing data and generate IGV based visualizations, as well as the input data to create GViz based figures as are present in the original manuscript. These two scripts are 1) LR_IGV_report_gen.py and 2) LR_contig_combination. 

LR_IGV_report_gen.py will generate an IGV-Reports HTML for a selected set of fusions. It will, for a specified set of Fusions of Interest (provided in the --fusion_list variable) extract the reads that were used to identify the fusion and generate an IGV-Reports file demonstrating alignment patterns. This can be used to analyze the breakpoints and investigate the structure of the fusions. It is run as follows:

python $LRTE_basedir/scripts/LR_IGV_report_gen.py \
  --Phase3_TELR_Fusion_dir $out_basedir/Phase3_TELR_Fusion/ \
  --Sample_prefix       GB2_p20_PCS114_20241114_py_unmasked \
  --fusion_list  $out_basedir/GB11_Fusions_of_Interest_all.txt \
  --max_LR_per_fusion N \
  --threads N \
  --genome_lib_dir $LRTE_basedir/TEgenome_prep \
  --util_dir       $FUTER_basedir/scripts/CTAT-LR-fusion.v1.1.1/util \
  --fi_util_dir    $FUTER_basedir/scripts/CTAT-LR-fusion.v1.1.1/FusionInspector/util \
  --basedir         $out_basedir/

Phase3_TELR_Fusion_dir
Sample_prefix
fusion_list
max_LR_per_fusion
threads
genome_lib_dir
util_dir

The next step, LR_contig_combination.py, analyses the reads in order to identify whether the reads contain a significant proportion of internal promoter sequences. This can be a sign of a fusion being a technical artifact arising from the PCR steps used in the ONT protocol. After this filtering, the script will generate consensus sequences for the fusion transcripts using both Medaka and RNA Bloom are used. These approaches polish the more error prone ONT reads and RNA bloom uses a scaffold free approach to generate a consensus sequence for the fusion transcript. The script will also align reads to the consensus sequences in order to provide the bam files necessary to create the visualizations found in the paper. LR_contig_combination.py is run as follows:

python $LRTE_basedir/scripts/LR_contig_combination.py \
  --base_dir \
  --TEgenome_dir $FUTER_basedir/TEgenome_prep/ \
  --fastq_dir [Path to long read fastq files] \
  --Script_dir $FUTER_basedir/scripts \
  --TEIF_dir $FUTER_basedir/TEIF \
  --genome_fasta $FUTER_basedir/TEgenome_prep/gencode.v44.genome.fa \
  --fusion_list [path to]/GB11_Fusions_of_Interest_all.txt \
  --Sample_prefix GB11 \
  --min_num_LR [Min number of long reads required to proceed] \
  --threads N \
  --rnabloom_extra [Options passed to RNA Bloom]

base_dir - 
TEgenome_dir - 
fastq_dir - Path to long read fastq files
Script_dir - 
TEIF_dir - 
genome_fasta - 
fusion_list - see above. 
Sample_prefix - 
min_num_LR - Min number of long reads required to proceed with analysis. Higher read counts ensure better coverage 
threads - Number of CPUs provided 
rnabloom_extra - Options passed to RNA Bloom
