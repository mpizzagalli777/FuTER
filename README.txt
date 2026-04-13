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
