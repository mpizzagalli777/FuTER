#!/usr/bin/env python3

import sys, os, re
import subprocess
import argparse
import logging
import pandas as pd
from Bio import SeqIO
from pathlib import Path


sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "utils/PyLib"])
    )
from Pipeliner import Pipeliner, Command

FORMAT = "%(asctime)-15s: %(levelname)s %(module)s.%(name)s.%(funcName)s %(message)s"
logger = logging.getLogger(__file__)
logging.basicConfig(stream=sys.stderr, format=FORMAT, level=logging.INFO)

def prep_minimap2_reference(genome_fa, mm2_db_dir, mm2_db_name, mm2_splice_file, ref_gtf, which, ref_splice_bed=None):

    chk = Path(f"{mm2_db_dir}/__mm2_prep_chkpts/{which}.build.ok")
    if chk.exists():
        logger.info("Skipping prep ref: found %s", chk)
        return

    pip = Pipeliner(checkpoint_dir=os.path.join(mm2_db_dir, "__mm2_prep_chkpts"))
    ctat_mm2_dir = os.path.join(os.path.dirname(__file__), "utils", "ctat-minimap2")

    # index genome
    cmd = f"{ctat_mm2_dir}/ctat-minimap2 -d {mm2_db_name} {genome_fa}"
    pip.add_commands([ Command(cmd, f"{which}_mm2_prep_genome.ok") ])
    
    # Making the splice files
    if which == "ref":
        cmd = (
            f"{ctat_mm2_dir}/misc/paftools.ctat.js gff2bed "
            f"{ref_gtf} > {mm2_splice_file}"
        )
        pip.add_commands([ Command(cmd, "make_ref_splice.ok") ])
        
    elif which == "comb":
        # use the already created ref splice as the starting point
        if not ref_splice_bed:
            raise ValueError("you must pass ref_splice_bed when which=='comb'")
        cmd = f"cp {ref_splice_bed} {mm2_splice_file}"
        pip.add_commands([ Command(cmd, "copy_ref_splice.ok") ])
        
        cmd = f"grep '^chrTE_' {ref_gtf} > {mm2_db_dir}/te_only.gtf"
        pip.add_commands([ Command(cmd, "dump_te_gtf.ok") ])
        
        cmd = (
            f"{ctat_mm2_dir}/misc/paftools.ctat.js gff2bed "
            f"{mm2_db_dir}/te_only.gtf >> {mm2_splice_file}"
        )    
        pip.add_commands([ Command(cmd, "append_TE_splice.ok") ])

    elif which == "comb_un":
        # use the already created ref splice as the starting point
        if not ref_splice_bed:
            raise ValueError("you must pass ref_splice_bed when which=='comb2'")
        cmd = f"cp {ref_splice_bed} {mm2_splice_file}"
        pip.add_commands([ Command(cmd, "copy_ref_splice2.ok") ])
        
        cmd = (
            f"{ctat_mm2_dir}/misc/paftools.ctat.js gff2bed "
            f"{mm2_db_dir}/te_only.gtf >> {mm2_splice_file}"
        )    
        pip.add_commands([ Command(cmd, "append_TE_splice2.ok") ])        

    pip.run()
    chk.touch()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--TE_db_fasta", dest="TE_db_fasta",
                        type=str, required=True,
                        help="Full path to TE database (fasta file)")

    parser.add_argument("--genome_fasta", dest="ref_genome_fasta",
                        type=str, required=True,
                        help="Full path to genome fasta file")
    
    parser.add_argument("--genome_gtf", dest="ref_gtf",
                        type=str, required=True,
                        help="Full path to genome gtf")
    
    parser.add_argument("--repeatmasker_gtf", dest="repeatmasker_gtf",
                        type=str, required=True,
                        help="Full path to the output gtf from repeat masker and One Code for the genome")
                        
    parser.add_argument("--CPU",
                        type=int,
                        default=4,
                        help="num threads to use")
                                                
    parser.add_argument("--TE_splice_acceptor", dest="splice_acceptor_gtf",
                        type=str, required=True,
                        help="SpliceAI acceptor-style GTF for TEs")
                        
    parser.add_argument("--TE_splice_donor",    dest="splice_donor_gtf",
                        type=str, required=True,
                        help="SpliceAI donor-style GTF for TEs")
                        
    args=parser.parse_args()
    
    TE_db_fasta      = args.TE_db_fasta
    RM_gtf			 = args.repeatmasker_gtf
    ref_genome_fasta = args.ref_genome_fasta
    ref_gtf          = args.ref_gtf
    num_threads      = args.CPU
    homedir 		 = os.path.dirname(os.path.realpath(__file__))
    	
    if not os.path.exists(ref_genome_fasta):
        exit("Error, not finding genome at: {}".format(ref_genome_fasta))

    if not os.path.exists(ref_gtf):
        exit("Error, not finding ref annotation at: {}".format(ref_gtf))


    TEIF_dir = os.path.join(homedir, "TEIF")
    if not os.path.exists(TEIF_dir):
        os.makedirs(TEIF_dir)

    checkpoints_dir = f"{TEIF_dir}/__checkpts.dir"
    pipeliner = Pipeliner(checkpoints_dir)
    
    # copy the TE db to the TEIF dir and index it.
    installed_TE_db = os.path.join(TEIF_dir, "TE_db.fasta")
    filt_TE_db = os.path.join(TEIF_dir, "TE_db_filt.fasta")
    logger.info("-installing TE db")
    pipeliner.add_commands([Command(f"cp {TE_db_fasta} {installed_TE_db}", "cp_TE_to_TEIF.ok")])
    pipeliner.run()
    
    pipeliner.add_commands([Command(f"cd-hit-est -i {installed_TE_db} -o {filt_TE_db} -c .98 -d 0 -M 4000 -T 0 ", "cd-hit-est.ok")])
    pipeliner.run()
    
    reheadered_TE_db = os.path.join(TEIF_dir, "TE_db_filt_renamed.fa")
    mapping = {}
    with open(reheadered_TE_db, "w") as outfh:
        for i, rec in enumerate(SeqIO.parse(filt_TE_db, "fasta"), start=1):
            new_name = f"chrTE_{i}"
            mapping[new_name] = rec.id
            rec.id = new_name
            rec.description = ""       # drop the old description
            SeqIO.write(rec, outfh, "fasta")
            
    inv_map = {orig: te_chr for te_chr, orig in mapping.items()}
    
    te_breaks = {}
    for gtf_path in (args.splice_acceptor_gtf, args.splice_donor_gtf):
         with open(gtf_path) as gf:
             for line in gf:
                 if line.startswith("#") or not line.strip():
                     continue
                 cols = line.split("\t")
                 orig_id = cols[0]   
                 pos     = int(cols[3])
                 te_chr  = inv_map.get(orig_id)
                 if not te_chr:
                     continue
                 te_breaks.setdefault(te_chr, set()).add(pos)
    
    pipeliner.add_commands([Command(f"samtools faidx {reheadered_TE_db}", "faidx_TEdb.ok")])
    pipeliner.run()    

    
    # prevent any exon masking (and expand exon regions by 10 bp on either side to provide a buffer)
    exons_bed = f"{homedir}/HG38_Genome_Data/gencode.v44.annotation_exons.bed"
    exons_exp_bed = f"{homedir}/HG38_Genome_Data/gencode.v44.annotation_exons_exp.bed"    
    cmd = (
        f"bedtools slop "
        f"-i {exons_bed} "
        f"-g {homedir}/HG38_Genome_Data/gencode.v44.genome.chrom.sizes -b 10 "
        f"> {exons_exp_bed}"
    )
    pipeliner.add_commands([ Command(cmd, "expand_exons.ok") ])
    

    #This section masks all non exonic TEs identified by RepeatMasker 
    RM_gtf_nonexonic = f"{TEIF_dir}/RM_nonexonic_TEs.gtf"
    cmd = (
        f"bedtools subtract "
        f"-a {RM_gtf} "
        f"-b {exons_exp_bed} "
        f"> {RM_gtf_nonexonic}"
    )
    pipeliner.add_commands([ Command(cmd, "subtract_exons_RM.ok") ])    
    
    
    # do the masking using bedtools
    TE_masked_genome = f"{TEIF_dir}/ref_genome.TE_masked.fa"
    cmdstr = f"maskFastaFromBed -fi {ref_genome_fasta} -fo {TE_masked_genome} -bed {RM_gtf_nonexonic}"
    pipeliner.add_commands([Command(cmdstr, "TEmaskingrefgenome.ok")])
    pipeliner.run()

	# before masking: 161,334,628
	# after RM masking: 1,729,240,559
    # after blast masking: 1,731,177,837
    # so, lots of new bases masked but minimal from additional layer that adds significant time.
    
    # create new fasta file including TEes and human genome together: 
    combined_genomes_fa = os.path.join(TEIF_dir, "ref_genome_plus_TE.fa")
    logger.info(f"-combining {TE_masked_genome} and {reheadered_TE_db} into {combined_genomes_fa}")
    cmdstr = f"cat {TE_masked_genome} {reheadered_TE_db} > {combined_genomes_fa}"
    pipeliner.add_commands([Command(cmdstr, "combineMaskedGenomeWithTEs.ok")])
    pipeliner.run()
        
    cmdstr = f"samtools faidx {combined_genomes_fa}"
    pipeliner.add_commands([Command(cmdstr, "combinedgenomes.faidx.ok")])
    pipeliner.run()
        
    combined_genomes_fa_unmasked = os.path.join(TEIF_dir, "ref_genome_unmasked_plus_TE.fa")
    logger.info(f"-combining {ref_genome_fasta} and {reheadered_TE_db} into {combined_genomes_fa_unmasked}")
    cmdstr = f"cat {ref_genome_fasta} {reheadered_TE_db} > {combined_genomes_fa_unmasked}"
    pipeliner.add_commands([Command(cmdstr, "combineUnmaskedGenomeWithTEs.ok")])
    pipeliner.run()
            
    cmdstr = f"samtools faidx {combined_genomes_fa_unmasked}"
    pipeliner.add_commands([Command(cmdstr, "combinedgenomes_unmask.faidx.ok")])
    pipeliner.run()
    
    # make the new masked genome gtf. Will just use the name of the TE as gene name and have them be one gene/transcript
    combined_gtf = os.path.join(TEIF_dir, "ref_genome_plus_TE.gtf")
    logger.info(f"-writing combined GTF to {combined_gtf}")
    with open(combined_gtf, "w") as out:
        with open(ref_gtf) as orig:
            for line in orig:
                out.write(line)

        # 2) append one gene/transcript plus exon(s) per TE
        for rec in SeqIO.parse(reheadered_TE_db, "fasta"):
             te_chr = rec.id
             orig_id = mapping[te_chr]
             length = len(rec.seq)
             # gene + transcript attrs
             gene_attr = f'gene_id "{te_chr}"; gene_name "{orig_id}"; gene_biotype "transposable_element";\n'
             tx_attr   = f'gene_id "{te_chr}"; transcript_id "{orig_id}.t1";\n'
 
             #  a) write gene & transcript
             out.write(f"{te_chr}\tTE_library\tgene\t1\t{length}\t.\t+\t.\t{gene_attr}")
             out.write(f"{te_chr}\tTE_library\ttranscript\t1\t{length}\t.\t+\t.\t{tx_attr}")
 
             #  b) build exon fragments
             bps = sorted(te_breaks.get(te_chr, []))
             prev = 0
             if not bps:
                 # no sites → single exon
                 out.write(
                     f"{te_chr}\tTE_library\texon\t1\t{length}\t.\t+\t.\t{tx_attr}"
                 )
             else:
                 for bp in bps:
                     start = prev + 1
                     end   = bp
                     out.write(
                         f"{te_chr}\tTE_library\texon\t{start}\t{end}\t.\t+\t.\t{tx_attr}"
                     )
                     prev = bp
                 # final tail
                 if prev < length:
                     out.write(
                         f"{te_chr}\tTE_library\texon\t{prev+1}\t{length}\t.\t+\t.\t{tx_attr}"
                     )
            
    
    # build minimap index and splice bed for genome
    refgenome_dir  = os.path.dirname(os.path.realpath(ref_genome_fasta))
    refmm2_db_name    = ref_genome_fasta + ".mm2"
    refmm2_splice_bed = ref_gtf + ".mm2.splice.bed"
    which = "ref"
    prep_minimap2_reference(ref_genome_fasta, refgenome_dir, refmm2_db_name, refmm2_splice_bed, ref_gtf, which, ref_splice_bed=None)
    
    # build minimap index and splice bed for combined genome
    mm2_db_name    = combined_genomes_fa + ".mm2"
    mm2_splice_bed = combined_gtf + ".mm2.splice.bed"
    which = "comb"
    prep_minimap2_reference(combined_genomes_fa, TEIF_dir, mm2_db_name, mm2_splice_bed, combined_gtf, which, ref_splice_bed=refmm2_splice_bed)
        
    # build minimap index and splice bed for unmasked combined genome
    mm2_db_name    = combined_genomes_fa_unmasked + ".mm2"
    mm2_splice_bed = combined_genomes_fa_unmasked + ".mm2.splice.bed"
    which = "comb_un"
    prep_minimap2_reference(combined_genomes_fa_unmasked, TEIF_dir, mm2_db_name, mm2_splice_bed, combined_gtf, which, ref_splice_bed=refmm2_splice_bed)
    
    cln_cmd = f"rm {installed_TE_db} {filt_TE_db} {filt_TE_db}.clstr {TE_masked_genome} "
    pipeliner.add_commands([Command(cln_cmd, "cleanup.ok")])
    pipeliner.run()

    
    logger.info("-Done.")

    sys.exit(0)


if __name__=='__main__':
    main()

