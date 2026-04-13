#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import subprocess
import re
import shlex
from pathlib import Path
from argparse import RawTextHelpFormatter

# make sure your Pipeliner lives on PYTHONPATH
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "utils/PyLib"))
from Pipeliner import Pipeliner, Command

logging.basicConfig(
    format="%(asctime)-15s: %(levelname)s %(message)s",
    level=logging.INFO
)
logger = logging.getLogger()

def append_suffix(p: Path, extra: str) -> Path:
    # If there’s an existing suffix, keep it and append the extra.
    return p.with_suffix(p.suffix + extra) if p.suffix else p.with_name(p.name + extra)

def _load_breakinator_labels(breakinator_qc_path):
    """
    Return dict: read_id -> label. If a read has multiple rows, prioritize:
    Foldback > Chimeric > Pass.
    """
    priority = {"Foldback": 3, "Chimeric": 2, "Pass": 1}
    labels = {}

    if not os.path.exists(breakinator_qc_path):
        return labels  # empty -> treated as NA later

    with open(breakinator_qc_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            
            read_id = cols[6]
            label   = cols[-1]
            if label not in priority:
                # If Breakinator ever emits something unexpected, store but low priority
                priority.setdefault(label, 0)

            if (read_id not in labels) or (priority[label] > priority[labels[read_id]]):
                labels[read_id] = label

    return labels

def classify_reads(paf_file, outdir, fastq_basename, fastq_truebasename, fastq_dir):
    """
    Parse PAF, compute per-read fraction aligned to any TE using max alignment length,
    and write ID lists: pure_TE, chimeric, no_TE.
    Also include any high-MAPQ (>40) TE alignments in chimeric group.
    """
    max_len = {}
    qlen    = {}
    max_mapq = {}

    breakinator_qc = fastq_dir / f"{fastq_truebasename}_breakinator_qc.txt"
    brk_label = _load_breakinator_labels(breakinator_qc)

    with open(paf_file) as fh:
        for line in fh:
            cols = line.split()
            if len(cols) < 12:
                continue
            qname   = cols[0]
            length  = int(cols[1])
            aln_len = int(cols[10])
            mapq    = int(cols[11])

            max_len[qname]  = max(aln_len, max_len.get(qname, 0))
            qlen[qname]     = length
            max_mapq[qname] = max(mapq, max_mapq.get(qname, 0))

    os.makedirs(outdir, exist_ok=True)
    f_pure = open(outdir / f"{fastq_basename}_pure_TE.ids",  "w")
    f_chim = open(outdir / f"{fastq_basename}_chimeric.ids", "w")
    f_none = open(outdir / f"{fastq_basename}_no_TE.ids",    "w")

    # 4-column table: read_id, frac, maxMAPQ, Breakinator_Type
    read_frac = open(outdir / f"{fastq_basename}_read_TEfraction.txt", "w")

    for r in qlen:
        frac = max_len.get(r, 0) / qlen[r]
        mapq = max_mapq.get(r, 0)
        btyp = brk_label.get(r, "NA") if brk_label is not None else "NA"

        # always record the summary row (including Foldbacks)
        read_frac.write(f"{r}\t{frac:.4f}\t{mapq}\t{btyp}\n")

        # If Breakinator called this a Foldback, stop it from advancing
        if btyp == "Foldback":
            continue

        # otherwise, proceed with original classification logic
        # Reads are used if they either have >2.5% TE OR very good mapping score. 
        # For reads 1500 bp, this means >37bps
        if frac >= 0.95:
            f_pure.write(r + "\n")
        if frac > 0 or mapq >= 50:
            f_chim.write(r + "\n")
        else:
            f_none.write(r + "\n")

    f_pure.close(); f_chim.close(); f_none.close(); read_frac.close()
    logger.info(
        f"Wrote ID lists: {fastq_basename}_pure_TE.ids, {fastq_basename}_chimeric.ids, "
        f"{fastq_basename}_no_TE.ids and per-read table with Breakinator labels."
    )

def prep_minimap2_reference(genome_fa, inputdir, mm2_db_name, mm2_splice_file, ref_gtf, fastq_basename):

    chk = Path(f"{inputdir}/__mm2_prep_chkpts/{fastq_basename}.build.ok")
    if chk.exists():
        logger.info("Skipping prep ref: found %s", chk)
        return

    pip = Pipeliner(checkpoint_dir= inputdir / "__mm2_prep_chkpts")
    ctat_mm2_dir =  Path(__file__).resolve().parent / "utils/ctat-minimap2"
    
    # index genome
    cmd = f"{ctat_mm2_dir}/ctat-minimap2 -d {mm2_db_name} {genome_fa}"
    pip.add_commands([ Command(cmd, f"{fastq_basename}_mm2_prep_fusion_genome.ok") ])
    
    # Making the splice files
    cmd = (
        f"{ctat_mm2_dir}/misc/paftools.ctat.js gff2bed "
        f"{ref_gtf} > {mm2_splice_file}"
        )
    pip.add_commands([ Command(cmd, f"{fastq_basename}_make_fusion_splice.ok") ])

    pip.run()
    chk.touch()


def get_total_read_count(total_read_count_filename: str) -> int:
    """
    Read a file containing exactly one integer (the total read count)
    and return it.  Raise an error if the file is missing, unreadable,
    or doesn’t contain a valid integer.
    """
    if not os.path.exists(total_read_count_filename):
        raise FileNotFoundError(
            f"Error, cannot locate file: {total_read_count_filename}"
        )
    try:
        with open(total_read_count_filename, "r") as fh:
            line = fh.read().strip()
    except Exception as e:
        raise IOError(f"Error, cannot open file: {total_read_count_filename}: {e}")
    if not line.isdigit():
        raise ValueError(
            f"Error, not interpreting total read count [{line}] as integer value"
        )
    return int(line)

    
def main():
    parser = argparse.ArgumentParser(
        description="Classify ONT reads by how much they align to TEs"
    )
    #General Inputs
    parser.add_argument("--reads_fastq", dest="reads_fastq",   required=True,
                        help="your ONT cDNA reads (FASTQ or FASTA). The name of this file will serve as the basis of all output files")
                        
    parser.add_argument("--genome_fasta", dest="genome_fasta", required=True,
                   help="Standard genome fasta file (Ex: gencode.v44.genome.fa)."
                   "Should be in TEgenome_prep and after TE-lib-integration was run")
                   
    parser.add_argument("--genome_gtf", dest="genome_gtf", required=True,
                   help="Standard genome annotation file (Ex: gencode.v44.annotation.gtf)")
   
    #Workflow and efficiency
    parser.add_argument("--outdir", dest="outdir",        default=None,
                        help="where to put results (default: in parent directory of scripts)")

    parser.add_argument("--CPU", dest="CPU",           type=int, default=4,
                        help="threads for minimap2. Default = 4")
    
    parser.add_argument("--pychopper",  action="store_true",
                        help="Will trigger pychopper to run on samples. Default filter set to mean base quality score of 15")
    
    parser.add_argument("--kit_number",   type=str, default = "PCS114",
                        help="Identity of ONT kit passed into pychopper. Used to identify primer sequences")
    
    parser.add_argument("--prep_reference_only",   action="store_true")

    parser.add_argument("--examine_coding_effect", action="store_true")

    parser.add_argument("--no_ctat_mm2",  action="store_true",
                   help="Use to utilize regular minimap (might cause errors)")
                   
    parser.add_argument("--chim_candidates_only",  action="store_true",
                   help="Stop pipeline after the candidate identification")
                   
    parser.add_argument("--no_annot_filter",       action="store_true",
                   help="Skip the annotation filter. "
                   "Will have minimal effect but could introduce HLA and other genes that are more likely to be false pos")

    #Filtering Settigns                                            
    parser.add_argument("--min_per_id",            type=int,   default=70,
                   help="Minimum percent match required to consider minimap gene alignments")
                   
    parser.add_argument("--min_per_id_TE",            type=int,   default=60,
                   help="Minimum percent match required to consider minimap TE alignments")
                                                                                                    
    parser.add_argument("--min_FFPM",              type=float, default=0.75,
                   help="Minimum FFPM required for fusion transcripts to be considered. Default = 0.75")
    
    parser.add_argument("--max_exon_delta", type=int, default=10,
                   help="max exon delta for retrieval")

    parser.add_argument("--MIN_FRACTION_DOMINANT_ISO", type=float, default=0.05,
                   help="minimum expression of dominant isoform. Default = 0.05")

    parser.add_argument("--max_intron_length", type=int, default=100000,
                   help="maximum intron length during minimap2 search. Default = 100000")

    parser.add_argument("--shrink_intron_max_length", type=int, default=1000,
                   help="length for shrinking long introns to during the FusionInspector-style alignments. Default = 1000")

    parser.add_argument("--no_shrink_introns", action="store_true",
                   help="disable intron shrinking during FusionInspector-style alignments.")

    parser.add_argument("--snap_dist", type=int, default=3,
                   help="if breakpoint is at most this distance from a reference exon boundary, position gets snapped to the splice site. Default = 3")

    parser.add_argument("--min_trans_overlap", type=int, default=100,
                   help="minimum read overlap length for each gene in the fusion pair. Default = 100")

    parser.add_argument("--top_candidates_only", default = None, type=int, metavar="N",
                   help="only build FI contigs for the top N fusion candidates")

    parser.add_argument("--MIN_NUM_LR",            type=int,   default=1,
                   help="min number of long reads supporting fusion including canonical splicing at breakpoint. Default = 1")

    parser.add_argument("--min_LR_novel_junction_support",          type=int,   default=2,
                   help="min number of long reads with support at noncanonical splice breakpoint (eg. maybe not spliced!). Default = 2")

    parser.add_argument("--min_J",                 type=int,   default=1,
                   help="SR filter: minimum number of junction frags. Default = 1")

    parser.add_argument("--min_sumJS",             type=int,   default=1,
                   help="SR filter: minimum sum (junction + spanning) frags. Default = 1")

    parser.add_argument("--min_novel_junction_support",    type=int,   default=1,
                   help="SR filter: minimum number of junction reads required for novel (non-reference) exon-exon junction support")
                       
    parser.add_argument("--left_fq",   default=None, help="Illumina /1")
    
    parser.add_argument("--right_fq",  default=None, help="Illumina /2")

    parser.add_argument("--vis",                   action="store_true")

    parser.add_argument("--extract_fusion_LR_fasta", help="path to write evidence reads")

    parser.add_argument("--version",               action="store_true")

    args = parser.parse_args()
    
    #Defining the various libraries that will be used
    script_dir 		= Path(__file__).resolve().parent
    outdir 			= Path(args.outdir).expanduser().resolve() if args.outdir else script_dir.parent
    genome_lib 		= script_dir / "HG38_Genome_Data"
    TE_db_dir		= script_dir / "TEIF"
    init_algn_dir	= outdir / "Phase1_init_alignments"
    LR_fus_cand_dir = outdir / "Phase2_LR_fusion_cand"
    intermediates	= outdir / "Phase3_TELR_Fusion"
    finaloutputdir	= outdir / "Final_output"

        
    reads_fastq		= Path(args.reads_fastq).expanduser().resolve()
    fastq_truebasename = re.sub(r'\.(fastq|fq|fasta|fa)(\.gz)?$', '', reads_fastq.name, flags=re.I)
    fastq_basename  = re.sub(r'\.(fastq|fq|fasta|fa)(\.gz)?$', '_unmasked', reads_fastq.name, flags=re.I)
    fastq_dir 		= reads_fastq.parent
    
    outdir.mkdir(parents=True, exist_ok=True)
    init_algn_dir.mkdir(parents=True, exist_ok=True)
    LR_fus_cand_dir.mkdir(parents=True, exist_ok=True)
    intermediates.mkdir(parents=True, exist_ok=True)
    finaloutputdir.mkdir(parents=True, exist_ok=True)
        
    #─────── setting paths to genomes ───────────────────────────────────────────
    ref_fa		= Path(args.genome_fasta).expanduser().resolve() 
    ref_gtf		= Path(args.genome_gtf).expanduser().resolve() 
    
    if not os.path.isfile(ref_fa):
        sys.exit(f"Error: cannot find FASTA at {ref_fa}")
    if not os.path.isfile(ref_gtf):
        sys.exit(f"Error: cannot find GTF at {ref_gtf}")
    
    TEref_gtf       	= TE_db_dir / "ref_genome_plus_TE.gtf"
    TEref_fa        	= TE_db_dir / "ref_genome_plus_TE.fa"
    TEref_fa_unmasked	= TE_db_dir / "ref_genome_unmasked_plus_TE.fa"
    if not os.path.isfile(TEref_gtf):
        sys.exit(f"Error: cannot find GTF at {TEref_gtf}")
        
    
    #─────── Define variables ───────────────────────────────────────────────────
    num_threads     		= args.CPU
    snap_dist				= args.snap_dist
    min_trans_overlap		= args.min_trans_overlap
    MIN_FRACTION_DOMINANT_ISO = args.MIN_FRACTION_DOMINANT_ISO
    
        
    mm2_prog 			= "minimap2" if args.no_ctat_mm2 else script_dir / "utils/ctat-minimap2/ctat-minimap2"

    normfusionfastq		= init_algn_dir / f"{fastq_basename}_no_TE_reads.fastq"
    TEfusionfastq		= init_algn_dir / f"{fastq_basename}_chimeric_reads.fastq"
    
    ref_mm2_index		= append_suffix(ref_fa, ".mm2")
    ref_mm2_splice_bed	= append_suffix(ref_gtf, ".mm2.splice.bed")
    
    TE_mm2_index		= append_suffix(TEref_fa, ".mm2")
    TE_mm2_splice_bed	= append_suffix(TEref_gtf, ".mm2.splice.bed")
    
    unmaskTE_mm2_index		= append_suffix(TEref_fa_unmasked, ".mm2")
    unmaskTE_mm2_splice_bed	= append_suffix(TEref_fa_unmasked, ".mm2.splice.bed")

    te_fasta 			= TE_db_dir / "TE_db_filt_renamed.fa"
        

    #─────── check indexes exist ───────────────────────────────────────────────────
    if not os.path.exists(ref_mm2_index):
        sys.exit(f"Error, missing minimap2 index: {ref_mm2_index}. This should be made by TE-lib-integration")
    if not os.path.exists(ref_mm2_splice_bed):
        sys.exit(f"Error, missing splice BED: {ref_mm2_splice_bed}. This should be made by TE-lib-integration")
    if not os.path.exists(TE_mm2_index):
        sys.exit(f"Error, missing TE minimap2 index: {TE_mm2_index}. This should be made by TE-lib-integration")
    if not os.path.exists(TE_mm2_splice_bed):
        sys.exit(f"Error, missing TE splice BED: {TE_mm2_splice_bed}. This should be made by TE-lib-integration")
    
    
    #─────── Starting the first step. Identification of TE Containing Reads ─────────
    
    os.chdir(outdir)

    # Make the checkpoint directory 
    print("Writing into:", init_algn_dir)
    pip = Pipeliner(outdir / "__checkpt")

    # Check for TE fasta file ───────────────────────────────
    if not os.path.exists(te_fasta):
        sys.exit(f"Error, missing collapsed TE fasta: {te_fasta}. This should be made by TE-lib-integration")
        
    #─────── trim and filter raw FASTQ  ────────────────────────────────────
    if args.pychopper:
    	pychopper_fastq = init_algn_dir / f"{fastq_basename}_py.fastq"
    	cmd = (
    	f"pychopper -k args.kit_number -Q 15 -t {num_threads} {reads_fastq} {pychopper_fastq}"
    	)
    	reads_fastq = pychopper_fastq

    #─────── align reads against TE DB  ────────────────────────────────────
    paf = init_algn_dir / f"{fastq_basename}_TE_align.paf"
    cmd = (
        f"{mm2_prog} "
        f"-x map-ont -B2 --secondary=no -N 1 "
        f"-t {num_threads} --paf-no-hit "
        f"{te_fasta} {reads_fastq} "
        f"> {paf}"
    )
    pip.add_commands([ Command(cmd, f"{fastq_basename}_mm2_align.ok") ])
    pip.run()

    #─────── classify by fraction aligned ───────────────────────────────
    classify_reads(paf, init_algn_dir, fastq_basename, fastq_truebasename, fastq_dir)
    
    touch_file = init_algn_dir / f"{fastq_basename}_reads_classified.ok"
    open(touch_file, "a").close()

    #─────── extract the three groups of reads ──────────────────────────
    for label in ["pure_TE", "chimeric", "no_TE"]:
        ids = init_algn_dir / f"{fastq_basename}_{label}.ids"
        outf = init_algn_dir / f"{fastq_basename}_{label}_reads.fastq"
        cmd = f"seqtk subseq {reads_fastq} {ids} > {outf}"
        pip.add_commands([ Command(cmd, f"{fastq_basename}_extract_{label}.ok") ])

    pip.run()
    
    logger.info("Phase 1 complete, next stage starting")
    
    #─────── Starting Phase 2: Identifying the Fusion Transcripts ────────────────────────   

    #─────── run ctat-minimap2 to find normal chimeric candidates ────────────────────────   
    
    mm2_pref = LR_fus_cand_dir / fastq_basename
    prelim_bam = append_suffix(mm2_pref, ".mm2.prelim.bam")
    final_bam  = append_suffix(mm2_pref, ".mm2.bam")

    cmd_p2 = (
		f"{mm2_prog} --only_chimeric "
		f"--sam-hit-only --junc-bed {ref_mm2_splice_bed} "
		f"-ax splice -ub -t {num_threads} {ref_mm2_index} {normfusionfastq} "
		f"| samtools view -Sb -o {prelim_bam}"
    )
    pip.add_commands([ Command(cmd_p2, f"{fastq_basename}_run_mm2.ok") ])

    cmd_p2 = (
        f"samtools view -@ {num_threads} -h -d SA {prelim_bam} "
        f"| samtools sort -@ {num_threads} -N -o {final_bam}"
    )
    pip.add_commands([ Command(cmd_p2, f"{fastq_basename}_extract_chim.ok") ])
    pip.run()
    
    #─────── run ctat-minimap2 to find TE related chimeric candidates ────────────────────
    
    TEmm2_pref = LR_fus_cand_dir / fastq_basename 
    TEprelim_bam = append_suffix(TEmm2_pref, "_TEFusion.mm2.prelim.bam")
    TEfinal_bam  = append_suffix(TEmm2_pref, "_TEFusion.mm2.bam")

    cmd_p2 = (
		f"{mm2_prog} --only_chimeric "
		f"--sam-hit-only --junc-bed {TE_mm2_splice_bed} "
		f"-ax splice:hq -B2 -ub -t {num_threads} {TEref_fa} {TEfusionfastq} "
		f"| samtools view -Sb -o {TEprelim_bam}"
    )
    pip.add_commands([ Command(cmd_p2, f"{fastq_basename}_run_TEmm2.ok") ])

    cmd_p2 = (
        f"samtools view -@ {num_threads} -h -d SA {TEprelim_bam} "
        f"| samtools sort -@ {num_threads} -N -o {TEfinal_bam}"
    )
    pip.add_commands([ Command(cmd_p2, f"{fastq_basename}_extract_TEchim.ok") ])
    pip.run()
        

    #─────── merge sam files and convert to GFF3 and summarize ────────────────────────────
    
    merge_cmd = (
        f"samtools merge -f {LR_fus_cand_dir}/{fastq_basename}_all_fusions.bam "
        f"{final_bam} {TEfinal_bam} && "
        f"samtools sort -@ {num_threads} -N -o {LR_fus_cand_dir}/{fastq_basename}_all_fusions.sorted.bam {LR_fus_cand_dir}/{fastq_basename}_all_fusions.bam "
    )
    pip.add_commands([ Command(merge_cmd, f"{fastq_basename}_merge_all_fusions.ok") ])
    pip.run()

    gff3 = LR_fus_cand_dir / f"{fastq_basename}_all_fusions.sorted.bam.gff3"
    cmd_p2 = f"{script_dir}/utils/scripts/SAM_to_gxf.pl --sam {LR_fus_cand_dir}/{fastq_basename}_all_fusions.sorted.bam --format gff3 > {gff3}"
    pip.add_commands([ Command(cmd_p2, f"{fastq_basename}_to_gff3.ok") ])

    chims_out = append_suffix(mm2_pref, ".chims_described")
    debug = append_suffix(mm2_pref, ".debug.txt")
    cmd_p2 = (
        f"{script_dir}/utils/scripts/genome_gff3_to_chim_summary.pl "
        f"--align_gff3 {gff3} --annot_gtf {TEref_gtf} --min_per_id_TE {args.min_per_id_TE} "
        f"--min_per_id {args.min_per_id} --TE_only > {chims_out} "
    )
    pip.add_commands([ Command(cmd_p2, f"{fastq_basename}_summarize_chims.ok") ])
    pip.run()


    logger.info("Phase 2 complete, next stage starting…")

    #─────── Starting Phase 3: Filtering the Fusion Transcripts ────────────────────────   


    #─────── convert to FASTA (if needed) and count total reads -----
    if reads_fastq.suffix.lower() == ".gz":
        uncompressed = reads_fastq.with_suffix("")
        cmd_p3 = f"gunzip -c {reads_fastq} > {uncompressed}"
        pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_gunzip_transcripts.ok") ])
        reads_fastq = uncompressed
        
    if reads_fastq.suffix.lower() in {".fastq", ".fq"}:
        fasta_path = reads_fastq.with_suffix(".fasta")        
        cmd_p3 = f"seqtk seq -a {reads_fastq} > {fasta_path}"
        pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_fastq_to_fasta_seqtk.ok") ])
        reads_fastq = fasta_path
    pip.run()         
        
    total_reads_file = intermediates / f"{fastq_basename}.LR_read_count.txt"
    
    cmd_p3 = f"seqtk comp {reads_fastq} | wc -l > {total_reads_file}"
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_count_reads.ok")])
    pip.run()         
    
    with open(total_reads_file) as fh:
        num_total_reads = int(fh.read().strip())
    print(f"{fastq_basename} Total long‐reads = {num_total_reads}")
    
    #─────── extract candidate fusion transcripts -----
    cand_prefix 	= intermediates / f"{fastq_basename}_chimeric_read_candidates"
    chim_fa     	= f"{cand_prefix}.transcripts.fa"
    FI_list     	= f"{cand_prefix}.FI_listing"
    FI_list_wreads	= f"{FI_list}.with_reads"
    
    cmd_p3 = (
        f"{script_dir}/utils/scripts/retrieve_fusion_transcript_candidates.pl "
        f"--trans_fasta {reads_fastq} "
        f"--chims_described {chims_out} "
        f"--max_exon_delta {args.max_exon_delta} "
        f"--num_total_reads {num_total_reads} "
        f"--min_FFPM {args.min_FFPM} "
        f"--output_prefix {cand_prefix}"
    )
    if args.chim_candidates_only:
        cmd_p3 += " --skip_read_extraction"
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_chim_candidates.ok")])
    pip.run()

#      #───────rezip the fasta file to conserve space
#     if reads_fastq.suffix.lower() == ".gz":
#         logger.info("Already gzipped: %s", reads_fastq)
#     else:
#         gz_path = reads_fastq.with_suffix(reads_fastq.suffix + ".gz")  
#         if gz_path.exists():
#             logger.info("Compressed file already exists, skipping: %s", gz_path)
#         else:
#             cmd = f"pigz --best -p {num_threads} {shlex.quote(str(reads_fastq))}"
#             pip.add_commands([Command(cmd, f"{fastq_basename}_rezip.ok")])
    
    #─────── Annotation of candidates -----
    cmd_p3 = f"{script_dir}/utils/scripts/FusionAnnotator --genome_lib_dir {genome_lib} --annotate {FI_list} > {FI_list}.wAnnot"
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_chim_candidates_fasta.FI_listing.annotate.ok")])
    
    FI_list    = f"{FI_list}.wAnnot"
    
    if not args.no_annot_filter:
        cmd_p3 = (
            f"{script_dir}/utils/scripts/filter_by_annotation_rules.pl "
            f"--fusions {FI_list} "
            f"--genome_lib_dir {genome_lib} "
            f"--custom_exclusions NEIGHBORS_OVERLAP"
        )
        
        pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_filter_by_annot_rules.ok") ])
        FI_list = f"{FI_list}.annot_filter.pass"
    pip.run()
    

    #─────── finalize fasta of candidates -----
    revised_fa = f"{chim_fa}.revised.fasta"
    cmd_p3 = (
        f"{script_dir}/utils/scripts/revise_fusion_reads_fasta.pl "
        f"{FI_list} {FI_list_wreads} {chim_fa} > {revised_fa}"
    )
    
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_revise_chimeric_reads_fasta.ok")])
    pip.run()
    
    
    #─────── Phase 3b: FusionInspector‐style contig build & precise breakpoint resolution ───


    #─────── Stage 1: generate the mini‐genome contigs for FusionInspector
    FI_outprefix = intermediates / f"{fastq_basename}_LR-FI_targets"
    cmd_p3 = (
        f"{script_dir}/utils/scripts/fusion_pair_to_mini_genome_join.pl "
        f"--fusions       {FI_list} "
        f"--gtf           {TEref_gtf} "
        f"--genome_fa     {TEref_fa_unmasked} "
        f"--out_prefix    {FI_outprefix}"
    )
    if getattr(args, "no_shrink_introns", False) is False and hasattr(args, "shrink_intron_max_length"):
        cmd_p3 += f" --shrink_introns --max_intron_length {args.shrink_intron_max_length}"
    if args.top_candidates_only is not None:
        cmd_p3 += f" --top_candidates_only {args.top_candidates_only}"

    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_FI_contigs.ok") ])
    pip.run()

    #─────── prepare minimap2 index + splice‐BED for those contigs
    FI_splice_bed   = intermediates / f"{fastq_basename}_LR-FI_targets.gtf.mm2.splice.bed"
    FI_mm2_db       = intermediates / f"{fastq_basename}_LR-FI_targets.fa.mm2"
    
    FI_contigs_fa   = intermediates / f"{fastq_basename}_LR-FI_targets.fa"
    FI_annots_gtf  = intermediates / f"{fastq_basename}_LR-FI_targets.gtf"

    prep_minimap2_reference(
        FI_contigs_fa,
        intermediates,
        FI_mm2_db,
        FI_splice_bed,
        FI_annots_gtf,
        fastq_basename
    )

    #─────── Stage 3: align your chimeric FASTA back to the new contigs
    LR_FI_mm2_bam 	= intermediates / f"{fastq_basename}_LR-FI.mm2.bam"
    cmd_p3 = (
        f"bash -c \"set -eou pipefail && {mm2_prog} --sam-hit-only "
        f"-ax splice -u b --junc-bed {FI_splice_bed} -t {num_threads} "
        f"{FI_mm2_db} {revised_fa} "
        f"| samtools view -Sb -o {LR_FI_mm2_bam}\" "
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_LR-FI.mm2.ok") ])
    pip.run()

    #─────── Stage 4: convert that BAM into GFF3 (include non‐primary alignments)
    LR_FI_gff3 = intermediates / f"{fastq_basename}_LR-FI.mm2.gff3"
    cmd_p3 = (
        f"{script_dir}/utils/scripts/SAM_to_gxf.pl "
        f"--sam {LR_FI_mm2_bam} "
        f"--format gff3 --allow_non_primary "
        f"> {LR_FI_gff3}"
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_LR-FI.mm2.sam_to_gff3.ok") ])
    pip.run()

    #─────── call the “seq‐similar” helper to annotate repetitive regions
    LR_FI_similar_seq_gff3 = intermediates / f"{fastq_basename}_LR-FI_targets.seqsimilar_regions.gff3"
    cmd_p3 = (
        f"{script_dir}/utils/scripts/get_seq_similar_region_FI_coordinates.pl "
        f"--finspector_gtf {FI_annots_gtf} "
        f"--genome_lib_dir {genome_lib} "
        f"> {LR_FI_similar_seq_gff3}"
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_FI_targets_seqsim_gff3.ok") ])
    pip.run()

    #─────── extract the fusion alignments into a breakpoint file
    LR_fusion_output_prefix = intermediates / f"{fastq_basename}_LR-FI.mm2.fusion_transcripts"
    cmd_p3 = (
        f"{script_dir}/utils/scripts/LR-FI_fusion_align_extractor.pl "
        f"--FI_gtf              {FI_annots_gtf} "
        f"--LR_gff3             {LR_FI_gff3} "
        f"--seq_similar_gff3    {LR_FI_similar_seq_gff3} "
        f"--output_prefix       {LR_fusion_output_prefix} "
        f"--snap_dist           {snap_dist} "
        f"--min_trans_overlap_length {min_trans_overlap} "
        f"> {LR_fusion_output_prefix}"
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_LR-FI.mm2.sam_to_gff3.extract_fusions.ok") ])
    pip.run()

    #─────── incorporate long‐read FFPM back into the breakpoint table
    LR_total_count = get_total_read_count(total_reads_file)
    fusions_filename =  append_suffix(LR_fusion_output_prefix, ".breakpoint_info.tsv")
    cmd_p3 = (
        f"{script_dir}/utils/scripts/incorporate_LR_FFPM.pl "
        f"--fusions      {fusions_filename} "
        f"--num_LR_total {LR_total_count} "
        f"--output_file {fusions_filename}.w_LR_FFPM"
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_added_LR_FFPM.ok") ])
    pip.run()
    
    #─────── update variable for subsequent steps
    fusions_filename = f"{fusions_filename}.w_LR_FFPM"
    
    ##########################
    ### Fusion Annotator   ###
    ##########################
    
    cmd_p3 = f"{script_dir}/utils/scripts/FusionAnnotator --genome_lib_dir {genome_lib} --annotate {fusions_filename} > {fusions_filename}.wAnnot"
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_annotate_fusions.ok")])
    pip.run()
    
    fusions_filename    = f"{fusions_filename}.wAnnot"

    #─────── Copy preliminary outputs to base directory
    prelim_report = finaloutputdir / f"{fastq_basename}_ctat-LR-fusion.fusion_predictions.preliminary.tsv"
    cmd_p3 = f"cp {fusions_filename} {prelim_report}"
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_copy_to_prelim.ok")])
    pip.run()
    
    #─────── Add abridged version without the evidence read names
    cmd_p3 = (
        f"{script_dir}/utils/scripts/column_exclusions.pl {prelim_report} "
        f" LR_accessions,JunctionReads,SpanningFrags,CounterFusionLeftReads,CounterFusionRightReads "
        f" > {prelim_report}.abridged.tsv "
    )
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_abridged_prelim_preds.ok")])
    pip.run()
    
    ##########################
    #### Fusion Filtering ####
    ##########################
    
    cmd_p3 = (
        f"{script_dir}/utils/scripts/blast_and_promiscuity_filter.pl "
        f" --fusion_preds {fusions_filename} --out_prefix {fusions_filename} --genome_lib_dir {genome_lib} --exclude_loci_overlap_check "
    )
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_blast_promisc_filter.ok")])
    pip.run()
    
    fusions_filename = f"{fusions_filename}.post_blast_and_promiscuity_filter"
    
    cmd_p3 = (
        f"{script_dir}/utils/scripts/filter_LR_fusions_by_evidence_abundance.py "
        f" --min_num_LR {args.MIN_NUM_LR} "
        f" --min_FFPM {args.min_FFPM} "
        f" --min_LR_novel_junction_support {args.min_LR_novel_junction_support} "
        f" --min_J {args.min_J} "
        f" --min_sumJS {args.min_sumJS} "
        f" --min_novel_junction_support {args.min_novel_junction_support} "
        f" --fusions_input {fusions_filename} "
        f" --filtered_fusions_output {fusions_filename}.filt_by_min_reads "
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_filter_LR_fusions_by_abundance.ok") ])
    
    fusions_filename = f"{fusions_filename}.filt_by_min_reads"
    
    if MIN_FRACTION_DOMINANT_ISO > 0:
        filt_out = f"{fusions_filename}.filt_by_min_dom_iso_frac"
        cmd_p3 = (
            f"{script_dir}/utils/scripts/filter_low_pct_dom_iso.py "
            f"--min_frac_dom_iso {MIN_FRACTION_DOMINANT_ISO} "
            f"--fusions_input {fusions_filename} "
            f"--filtered_fusions_output {filt_out}"
        )
        pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_filter_by_min_dom_iso_frac.ok") ])
        fusions_filename = filt_out

    if args.examine_coding_effect:
        coding_out = f"{fusions_filename}.w_coding_effect"
        cmd_p3 = (
            f"{script_dir}/utils/scripts/fusion_to_coding_region_effect.pl "
            f"--fusions {fusions_filename} "
            f"--genome_lib_dir {genome_lib} "
            f"> {coding_out}"
        )
        pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_coding_eff.ok") ])
        fusions_filename = coding_out


    #─────── Copy final predictions into the deliverable
    final_preds = finaloutputdir / f"{fastq_basename}_ctat-LR-fusion.fusion_predictions.tsv"
    cmd_p3 = f"cp {fusions_filename} {final_preds}"
    
    pip.add_commands([Command(cmd_p3, f"{fastq_basename}_copy_final_predictions_to_deliverable.{int(args.examine_coding_effect)}.ok")])

    #─────── Make an abridged version (dropping heavy evidence columns)
    abridged = finaloutputdir / f"{fastq_basename}_ctat-LR-fusion.fusion_predictions.abridged.tsv"
    cmd_p3 = (
        f"{script_dir}/utils/scripts/column_exclusions.pl {final_preds} "
        f"LR_accessions,JunctionReads,SpanningFrags,CounterFusionLeftReads,CounterFusionRightReads "
        f"> {abridged}"
    )
    pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_abridged_final_preds.{int(args.examine_coding_effect)}.ok") ])

    #─────── (optionally) extract the fusion‐supporting FASTA reads
    fusion_fasta_reads= finaloutputdir / f"{fastq_basename}_fusion_reads.fasta"
    if args.extract_fusion_LR_fasta:
        cmd_p3 = (
            f"{script_dir}/utils/scripts/extract_fusion_evidence_reads.pl "
            f"--fusions {final_preds} "
            f"--reads_fasta {revised_fa} "
            f"--reads_output {fusion_fasta_reads}"
        )
        pip.add_commands([ Command(cmd_p3, f"{fastq_basename}_extract_fusion_reads.ok") ])

    #─────── finally run all of the above
    pip.run()
    
    logger.info("Phase3 complete. Outputs in %s", outdir)


if __name__ == "__main__":
    main()
