#!/usr/bin/env python3
import argparse, csv, os, sys, shutil, json, gzip
from typing import List, Optional, Dict, Tuple
from pathlib import Path
from collections import defaultdict, Counter
import subprocess
import re

try:
	csv.field_size_limit(sys.maxsize)
except OverflowError:
	csv.field_size_limit(10**9)

# ------------------ utilities ------------------

def run(cmd, **kw):
	print(f"[cmd] {cmd}", flush=True)
	subprocess.run(cmd, shell=True, check=True, **kw)

def ensure(p: Path) -> Path:
	p.mkdir(parents=True, exist_ok=True)
	return p

def slug(s: str) -> str:
	"""Safe file-name token: keep letters, digits: and replace others with _"""
	s2 = re.sub(r"[^A-Za-z0-9]+", "_", s)
	return re.sub(r"_+", "_", s2).strip("_")
	

def find_fastq(sample_base: str, fastq_dir: Path) -> Path | None:
	# Try common patterns to find fastq file
	for suf in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
		p = fastq_dir / f"{sample_base}{suf}"
		if p.exists():
			return p


def parse_selected_fusions(tsv: Path, min_num_lr: int):
	"""
	Yield dicts for rows meeting num_LR >= min_num_lr.
	Robust to header variants (e.g., '#FusionName' vs 'FusionName',
	"""
	rows = []
	with open(tsv, newline="") as fh:
		reader = csv.DictReader(fh, delimiter="\t")
		if reader.fieldnames is None:
			raise RuntimeError(f"{tsv} has no header?")
		hdr = reader.fieldnames

		# column keys
		key_name = next((k for k in ("#FusionName", "FusionName", "Fusion") if k in hdr), None)
		key_nlr	 = next((k for k in ("num_LR", "numLR") if k in hdr), None)
		key_lb	 = next((k for k in ("LeftBreakpoint","left_breakpoint") if k in hdr), None)
		key_rb	 = next((k for k in ("RightBreakpoint","right_breakpoint") if k in hdr), None)
		key_reads = next((k for k in hdr if k.lower().endswith("lr_accessions")), None)


		if not all([key_name, key_nlr, key_lb, key_rb, key_reads]):
			raise RuntimeError(
				f"{tsv} header not recognized.\n"
				f"Have: {hdr}\n"
				f"Need something like: #FusionName / num_LR / ... / LeftBreakpoint / ... / RightBreakpoint / ... / LR_accessions"
			)

		for r in reader:
			try:
				nlr = int(str(r.get(key_nlr, "0")).strip() or 0)
			except ValueError:
				nlr = 0
			if nlr < min_num_lr:
				continue
			name = str(r.get(key_name, "")).strip()
			lb	 = str(r.get(key_lb, "")).strip()
			rb	 = str(r.get(key_rb, "")).strip()
			accs = str(r.get(key_reads, "")).strip()
			if not name or not accs:
				continue
			# split comma-separated list; strip spaces
			read_ids = [a.strip() for a in accs.split(",") if a.strip()]
			rows.append({
				"fusion": name,
				"left_bp": lb,
				"right_bp": rb,
				"num_LR": nlr,
				"read_ids": read_ids
			})
	return rows

def find_group_dirs(base_dir: Path, cell: str, fusion_name: str) -> List[Path]:
	"""
	Inside base_dir, find subdirs like: {cell}__{fusion_slug}__LBP_*
	"""
	fusion_slug = slug(fusion_name)
	prefix = f"{cell}__{fusion_slug}__LBP_"
	return sorted([p for p in base_dir.iterdir() if p.is_dir() and p.name.startswith(prefix)])

def pick_reads_fastq(d: Path) -> Optional[Path]:
	"""
	Prefer filtered reads; else fall back to raw. If multiple, pick largest.
	"""
	cands = list(d.glob("*_filt_reads.fastq.gz"))
	if not cands:
		cands = list(d.glob("*.reads.fastq.gz"))
	if not cands:
		return None
	if len(cands) == 1:
		return cands[0]
	return max(cands, key=lambda p: p.stat().st_size)
	
# ------------------ pipeline steps ------------------

def extract_reads_seqtk(ids: list[str], fastq: Path, out_fastq_gz: Path, threads: int):
	if not ids:
		return
	ensure(out_fastq_gz.parent)
	ids_txt = out_fastq_gz.with_suffix(".ids.txt")
	with open(ids_txt, "w") as ofh:
		for rid in sorted(set(ids)):
			print(rid, file=ofh)
	
	run(f"seqkit grep -j {threads} -f {ids_txt} {fastq} > {out_fastq_gz}")
	ids_txt.unlink(missing_ok=True)

def cat_gz(inputs: list[Path], out_gz: Path):
	ensure(out_gz.parent)
	parts = " ".join(str(p) for p in inputs)
	run(f"cat {parts} > {out_gz}")

def rnabloom_assemble(reads_gz: Path, outdir: Path, threads: int, tag: str, extra: str|None=None) -> Path:
	"""
	Run RNABloom (long-read mode). Returns path to assembled transcripts (FASTA).
	"""
	ensure(outdir)
	outdirprefix = outdir 
	outname = f"{tag}_rnabloom"
	extra = extra or ""

	cmd = f"rnabloom -long {reads_gz} -o {outdirprefix} -n {outname} -t {threads} {extra}".strip()
	run(cmd)
	fa = outdirprefix / f"{outname}.transcripts.fa"
	return fa

def minimap2_align(target_fa: Path, reads_gz: Path, paf_out: Path, threads: int):
	run(f"minimap2 -t {threads} -x lr:hq {target_fa} {reads_gz} > {paf_out}")

def medaka_consensus(reads_gz: Path, draft_fa: Path, outdir: Path, threads: int, model: str, medaka_bin="medaka_consensus") -> Path:
	ensure(outdir)
	run(f"{medaka_bin} -i {reads_gz} -d {draft_fa} -o {outdir} -t {threads} -M 20 -m {model}")
	cons = outdir / "consensus.fasta"
	if not cons.exists():
		raise RuntimeError(f"Medaka did not create {cons}")
	return cons

# ------------------ main ------------------

def main():
	ap = argparse.ArgumentParser(
		description="Assemble and polish fusion transcripts per fusion breakpoint across replicates."
	)
	ap.add_argument("--base_dir", required=True,
					help="Base directory containing Phase3_TELR_Fusion, IGV_prep, Final_output, etc. "
						 "Output 'Fusion_Alignments' will be created here.")
	ap.add_argument("--fusion_list", required=True,			 
					help="Text file: one fusion name per line (same as for IGV script) (e.g. ENSG...--chrTE_XXX)")		 
	ap.add_argument("--fastq_dir", required=True,
					help="Directory with long-read FASTQ/FASTQ.GZ for each sample (same basename as *_selected_fusions.tsv).")
	ap.add_argument("--Sample_prefix", required=True,
					help="Prefix of samples to include (e.g., GB24). Matches *_selected_fusions.tsv and corresponding FASTQs.")
	ap.add_argument("--min_num_LR", type=int, default=10,
					help="Minimum LR reads required per fusion-breakpoint (default: 10).")
	ap.add_argument("--threads", type=int, default=8,
					help="Threads for align/polish tools (default: 8).")

	# tool options / misc
	ap.add_argument("--genome_fasta", default=None,
					help="Path to genome fasta file (default is human genome in TEgenome_prep)")
	ap.add_argument("--all_breakpoints", action="store_true",
					help="Run RNABloom and medaka on all breakpoint samples individually in addition to the combined reads.")
	ap.add_argument("--rnabloom_extra", default=None,
					help="Additional RNABloom flags (optional, passed as-is).")
	ap.add_argument("--medaka_model", default=None,
					help="Medaka model name (e.g., r1041_e82_400bps_sup). If omitted, medaka step is skipped."
					"For more information see: https://github.com/nanoporetech/medaka"
					"Can also run medaka tools resolve_model --auto_model consensus {target}.fastq")
	ap.add_argument("--medaka_bin", default="medaka_consensus",
					help="Medaka consensus executable (default: medaka_consensus).")
	ap.add_argument("--keep_intermediate", action="store_true",
					help="Keep per-sample extracted FASTQs and intermediate files.")

	args = ap.parse_args()

	base_dir	= Path(args.base_dir).expanduser().resolve()
	fastq_dir	= Path(args.fastq_dir).expanduser().resolve()
	script_dir	= Path(__file__).resolve().parent
	util_dir	= ensure(script_dir / "utils")
	igv_dir		= ensure(script_dir / "IGV_prep")
	TEIF_dir	= ensure(script_dir / "TEIF")
	tegenome_dir= ensure(script_dir / "TEgenome_prep")
	out_root	= ensure(base_dir / "Fusion_Alignments")
	
	fusion_list = Path(args.fusion_list).expanduser().resolve()
	genome_fasta= Path(args.genome_fasta).expanduser().resolve() if args.genome_fasta else tegenome_dir / "gencode.v44.genome.fa"
	prefix		= args.Sample_prefix
	threads		= args.threads
	min_num_lr	= args.min_num_LR
	adapter_ont = tegenome_dir / "ONT_adapter_fromSequali_all.fa"  

	print("Starting data processing...")

	tsvs = sorted(igv_dir.glob(f"{prefix}_*/{prefix}_*_selected_fusions.tsv"))
	if not tsvs:
		sys.exit(f"ERROR: no files like {igv_dir}/{prefix}*/{prefix}*_selected_fusions.tsv") #These files are made by the previous script

	# Collecting read id for fusions with read support over the minimum
	# key = (fusion_name, left_bp, right_bp)
	groups: dict[tuple[str,str,str], dict] = {}
	for tsv in tsvs:
		sample_base = tsv.name.replace("_unmasked_selected_fusions.tsv", "")
		
		#Searching for fastq to ensure it exists prior to proceeding
		fq = find_fastq(sample_base, fastq_dir)
		if not fq:
			sys.exit(f"ERROR: cannot find FASTQ for sample {sample_base} in {fastq_dir}")
			
		for row in parse_selected_fusions(tsv, min_num_lr):
			key = (row["fusion"], row["left_bp"], row["right_bp"])
			g = groups.setdefault(key, {"samples": defaultdict(set), "total_reads": 0, "sample_bases": set()})
			g["samples"][sample_base].update(row["read_ids"])
			g["total_reads"] += len(row["read_ids"])
			g["sample_bases"].add(sample_base)

	if not groups:
		sys.exit("ERROR: after filtering, no fusion-breakpoints meet min_num_LR")

	# Identifying if the adapter sequences are in the fusion contig itself. If so, then filtering below should be skipped
	fusion_contigs_prefix = out_root / f"{prefix}_fusioncontigs"
	run(f"perl {util_dir}/scripts/fusion_pair_to_mini_genome_join.pl --fusions {fusion_list} " 
		f"--gtf {TEIF_dir}/ref_genome_unmasked_plus_TE.gtf --genome_fa {TEIF_dir}/ref_genome_unmasked_plus_TE.fa "
		f"--shrink_introns --max_intron_length 1000 --out_prefix {fusion_contigs_prefix} ")
		
	# Identifies fusions that have an ONT adapter internally. These "hits" are more likely to arise from PCR artifacts 
	filtfusion_contigs = out_root / f"{prefix}_fusioncontigs_filtered.fa"
	run(f"cutadapt {out_root}/{prefix}_fusioncontigs.fa -g file:{adapter_ont} --discard "
		f"--rc --overlap 99 --error-rate 0.2 --max-aer 0.15 -j 0 -o {filtfusion_contigs} ")
	
	fusion_contigs	= out_root / f"{prefix}_fusioncontigs.fa"
	fusion_contigsgtf = out_root / f"{prefix}_fusioncontigs.gtf"
	orig_ids		= out_root / f"{prefix}_orig.ids"
	kept_ids		= out_root / f"{prefix}_kept.ids"
	removed_ids_path = out_root / f"{prefix}_contig_with_adapters.ids"
	comb_genome		= out_root / f"{prefix}_comb_genome.fa"
	comb_genome_gtf	= out_root / f"{prefix}_comb_genome.gtf"
	
	run(f"seqkit fx2tab -j {threads} -n {fusion_contigs} | awk '{{print $1}}' | sort -u	> {orig_ids} ")
	run(f"seqkit fx2tab -n -j {threads} {filtfusion_contigs} | awk '{{print $1}}' | sort -u > {kept_ids} ")
	run(f"comm -23 {orig_ids} {kept_ids} > {removed_ids_path} ")
	
	run(f"cat {genome_fasta} {fusion_contigs} > {comb_genome} ")
	run(f"cat {tegenome_dir}/gencode.v44.annotation.gtf {fusion_contigsgtf} > {comb_genome_gtf} ")
	run(f"samtools faidx {comb_genome} ")
	
	removed_ids = set()
	if removed_ids_path.exists():
		with open(removed_ids_path) as fh:
			for line in fh:
				rid = line.strip()
				if rid:
					removed_ids.add(rid)
					
	fusion_to_group_dirs: Dict[str, List[Path]] = defaultdict(list)

	for (fusion, lbp, rbp), info in groups.items():
		fusion_slug = slug(fusion)
		bp_slug		= f"LBP_{slug(lbp)}_RBP_{slug(rbp)}"
		tag			= f"{prefix}__{fusion_slug}__{bp_slug}"
		fusion_dir	= ensure(out_root / f"{prefix}__{fusion_slug}")
		group_dir	= ensure(fusion_dir / tag)
		
		out_bam = fusion_dir / f"{prefix}__{fusion_slug}_sorted.bam" 
		
		if out_bam.exists(): 
			print(f"[resume] Found existing final BAM+index for {fusion} → {out_bam}. Skipping assembly/polish/alignment.") 
			continue
		
		fusion_to_group_dirs[fusion].append(group_dir)

		# Will extract all of the fastq reads that match for the specific fusion breakpoint 
		per_sample_fastqs = []
		for sample_base, ids in sorted(info["samples"].items()):
			fq = find_fastq(sample_base, fastq_dir)
			if fq is None:
				sys.exit(f"ERROR: missing FASTQ for {sample_base}")
				
			out_fq_gz = group_dir / f"{sample_base}.reads.fastq.gz"
			extract_reads_seqtk(list(ids), fq, out_fq_gz, threads)
			per_sample_fastqs.append(out_fq_gz)

		combined_gz = group_dir / f"{tag}.reads.fastq.gz"
		cat_gz(per_sample_fastqs, combined_gz)
		
		
		if not args.keep_intermediate:
			for p in per_sample_fastqs:
				p.unlink(missing_ok=True)
		
		filtered_comb = combined_gz
		
		if fusion not in removed_ids:
			#Scan reads for presence of adapters in the middle
			filtered_comb = group_dir / f"{tag}_filt_reads.fastq.gz"
			info_file = group_dir / f"{tag}_cutadapt_hits.tsv"
		
			run(f"cutadapt -g file:{adapter_ont} --discard --info-file {info_file} --overlap 99 --error-rate 0.2 --max-aer 0.15 -j 0 -o {filtered_comb} {combined_gz} ")
		
			kept = int(subprocess.check_output(
				f"seqkit stats -T {filtered_comb} | awk 'NR==2{{print $4}}'",
				shell=True, text=True).strip() or 0)
			
			orig = int(subprocess.check_output(
				f"seqkit stats -T {combined_gz} | awk 'NR==2{{print $4}}'",
				shell=True, text=True).strip() or 0)

			removed = max(0, orig - kept)
			pct = (100.0 * removed / orig) if orig else 0.0

			# write a summary
			with open(group_dir / f"{tag}_adapter_filter_summary.txt", "a") as ofh:
				ofh.write(f"Removed {removed}/{orig} ({pct:.1f}%)  kept={kept}	total={orig}\n")

			# stop if >50% removed
			if orig and removed > orig / 2:
				(group_dir / "RT_ARTIFACT_LIKELY").write_text(
					f">50% adapter/RT artifacts: removed {removed}/{orig} ({pct:.1f}%).\n"
					f"Input:	{combined_gz.name}\nFiltered: {filtered_comb.name}\n"
				)
				continue
		
		if args.all_breakpoints:
			# RNABloom assembly
			print("Starting RNABloom Assembly")

			rnab_dir = ensure(group_dir / "rnabloom")
			draft = rnabloom_assemble(filtered_comb, rnab_dir, threads, tag, extra=args.rnabloom_extra)

			# medaka (optional)
			print("Starting medaka polishing")

			medaka_path = None
			if args.medaka_model:
				medaka_dir = ensure(group_dir / "medaka")
				medaka_path = medaka_consensus(filtered_comb, draft, medaka_dir, threads,
										   model=args.medaka_model, medaka_bin=args.medaka_bin)
	
	
	for fusion, group_dirs in fusion_to_group_dirs.items():
		fusion_slug = slug(fusion)
		fusion_dir	= ensure(out_root / f"{prefix}__{fusion_slug}")
		tag			= f"{prefix}__{fusion_slug}"
		
		out_bam = fusion_dir / f"{tag}_sorted.bam"
		if out_bam.exists():
			print(f"[resume] Found existing final BAM+index for {fusion} → {out_bam}. Skipping assembly/polish/alignment.")
			continue		
		
		#Find all subdirectories with reads and combine them
		group_dirs = find_group_dirs(fusion_dir, prefix, fusion_slug)
		fq_list: List[Path] = []
		for gd in group_dirs:
			fq = pick_reads_fastq(gd)
			if fq:
				fq_list.append(fq)
		
		if not fq_list:
			print(f"[warn] No reads found for fusion '{fusion}'. Skipping.", file=sys.stderr)
			continue

		
		allcombined_gz = fusion_dir / f"{prefix}__{fusion_slug}.allreads.fastq.gz"
		dedupcombined_gz = fusion_dir / f"{prefix}__{fusion_slug}.dedupreads.fastq.gz"

		cat_gz(fq_list, allcombined_gz)
	
		#Deduplicate the reads to avoid any biasing of consensus
		run(f"seqkit rmdup -s -i -j {threads} -o {dedupcombined_gz} {allcombined_gz}")
		
		if not args.keep_intermediate:
			try: allcombined_gz.unlink()
			except: pass

	
		# RNABloom assembly
		print("Starting full RNABloom Assembly")
		rnab_dir = ensure(fusion_dir / "rnabloom_All")
		draft = rnabloom_assemble(dedupcombined_gz, rnab_dir, threads, tag, extra=args.rnabloom_extra)
	
		# medaka (optional)
		medaka_path = None
		if args.medaka_model:
			print("Starting medaka polishing")
			medaka_dir = ensure(fusion_dir / "medaka_All")
			medaka_path = medaka_consensus(dedupcombined_gz, draft, medaka_dir, threads,
										   model=args.medaka_model, medaka_bin=args.medaka_bin)
		
		out_bam = fusion_dir / f"{tag}_sorted.bam"
		# alignment. Since RNABloom and Medaka put out transcript structures, we use the splicing alignment 
		if args.medaka_model:
			run(f"minimap2 -ax splice:hq -N 1 -t {threads} {comb_genome} {medaka_path} | samtools sort -O BAM -o {out_bam}")
			run(f"samtools index {out_bam}")
		
		#Aligns reads to the fusion contigs only. 
		out_bam_reads = fusion_dir / f"{tag}_sorted_contig_reads.bam"
		run(f"minimap2 -ax splice:hq --secondary=no -t {threads} {fusion_contigs} {dedupcombined_gz} | samtools sort -O BAM -o {out_bam_reads}")

		#Aligns reads to the full combined genome. This will result in fusions that arise from internal TEs to map to those areas of the genome. 
		out_bam_genomereads = fusion_dir / f"{tag}_sorted_genome_reads.bam"
		run(f"minimap2 -ax splice:hq --secondary=no -t {threads} {comb_genome} {dedupcombined_gz} | samtools sort -O BAM -o {out_bam_genomereads}")
	
	# Clean up and remove unnecessary files		   
	run(f"rm {orig_ids} {kept_ids} {removed_ids_path} {filtfusion_contigs} ")


if __name__ == "__main__":
	main()