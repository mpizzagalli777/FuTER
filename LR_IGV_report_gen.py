#!/usr/bin/env python3

import argparse, csv, os, sys, subprocess, shutil, json
from pathlib import Path

try:
    csv.field_size_limit(sys.maxsize)
except OverflowError:
    csv.field_size_limit(10**9)

def run(cmd, **kw):
    print(f"[cmd] {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True, **kw)

def ensure(parent: Path):
    parent.mkdir(parents=True, exist_ok=True)
    return parent

def read_fusion_names(fusion_list_path, fusion_args):
    names = set()
    if fusion_list_path:
        with open(fusion_list_path) as fh:
            for line in fh:
                s = line.strip()
                if s and not s.startswith("#"):
                    names.add(s)
    for s in (fusion_args or []):
        names.add(s)
    if not names:
        sys.exit("ERROR: no fusion names provided.")
    return names

def filter_fusions_table(in_tsv: Path, keep_names: set, out_tsv: Path):
    # Header is typically "#FusionName"
    key_candidates = ["#FusionName", "FusionName", "Fusion", "FusionName"]
    with open(in_tsv, newline="") as infh:
        reader = csv.DictReader(infh, delimiter="\t")
        hdr = reader.fieldnames or []
        key = next((k for k in key_candidates if k in hdr), None)
        if key is None:
            sys.exit(f"ERROR: cannot find a fusion-name column in {in_tsv}. "
                     f"Looked for {key_candidates}. Found: {hdr}")
        rows = [r for r in reader if r.get(key) in keep_names]
    if not rows:
        sys.exit("ERROR: none of the requested fusion names were found in the fusion table.")
    with open(out_tsv, "w", newline="") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=hdr, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return out_tsv
    
def filter_fusions_json(json_path: Path, keep_names: set):
    """
    Keep only entries whose fusion name matches keep_names.
    Handles both list- and dict-wrapped schemas from CTAT utilities.
    """
    with open(json_path) as fh:
        data = json.load(fh)

    def fusion_name(obj):
        for k in ("fusion_name", "FusionName", "#FusionName", "fusion", "name", "Fusion"):
            if isinstance(obj, dict) and k in obj:
                return obj[k]
        return None

    def filter_list(lst):
        return [o for o in lst if fusion_name(o) in keep_names]

    changed = False
    if isinstance(data, list):
        data = filter_list(data); changed = True
    elif isinstance(data, dict):
        for key in ("fusions", "items", "rows", "features"):
            if key in data and isinstance(data[key], list):
                data[key] = filter_list(data[key]); changed = True

    if changed:
        with open(json_path, "w") as out:
            json.dump(data, out, indent=2)
    else:
        print(f"[warn] {json_path} did not look like a fusion list JSON; left unfiltered.", file=sys.stderr)

def filter_roi_bed(roi_bed: Path, keep_names: set):
    """
    BED from CTAT has fusion name as column 4; keep only those in keep_names.
    Safe to no-op if format differs.
    """
    if not roi_bed.exists():
        return
    out = roi_bed.with_suffix(".filtered.bed")
    kept = []
    with open(roi_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            name = fields[0] if len(fields) >= 3 else None
            if name in keep_names:
                kept.append(line)
    if kept:
        with open(out, "w") as ofh:
            ofh.writelines(kept)
        roi_bed.unlink(missing_ok=True)
        out.rename(roi_bed)
    else:
        print(f"[warn] No ROI entries matched keep list in {roi_bed}; leaving as-is.", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description="Generate an IGV-Reports HTML for a selected set of fusions "
                    "(mirrors CTAT-LR-fusion include_IGV_REPORTS, but filtered).")
    ap.add_argument("--Phase3_TELR_Fusion_dir", required=True, help="Path to directory containing *_LR-FI_targets.gtf, *_LR-FI.mm2.bam, etc")
    ap.add_argument("--Sample_prefix",        required=True, help="Prefix of sample of interest. Used to find (*_LR-FI_targets.gtf, *_LR-FI.mm2.bam, etc)")

    ap.add_argument("--fusion_list",            help="Text file: one fusion name per line (e.g. ENSG...--chrTE_XXX)")
    ap.add_argument("--fusion", action="append", help="Fusion name (may be given multiple times)")

    ap.add_argument("--max_LR_per_fusion", type=int, default=20, help="Cap per-fusion LR read count (default: 20)")
    ap.add_argument("--threads", type=int, default=8, help="Threads for samtools (default: 8)")

    ap.add_argument("--outdir", required=True, help="Output directory. This will create IGV_prep inside that contains intermediate files and a Final_output folder. If the same directory is used as during the original FuTER run, all final outputs will be in the same folder.")

    # Optional short-read tracks (if you ran FI on Illumina)
    ap.add_argument("--fi-shortreads-dir", help="FusionInspector fi_workdir (will look for *fusion_junc_reads.sam / *fusion_span_reads.sam). This is currently not supported.")

    # Optional explicit tracks.json
    ap.add_argument("--tracks-json", help="Custom tracks.json to use instead of defaults")

    args = ap.parse_args()
    
    sample_prefix = args.Sample_prefix
    phase3_dir	  = Path(args.Phase3_TELR_Fusion_dir).expanduser().resolve()
    
    fi_fa   = phase3_dir / f"{sample_prefix}_LR-FI_targets.fa"
    fi_gtf  = phase3_dir / f"{sample_prefix}_LR-FI_targets.gtf"
    lr_bam  = phase3_dir / f"{sample_prefix}_LR-FI.mm2.bam"
    fusions = phase3_dir / f"{sample_prefix}_LR-FI.mm2.fusion_transcripts.breakpoint_info.tsv.w_LR_FFPM"

    script_dir = Path(__file__).resolve().parent #Gives path to FuTER base directory

    util_dir    	= ensure(script_dir / "utils/scripts")
    genome_lib_dir 	= ensure(script_dir / "TEgenome_prep")
    outdir     		= Path(args.outdir).expanduser().resolve()
    igv_dir    		= ensure(outdir / "IGV_prep")
    igv_samp_dir	= ensure(igv_dir / f"{sample_prefix}")
    Final_dir     	= ensure(outdir / "Final_output")
    
    
    igv_dir.mkdir(parents=True, exist_ok=True)
    igv_samp_dir.mkdir(parents=True, exist_ok=True)
	
    # 0) pick fusion names
    keep = read_fusion_names(args.fusion_list, args.fusion)
    filtered_tsv = igv_samp_dir / f"{sample_prefix}_selected_fusions.tsv"
    filter_fusions_table(fusions, keep, filtered_tsv)

    # 1) include the FI contigs + annotations
    (igv_samp_dir / "igv.genome.fa").unlink(missing_ok=True)
    (igv_samp_dir / "igv.genome.fa.fai").unlink(missing_ok=True)
    (igv_samp_dir / "igv.annot.gtf").unlink(missing_ok=True)

    os.symlink(fi_fa,  igv_samp_dir / "igv.genome.fa")
    # ensure FAI exists
    fai = fi_fa.with_name(fi_fa.name + ".fai")
    if not fai.exists():
        run(f"samtools faidx {fi_fa}")
    os.symlink(fai,    igv_samp_dir / "igv.genome.fa.fai")
    os.symlink(fi_gtf, igv_samp_dir / "igv.annot.gtf")

    # gtf -> bed for track convenience
    run(f"{util_dir}/gtf_gene_to_bed.pl {igv_samp_dir}/igv.annot.gtf > {igv_samp_dir}/igv.annot.bed")

    # 2) select LR reads supporting just these fusions
    sel_sam = igv_samp_dir / f"{sample_prefix}_LR-FI.mm2.max_per_fusion-{args.max_LR_per_fusion}.sam"
    run(
        f"{util_dir}/LR_sam_fusion_read_extractor.pl "
        f"--FI_LR_sam {lr_bam} "
        f"--LR_fusion_report {filtered_tsv} "
        f"--max_alignments_per_fusion {args.max_LR_per_fusion} "
        f"> {sel_sam}"
    )

    # 3) convert to BAM, sort, index
    igv_lr_bam = igv_samp_dir / f"igv.LR.bam"
    igv_lr_sorted = igv_samp_dir / f"igv.LR.sorted.bam"
    run(f"samtools view -@ {args.threads} -Sb {sel_sam} -o {igv_lr_bam}")
    run(f"samtools sort -@ {args.threads} {igv_lr_bam} -o {igv_lr_sorted}")
    run(f"samtools index {igv_lr_sorted}")

    # 4) Optional Illumina tracks
    have_illumina = False
    if args.fi_shortreads_dir:
        sr = Path(args.fi_shortreads_dir)
        junc = sr / "finspector.star.cSorted.dupsMarked.bam.fusion_junc_reads.sam"
        span = sr / "finspector.star.cSorted.dupsMarked.bam.fusion_span_reads.sam"
        if junc.exists():
            out = igv_samp_dir / "igv.illumina.junction_reads.bam"
            run(f"samtools view -Sb {junc} -T {fi_fa} -o {out} && "
                f"samtools sort {out} -o {igv_samp_dir}/igv.illumina.junction_reads.sorted.bam && "
                f"samtools index {igv_samp_dir}/igv.illumina.junction_reads.sorted.bam && "
                f"rm -f {out}")
            have_illumina = True
        if span.exists():
            out = igv_samp_dir / "igv.illumina.spanning_frags.bam"
            run(f"samtools view -Sb {span} -T {fi_fa} -o {out} && "
                f"samtools sort {out} -o {igv_samp_dir}/igv.illumina.spanning_frags.sorted.bam && "
                f"samtools index {igv_samp_dir}/igv.illumina.spanning_frags.sorted.bam && "
                f"rm -f {out}")
            have_illumina = True

    # 5) Pfam & seq-similar region tracks
    pfam_gff3 = igv_samp_dir / f"igv.pfam.gff3"
    run(f"{util_dir}/get_pfam_domain_info.pl "
        f"--finspector_gtf {fi_gtf} "
        f"--genome_lib_dir {genome_lib_dir} "
        f"> {pfam_gff3}")
    run(f"{util_dir}/transcript_gff3_to_bed.pl {pfam_gff3} > {igv_samp_dir}/igv.pfam.bed")

    seqsim_gff3 = igv_samp_dir / f"{sample_prefix}_igv.seqsimilar.gff3"
    run(f"{util_dir}/get_seq_similar_region_FI_coordinates.pl "
        f"--finspector_gtf {fi_gtf} "
        f"--genome_lib_dir {genome_lib_dir} "
        f"> {seqsim_gff3}")
    run(f"{util_dir}/transcript_gff3_to_bed.pl {seqsim_gff3} > {igv_samp_dir}/igv.seqsimilar.bed")

    # 6) Build the fusion JSON & ROI for igv-reports
    fusions_json = igv_samp_dir / f"igv.fusion_inspector_web.json"
    roi_bed      = igv_samp_dir / f"igv.LR.breakoint.roi.bed"
    run(
        f"{args.util_dir}/create_ctat-LR-fusion_inspector_igvjs.py "
        f"--fusion_inspector_directory {Final_dir}/Final_output "
        f"--json_outfile {fusions_json} "
        f"--roi_outfile {roi_bed} "
        f"--sample_prefix {sample_prefix} "        
        f"--file_prefix ctat-LR-fusion"
    )
    
    filter_fusions_json(fusions_json, keep)
    filter_roi_bed(roi_bed, keep)

    # 7) tracks.json (default with/without Illumina)
    tracks_src = None
    if args.tracks_json:
        tracks_src = Path(args.tracks_json)
    else:
        base = ensure(util_dir / "fusion_report_html_template")
        tracks_src = base / ("tracks.wIllumina.json" if have_illumina else "tracks.json")
    tracks_dst = igv_samp_dir / "tracks.json"
    shutil.copyfile(tracks_src, tracks_dst)

    # 8) Render IGV-Reports HTML
    out_html = Final_dir / f"{sample_prefix}.fusion_inspector_web.html"
    run(
        f"cd {igv_samp_dir} && "
        f"create_report {fusions_json} igv.genome.fa "
        f"--type fusion --track-config {tracks_dst} "
        f"--output {out_html}"
    )

    print(f"\n[done] IGV report: {out_html}\n")
    

    
if __name__ == "__main__":
    main()
