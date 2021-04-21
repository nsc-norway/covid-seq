#!/usr/bin/env python

import glob
import sys
import os
import multiprocessing
import vcf
import pysam
import pandas as pd
import functools

# Variant list file format:
# Name,Type,Pos,Ref,Alt,DelLength

#---
CONTIG_NAME = "NC_045512.2"
MIN_AF_CUT = 0.75
MIN_DEPTH_CUT = 10
DEPTH_BASE_QUAL_THRESHOLD = 20

class MutationStatus:
    """Class to hold the detection and QC status of a variant, SNP or deletion."""

    def __init__(self, detected=None, qc_problem=None, failure=False,
            qc_status_message=None, total_depth=None, alt_depth=None
            ):
        self.detected = detected
        self.qc_problem = qc_problem
        self.failure = failure
        self.qc_status_message = qc_status_message
        self.total_depth = total_depth
        self.alt_depth = alt_depth

    def get_table_message(self):
        if self.failure:
            return "FAIL: {}".format(self.qc_status_message)
        elif self.detected:
            if self.qc_problem: return "CHECK: {}".format(self.qc_status_message)
            else:               return "DETECTED"
        else:
            if self.qc_problem: return "FAIL: {}".format(self.qc_status_message)
            else:               return ""
    


def check_snp_at_pos(pos, ref, alts, ignore_alts, vcf_records, bam):
    """Check SNP at specified reference position, reference allele, alternative allele (mutation) in the
    given VCF data.
    
    Can take multiple expected ALT alleles, separated by pipe (|)."""

    if not bam:
        return MutationStatus(failure=True, qc_status_message="No BAM file.")

    # Get a tuple of coverages for A, C, G, T, each containing an array (of length same as our range, 1)
    acgt_covs = bam.count_coverage(CONTIG_NAME, pos-1, pos, quality_threshold=DEPTH_BASE_QUAL_THRESHOLD)
    total_depth = sum(nc[0] for nc in acgt_covs)
    alt_depth = sum(acgt_covs["ACGT".index(alt_nuc)][0] for alt_nuc in alts)
    alt_freq = alt_depth / max(total_depth,1)
    found_alt_alleles = set()
    if vcf_records is None:
        status = "MISSING"
    else:
        status = "UNKNOWN"
        for rec in vcf_records:
            if rec.POS == pos:
                if rec.is_snp and (rec.FILTER == [] or rec.FILTER == ["PASS"]):
                    if rec.REF != ref: raise ValueError("Reference allele in VCF '{}' doesn't match "
                                                        "reference allele '{}' in variant list (at pos {}).".format(
                                                            rec.REF, ref, pos
                                                        ))

                    found_alt_alleles.add(str(rec.ALT[0]))
                    status = "DETECTED"
                    break

    if status == "DETECTED":
        # SNP FOUND -- Check QC status
        if not (found_alt_alleles & alts): # No overlap of found alleles and expected alleles
            if found_alt_alleles.issubset(ignore_alts): # .. but we are supposed to ignore this one
                return MutationStatus(detected=False, qc_problem=False, total_depth=total_depth,
                    alt_depth=alt_depth)
            else:
                return MutationStatus(detected=True, qc_problem=True,
                    qc_status_message="Unexpected {} at {} (REF={},Expect={}).".format(
                            ",".join(found_alt_alleles), pos, ref, ",".join(found_alt_alleles)),
                    total_depth=total_depth, alt_depth=alt_depth
                )
        elif alt_freq < MIN_AF_CUT:
            return MutationStatus(detected=True, qc_problem=True,
                qc_status_message="Low fraction of reads: {}.".format(alt_freq),
                total_depth=total_depth, alt_depth=alt_depth
            )
        elif total_depth < MIN_DEPTH_CUT:
            return MutationStatus(detected=True, qc_problem=True,
                qc_status_message="Low read depth: {} (SNP found).".format(total_depth),
                total_depth=total_depth, alt_depth=alt_depth
            )
        else:
            return MutationStatus(detected=True, qc_problem=False,
                total_depth=total_depth, alt_depth=alt_depth
            )
    else:
        # NO SNP FOUND
        if status == "MISSING":
            return MutationStatus(failure=True, qc_status_message="Missing VCF file.")
        elif total_depth < MIN_DEPTH_CUT:
            return MutationStatus(detected=False, qc_problem=True,
                qc_status_message="Low read depth: {}.".format(total_depth),
                total_depth=total_depth, alt_depth=0
            )
        else:
            return MutationStatus(detected=False, qc_problem=False, total_depth=total_depth, alt_depth=0)


def check_deletion(startrange, endrange, del_length, vcf_records, bam):
    """Look for deletion in coordinate range."""

    if not bam:
        return MutationStatus(failure=True, qc_status_message="No BAM file.")

    # Check coverage in a range that extends 1 nucleotide outside the start/end range used to check for
    # deletion. Then we can get the coverage.
    acgt_covs = bam.count_coverage(CONTIG_NAME, startrange-3, endrange+1, quality_threshold=DEPTH_BASE_QUAL_THRESHOLD)
    depth_array = [sum(nd) for nd in zip(*acgt_covs)]
    total_depth = min(depth_array[0], depth_array[-1])
    ref_depth_at_deletion = min(depth_array[1:-1])
    alt_depth = max([total_depth - ref_depth_at_deletion, 0])
    alt_freq = alt_depth / max(1, total_depth)

    if max(depth_array) * 0.3 > total_depth:
        # Somehow reads at the start or end of interval are lower. THis could happen if there's a deletion at
        # the edge of the interval. The user may need to adjust search bounds for the deletion.
        return MutationStatus(failure=True,
                qc_status_message="Low coverage at edge of the range.",
                total_depth=total_depth, alt_depth=alt_depth)

    deletion_lengths = set()
    if vcf_records is None:
        return MutationStatus(failure=True, qc_status_message="Missing VCF file.")
    else:
        status = False
        for rec in vcf_records:
            if rec.affected_start >= startrange-1 and rec.affected_end <= endrange+1 and \
                    rec.is_deletion and (rec.FILTER == [] or rec.FILTER == ["PASS"]):
                deletion_lengths.add(rec.affected_end - rec.affected_start)
                status = True
                break
    if status:
        if deletion_lengths == set([del_length]):
            if alt_freq < MIN_AF_CUT:
                return MutationStatus(detected=True, qc_problem=True,
                    qc_status_message="Low fraction of reads: {}.".format(alt_freq),
                    total_depth=total_depth, alt_depth=alt_depth)
            elif total_depth < MIN_DEPTH_CUT:
                return MutationStatus(detected=True, qc_problem=True,
                        qc_status_message="Low read depth: {}.".format(total_depth),
                        total_depth=total_depth, alt_depth=alt_depth)
            else:
                return MutationStatus(detected=True, qc_problem=False,
                    total_depth=total_depth, alt_depth=alt_depth)
        else:
            return MutationStatus(detected=True, qc_problem=True,
                    qc_status_message="Deletion length(s): {}.".format(list(deletion_lengths)),
                    total_depth=total_depth, alt_depth=alt_depth)
    elif total_depth < MIN_DEPTH_CUT:
        return MutationStatus(detected=False, qc_problem=True,
                qc_status_message="Low read depth: {}.".format(total_depth),
                total_depth=total_depth, alt_depth=alt_depth)
    else:
        return MutationStatus(detected=False, total_depth=total_depth, alt_depth=alt_depth)


def get_vcfs(vcf_dir, compressed):
    """ Searches through specified VCF dirs, reads all vcf files and stores them in a dict indexed
    by variant caller and sample name.
    
    Returns a dict, or None if no VCF files are found."""

    if compressed: suffix = ".gz"
    else:          suffix = ""

    vcf_data = {}
    found = False
    for vcf_path in glob.glob(os.path.join(vcf_dir, "*_*.vcf" + suffix)):
        bas = os.path.basename(vcf_path)
        sample, _ = bas.split("_", maxsplit=1)
        reader = vcf.Reader(filename=vcf_path, compressed=compressed)
        vcf_data[sample] = list(reader)
        found = True
    if not found:   return None
    else:           return vcf_data


def process_sample(variants, bam_dir, vcf_data, sample_row):
    """Process a sample (called in parallel)"""

    out_row = {'Sample': sample_row['sample'], 'Well': sample_row['Well']}
    full_out_common = dict(out_row)
    full_out_rows = []

    bam_paths = glob.glob("{}/{}.*bam".format(bam_dir, sample_row['sample']))
    if bam_paths:
        bam = pysam.AlignmentFile(bam_paths[0])
    else:
        # The check_snp/deletion functions accept a None value for bam, and then return a QC failure status.
        bam = None

    # Get the VCF data for the current sample. If the VCF file was missing, the list will contain None.
    sample_vcf = vcf_data.get(sample_row['sample'])
    for _, variant in variants.iterrows():
        full_out_row = dict(full_out_common)
        full_out_row['Variant'] = variant.Name
        if variant.Type == "SNP":
            alts = set(variant.Alt.split('|')) if not pd.isna(variant.Alt) else set()
            ignore_alts = set(variant.IgnoreAlt.split('|')) if not pd.isna(variant.IgnoreAlt) else set()
            status = check_snp_at_pos(variant.Pos, variant.Ref, alts, ignore_alts, sample_vcf, bam)
            out_row[variant.Name] = status.get_table_message()
            for name, value in vars(status).items():
                full_out_row[name] = value
        elif variant.Type == "DELETION":
            status = check_deletion(variant.Pos, variant.DelRangeEnd, variant.DelLength, sample_vcf, bam)
            out_row[variant.Name] = status.get_table_message()
            for name, value in vars(status).items():
                full_out_row[name] = value
        else:
            print("Error: Unknown variant type '{}' for variant '{}'.".format(variant.Typem, variant.Name))
            sys.exit(1)
        full_out_rows.append(full_out_row)
    if bam:
        bam.close()
    return (out_row, full_out_rows)


def main(variants_list_path, sample_list_path, bam_dir, vcf_dir):
    
    samples = pd.read_csv(sample_list_path)
    variants = pd.read_csv(variants_list_path)

    vcf_data = get_vcfs(vcf_dir, compressed=True)
    if vcf_data is None:
        vcf_data = get_vcfs(vcf_dir, compressed=False)
    if vcf_data is None:
        raise RuntimeError("No VCF files found in {}.".format(vcf_dir))

    process_sample_p = functools.partial(
        process_sample, variants, bam_dir, vcf_data
    )

    pool = multiprocessing.Pool(len(os.sched_getaffinity(0)))
    results = pool.map(process_sample_p, [row for _, row in samples.iterrows()])
    #results = map(process_sample_p, [row for _, row in samples.iterrows()])
    outdata = [out_row for out_row, _ in results]
    pd.DataFrame(outdata).to_csv("variant_table.csv", index=False, sep=";")
    full_outdata = [fo_row for _, fo_rows in results for fo_row in fo_rows]
    pd.DataFrame(full_outdata).to_csv("detailed_variant_table.csv", index=False, sep=",")

if __name__ == "__main__":
    if len(sys.argv) > 4:
        main(*sys.argv[1:])
    else:
        print("Usage:")
        print(" check_variants.py VariantsList SampleList BamDir {VcfRootDirs}")
