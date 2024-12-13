import argparse as ap
from pathlib import Path

import pandas as pd

import re

from Bio.Seq import Seq

import hgvs.parser
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import hgvs.exceptions

parser = ap.ArgumentParser(
    description="Check all positions in a given transcript for genome-transcript discrepancies."
)
parser.add_argument(
    "infile",
    type=str,
    nargs=1,
    help="TSV file of genes and reference transcripts to analyze. First 3 columns must be: [id string] [transcript acc] [chr acc] ...",
)


def uta_tx_mapping_options_df(hdp, tx_ac):
    """
    Get the mapping options for a transcript accession from UTA and return as a DataFrame.

    If querying UTA directly, similar information would be retrieved as follows:

    import psycopg2, psycopg2.extras
    res = psycopg2.extras.DictCursor()
    conn = psycopg2.connect(
        host="localhost",
        dbname="uta",
        user="anonymous",
    )
    cur.execute(
         "select * from uta_20210129b.associated_accessions where tx_ac = %s", (tx_ac,)
    )
    res = cur.fetchall()

    """
    mapopts = hdp.get_tx_mapping_options(tx_ac)
    return pd.DataFrame(mapopts, columns=["tx_ac", "alt_ac", "method"])


def uta_get_similar_tx_df(hdp, tx_ac):
    """
    Check similar transcripts to see if UTA has a comparable transcript https://hgvs.readthedocs.io/en/1.4.0/modules/dataproviders.html#hgvs.dataproviders.uta.UTABase.get_similar_transcripts and return as a DataFrame.

    If querying UTA directly, similar information would be retrieved as follows:

    import psycopg2, psycopg2.extras
    res = psycopg2.extras.DictCursor()
    conn = psycopg2.connect(
        host="localhost",
        dbname="uta",
        user="anonymous",
    )
    cur.execute(
         "select * from tx_similarity_v where tx_ac1 ~ %s and es_fp_eq = 't'", (tx_ac,)
    )
    res = cur.fetchall()

    """
    similar_tx_res = hdp.get_similar_transcripts(tx_ac)
    return pd.DataFrame(
        similar_tx_res,
        columns=[
            "tx_ac1",
            "tx_ac2",
            "hgnc_eq",
            "cds_eq",
            "es_fp_eq",
            "cds_es_fp_eq",
            "cds_exon_lengths_fp_eq",
        ],
    )


def uta_get_tx_exons_df(hdp, tx_ac, alt_ac, alt_aln_method):
    """
    Get the exons for a transcript accession and alternate accession from UTA and return as a DataFrame.


    If querying UTA directly, similar information would be retrieved as follows:

    import psycopg2, psycopg2.extras
    res = psycopg2.extras.DictCursor()
    conn = psycopg2.connect(
        host="localhost",
        dbname="uta",
        user="anonymous",
    )
    cur.execute(
        "select * from uta_20210129b.tx_exon_aln_v where tx_ac = %s and alt_ac = %s",
        (tx_ac, chr_ac),
    )
    res = cur.fetchall()

    """
    txex = hdp.get_tx_exons(tx_ac, alt_ac, alt_aln_method)
    return pd.DataFrame(
        txex,
        columns=[
            "gene",
            "tx_ac",
            "alt_ac",
            "alt_aln_method",
            "alt_strand",
            "ord",
            "tx_start_i",
            "tx_end_i",
            "alt_start_i",
            "alt_end_i",
            "cigar",
            "unknown1",
            "unknown2",
            "tes_exon_set_id",
            "aes_exon_set_id",
            "tx_exon_id",
            "alt_exon_id",
            "unknown3",
        ],
    )


def uta_cigar_to_mismatch_vcf(hdp, id, row):
    # Initialize dataframe of detected mismatches
    mm = pd.DataFrame(
        columns=[
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "INFO",
        ],
    )

    # Get all match groups
    alngrps = [
        {
            "cigar_str_start": m.start(),
            "cigar_str_end": m.end(),
            "cigar_len": m.group(1),
            "cigar_op": m.group(2),
        }
        for m in re.finditer(r"([0-9]+)([A-Z=])", row["cigar"])
    ]
    tx_ac = row["tx_ac"]
    chr_ac = row["alt_ac"]
    alt_aln_method = row["alt_aln_method"]
    tx_cursor_i = row["tx_start_i"]
    chr_cursor_i = row["alt_start_i"] if row["alt_strand"] == 1 else row["alt_end_i"]
    # Iterate through the alignment groups. For each group:
    contiguous_delins = False
    for m in alngrps:
        if m["cigar_op"] == "=" or m["cigar_op"] == "M":
            contiguous_delins = False
            tx_cursor_i += int(m["cigar_len"])
            chr_cursor_i += (
                int(m["cigar_len"])
                if row["alt_strand"] == 1
                else (-int(m["cigar_len"]))
            )
            continue
        else:
            # If the previous iteration was a mismatch, skip processing the rest of the exon because the VCF will be incorrect
            if contiguous_delins:
                print(
                    f"Can't derive VCF for exon {row["ord"]} (CIGAR {row['cigar']}) for tx {tx_ac}, chr {chr_ac}, alt_aln_method {alt_aln_method}"
                )
                break
            if m["cigar_op"] == "X":
                contiguous_delins = False
                tx_cursor_i_new = tx_cursor_i + int(m["cigar_len"])
                chr_cursor_i_new = chr_cursor_i + (
                    int(m["cigar_len"])
                    if row["alt_strand"] == 1
                    else (-int(m["cigar_len"]))
                )
                tx_anchor_offset = 0
                chr_anchor_offset = 0
                chr_cursor_vcf_pos_3p_offset = 1
            elif m["cigar_op"] == "I":
                contiguous_delins = True
                tx_cursor_i_new = tx_cursor_i
                chr_cursor_i_new = chr_cursor_i + (
                    int(m["cigar_len"])
                    if row["alt_strand"] == 1
                    else (-int(m["cigar_len"]))
                )
                tx_anchor_offset = 1
                chr_anchor_offset = 1
                chr_cursor_vcf_pos_3p_offset = 0
            elif m["cigar_op"] == "D":
                contiguous_delins = True
                tx_cursor_i_new = tx_cursor_i + int(m["cigar_len"])
                chr_cursor_i_new = chr_cursor_i
                tx_anchor_offset = 1
                chr_anchor_offset = 1
                chr_cursor_vcf_pos_3p_offset = 0
            else:
                raise ValueError(
                    f"Unexpected CIGAR operation: {m['cigar_op']} for tx {tx_ac}, chr {chr_ac}, alt_aln_method {alt_aln_method}"
                )
            tx_mm_seq = (
                hdp.get_seq(tx_ac, tx_cursor_i - tx_anchor_offset, tx_cursor_i_new)
                if row["alt_strand"] == 1
                else hdp.get_seq(tx_ac, tx_cursor_i, tx_cursor_i_new + tx_anchor_offset)
            )
            chr_mm_seq = (
                hdp.get_seq(
                    chr_ac,
                    chr_cursor_i - chr_anchor_offset,
                    chr_cursor_i_new,
                )
                if row["alt_strand"] == 1
                else hdp.get_seq(
                    chr_ac,
                    chr_cursor_i_new - chr_anchor_offset,
                    chr_cursor_i,
                )
            )
            tx_pos = tx_cursor_i
            vcf_pos = (
                chr_cursor_i + chr_cursor_vcf_pos_3p_offset
                if row["alt_strand"] == 1
                else chr_cursor_i_new + chr_cursor_vcf_pos_3p_offset
            )
            vcf_ref = chr_mm_seq
            vcf_alt = str(
                tx_mm_seq
                if row["alt_strand"] == 1
                else Seq(tx_mm_seq).reverse_complement()
            )
            mismatch = pd.DataFrame(
                {
                    "#CHROM": chr_ac,
                    "POS": vcf_pos,
                    "ID": f"{chr_ac}|{vcf_pos}{vcf_ref}>{vcf_alt}|{tx_ac}|{alt_aln_method}",
                    "REF": vcf_ref,
                    "ALT": vcf_alt,
                    "INFO": f"tx_ac={tx_ac};cigar='{row["cigar"]}';alt_aln_method={alt_aln_method};uta_tx_exon_ord={row['ord']};uta_tx_exon_id={row['tx_exon_id']};uta_alt_exon_id={row['alt_exon_id']};uta_tx_start_i={row['tx_start_i']};uta_tx_end_i={row['tx_end_i']};tx_pos={tx_pos};uta_alt_start_i={row['alt_start_i']};uta_alt_end_i={row['alt_end_i']};strand={row['alt_strand']}",
                },
                index=[id],
            )
            mm = pd.concat(
                [
                    mm,
                    mismatch,
                ],
                ignore_index=False,
            )
        # Advance the cursor for the next iteration
        tx_cursor_i = tx_cursor_i_new
        chr_cursor_i = chr_cursor_i_new
    return mm


def main():
    args = parser.parse_args()

    infile = args.infile[0]  # 'mane_grch38_txlist.tsv'
    outfilebase = Path(infile).stem

    outvcf = outfilebase + ".mismatches.vcf"
    outfile = outfilebase + ".mismatches.tsv"

    txlist = pd.read_csv(infile, sep="\t", index_col=0)

    txlist = txlist.rename(
        {
            txlist.columns[0]: "tx_ac",
            txlist.columns[1]: "chr_ac",
        },
        axis=1,
    )

    # Initialize columns for results
    txlist["has_aln"] = None
    txlist["mismatch_exons"] = None
    txlist["mismatches"] = None
    txlist["errors"] = None

    # Initialize UTA connection
    hdp = hgvs.dataproviders.uta.connect()

    # Initialize dataframe of detected mismatches
    mm = pd.DataFrame(
        columns=[
            "#CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "INFO",
        ],
    )

    for id, row in txlist.iterrows():
        print(f"Now processing: {id}")
        tx_ac = row["tx_ac"]
        chr_ac = row["chr_ac"]
        # Check to see if target transcript is in UTA
        mapoptsdf = uta_tx_mapping_options_df(hdp, tx_ac)
        # if not res:
        if chr_ac not in mapoptsdf["alt_ac"].values:
            txlist.loc[id, "has_aln"] = False
            continue
        txlist.loc[id, "has_aln"] = True
        # Get transcript sequence
        for alt_aln_method in mapoptsdf[mapoptsdf["alt_ac"] == chr_ac]["method"]:
            txexdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
            # Get rows that don't have a perfect alignment according to CIGAR string
            if txexdf.empty:
                print(
                    f"No exons found for tx {tx_ac}, chr {chr_ac}, alt_aln_method {alt_aln_method}"
                )
                continue
            # For each exon that isn't perfectly aligned to the reference
            for i, row in txexdf[
                ~txexdf["cigar"].str.fullmatch("^[0-9]+=$")
            ].iterrows():
                mm = pd.concat(
                    [
                        mm,
                        uta_cigar_to_mismatch_vcf(hdp, id, row),
                    ],
                    ignore_index=False,
                )
    if not mm.empty:
        # Get the unique 0-based exon ordinal numbers for each mismatch in this transcript as ';'-delimited string
        txlist["mismatch_exons"] = (
            mm.groupby(mm.index)["INFO"]
            .agg("|".join)
            .str.extractall("uta_tx_exon_ord=([0-9]+);")
            .reset_index(level=0, names=["id"])
            .drop_duplicates(["id", 0])
            .sort_values(["id", 0])
            .groupby("id")
            .agg(";".join)
        )
        # Get the unique mismatch IDs for each mismatch in this transcript as ';'-delimited string
        txlist["mismatches"] = mm.groupby(mm.index)["ID"].agg(";".join)
    # Write output
    mm.to_csv(outvcf, sep="\t", index=False)
    txlist.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    main()
