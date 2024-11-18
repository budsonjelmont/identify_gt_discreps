# import pytest
from find_mismatch_positions_v2 import uta_cigar_to_mismatch_vcf, uta_get_tx_exons_df

import hgvs.dataproviders.uta

# Initialize UTA connection
hdp = hgvs.dataproviders.uta.connect()  # uta_20210129b used in development


def test_single_mismatch_pos_strand():
    """
    Scenario: Single mismatch
    CIGAR: 47=1X195=
    """
    tx_ac = "NM_000090.3"
    chr_ac = "NC_000002.12"
    alt_aln_method = "splign"
    ex_ord = 49
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 189010695
    expect_ref = "T"
    expect_alt = "G"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_single_mismatch_pos_strand_2():
    """
    Scenario: Single mismatch
    CIGAR: 389=1X38=
    """
    tx_ac = "NM_000059.3"
    chr_ac = "NC_000013.11"
    alt_aln_method = "splign"
    ex_ord = 13
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 32355250
    expect_ref = "T"
    expect_alt = "C"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_single_mismatch_min_strand():
    """
    Scenario: Single mismatch
    CIGAR: 2720=1X4851=
    """
    tx_ac = "NM_000384.2"
    chr_ac = "NC_000002.12"
    alt_aln_method = "splign"
    ex_ord = 25
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 21009931
    expect_ref = "T"
    expect_alt = "C"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_single_mismatch_min_strand_2():
    """
    Scenario: Single mismatch
    CIGAR: 204=1X10=
    """
    tx_ac = "NM_000130.4"
    chr_ac = "NC_000001.10"
    alt_aln_method = "splign"
    ex_ord = 9
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 169519049
    expect_ref = "T"
    expect_alt = "C"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_two_contig_mismatches_min_strand():
    """
    TODO - Implement
    Scenario: Two contiguous single bp events
    NM_033487.3
    NC_000001.10:g.1586821_1586822delACinsCT (NM_033487.3:c.-281_-280delAGinsGT)
    tx_ac = "NM_033487.3"
    chr_ac = "NC_000001.10"

    See also:
    NC_000004.11:g.76676624_76676625delCTinsTA	NM_003715.4:c.206_207delTAinsCT
    NC_000004.11:g.76676624_76676625delCTinsTA	NM_001290049.2:c.206_207delTAinsCT
    """
    assert 1 == 1


def test_single_bp_del_pos_strand():
    """
    Scenario: Single mismatch
    CIGAR: 980=1D2=
    """
    tx_ac = "NM_014487.4"
    chr_ac = "NC_000004.11"
    alt_aln_method = "blat"
    ex_ord = 9
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 142155848
    expect_ref = "T"
    expect_alt = "TA"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_single_bp_del_min_strand():
    """
    # Scenario: Single bp deletion
    # CIGAR: 2410=1D2=
    """
    tx_ac = "NM_004657.5"
    chr_ac = "NC_000002.11"
    alt_aln_method = "blat"
    ex_ord = 1
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 192699033
    expect_ref = "T"
    expect_alt = "TT"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_multi_bp_del_pos_strand():
    """
    Scenario: Multiple bp deletion
    CIGAR: 4=9D149=
    g. HGVS: NC_000014.8:g.94582133_94582134insCATGGCGGC
    c. HGVS: NM_001366994.1:c.129_137delCATGGCGGC
    """
    tx_ac = "NM_001366994.1"
    chr_ac = "NC_000014.8"
    alt_aln_method = "splign"
    ex_ord = 3
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 94582130
    expect_ref = "T"
    expect_alt = "TGGCCATGGC"  # Variant Recoder has "TGGCTGTGCC" as expected alt allele given g. HGVS above, but the discrepancy is likely due to normalization

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_multi_bp_del_pos_strand_2():
    """
    Scenario: Another multiple bp deletion
    CIGAR: 1453=3D2=
    """
    tx_ac = "NM_001256326.1"
    chr_ac = "NC_000017.10"
    alt_aln_method = "blat"
    ex_ord = 35
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 48704830
    expect_ref = "T"
    expect_alt = "TAAA"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_multi_bp_del_min_strand():
    """
    Scenario: Multiple bp deletion
    CIGAR: 459=14D1318=
    HGVSg: NC_000017.10:g.73234011_73234012insTCCCACCCCCCACC
    HGVSc: NM_001291642.2:c.*350_*363delGTGGGGGGTGGGAG
    """
    tx_ac = "NM_001291642.2"
    chr_ac = "NC_000017.10"
    ex_ord = 13
    alt_aln_method = "splign"
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 73234011
    expect_ref = "C"
    expect_alt = "CTCCCACCCCCCACC"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_single_bp_ins_pos_strand():
    """
    Scenario: Single bp insertion
    CIGAR: 136=1I129=
    HGVSg: NC_000004.11:g.39046587delA
    HGVSc: NM_001171654.1:c.-170_-169insA
    """
    tx_ac = "NM_001171654.1"
    chr_ac = "NC_000004.11"
    alt_aln_method = "splign"
    ex_ord = 0
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 39046586
    expect_ref = "CA"
    expect_alt = "C"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_multi_bp_ins_pos_strand():
    """
    Scenario: Multiple bp insertion
    CIGAR: 14=1X11=6I26=
    HGVSg: NC_000001.10:g.145328407_145328412delGAAGAC (note the 3' normalization)
    HGVSc: NM_001039703.6:c.4474_4479dupGAAGAC
    """
    tx_ac = "NM_001039703.6"
    chr_ac = "NC_000001.10"
    alt_aln_method = "splign"
    ex_ord = 34
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 145328400
    expect_ref = "AGAAGAC"
    expect_alt = "A"

    assert resultdf.shape[0] == 2
    resultdf = resultdf.iloc[1]
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_multi_bp_ins_pos_strand_2():  # Old & maybe wrong
    """
    Scenario: Another multiple bp insertion
    CIGAR: 284=3I1=
    """
    tx_ac = "NM_001788.5"
    chr_ac = "NC_000007.13"
    alt_aln_method = "blat"
    ex_ord = 0
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 35840879
    expect_ref = "GGGT"
    expect_alt = "G"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_single_bp_ins_neg_strand():
    """
    Scenario: Single bp insertion
    CIGAR: 428=1I76=
    """
    tx_ac = "NM_001160329.1"
    chr_ac = "NC_000019.9"
    alt_aln_method = "splign"
    ex_ord = 10
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 51125309
    expect_ref = "GG"
    expect_alt = "G"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_multi_bp_ins_neg_strand():
    """
    Scenario: Multiple bp insertion
    CIGAR: 52=6I14=
    """
    tx_ac = "NM_001290207.2"
    chr_ac = "NC_000003.11"
    alt_aln_method = "splign"
    ex_ord = 5
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_pos = 53324828
    expect_ref = "GCACTGG"
    expect_alt = "G"

    assert resultdf.shape[0] == 1
    resultdf = resultdf.squeeze()
    assert resultdf["POS"] == expect_pos
    assert resultdf["REF"] == expect_ref
    assert resultdf["ALT"] == expect_alt


def test_two_noncontig_single_bp_events():
    """
    Scenario: Two non-contiguous events, an insertion and a mismatch
    CIGAR: 666=1I39=1X404=
    """
    tx_ac = "NM_000314.4"  # PTEN
    chr_ac = "NC_000010.10"
    alt_aln_method = "splign"
    ex_ord = 0
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_ins_pos = 89623860
    expect_ins_ref = "CT"
    expect_ins_alt = "C"

    expect_mm_pos = 89623901
    expect_mm_ref = "G"
    expect_mm_alt = "C"

    assert resultdf.shape[0] == 2
    assert resultdf.iloc[0]["POS"] == expect_ins_pos
    assert resultdf.iloc[0]["REF"] == expect_ins_ref
    assert resultdf.iloc[0]["ALT"] == expect_ins_alt

    assert resultdf.iloc[1]["POS"] == expect_mm_pos
    assert resultdf.iloc[1]["REF"] == expect_mm_ref
    assert resultdf.iloc[1]["ALT"] == expect_mm_alt


def test_multiple_indels_min_strand():
    """
    TODO - Implement
    Scenario: Multiple indels
    CIGAR: 498=1D37=3I1809=
    HGVSg: NC_000001.10:g.150192526_150192528delAAA and NC_000001.10:g.150192565_150192566insA
    HGVSc: NM_001280560.2:c.*570_*572dupTTT and NM_001280560.2:c.*516delT
    """
    tx_ac = "NM_001280560.2"
    chr_ac = "NC_000001.10"
    alt_aln_method = "splign"
    ex_ord = 4
    exdf = uta_get_tx_exons_df(hdp, tx_ac, chr_ac, alt_aln_method)
    row = exdf.loc[exdf["ord"] == ex_ord].squeeze()
    resultdf = uta_cigar_to_mismatch_vcf(hdp, "ABC", row)

    expect_del_pos = 150192565
    expect_del_ref = "C"
    expect_del_alt = "CA"

    expect_ins_pos = 150192065
    expect_ins_ref = "CAAA"
    expect_ins_alt = "C"

    assert resultdf.shape[0] == 2
    assert resultdf.iloc[0]["POS"] == expect_del_pos
    assert resultdf.iloc[0]["REF"] == expect_del_ref
    assert resultdf.iloc[0]["ALT"] == expect_del_alt

    assert resultdf.iloc[1]["POS"] == expect_ins_pos
    assert resultdf.iloc[1]["REF"] == expect_ins_ref
    assert resultdf.iloc[1]["ALT"] == expect_ins_alt


if __name__ == "__main__":
    test_single_bp_del_min_strand()
