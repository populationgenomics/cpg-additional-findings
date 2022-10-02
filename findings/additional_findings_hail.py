"""
a hail classifier specific to the additional findings process
"""

from argparse import ArgumentParser

import hail as hl

from cpg_utils.hail_batch import output_path
from findings.utils import read_json_from_path
from findings.hail_utils import (
    checkpoint_and_repartition,
    filter_to_population_rare,
    filter_on_quality_flags,
    split_rows_by_gene_and_filter_to_green,
    PATHOGENIC,
    CONFLICTING,
    ONE_INT,
    MISSING_INT,
    write_matrix_to_vcf,
)


# this will be a dictionary mapping the exact consequences with the VEP CSQs
# currently the only category that exists is 'truncating'
VARIANT_CONSEQUENCES = {
    'truncating': hl.set(
        [
            'frameshift_variant',
            'splice_acceptor_variant',
            'splice_donor_variant',
            'stop_gained',
        ]
    )
}


# not sure how best to structure this whole thing...
# either we take the MT, annotate it, and return it
# which seems increasingly wasteful as the search goals get hyper-specific
# or process the MT, get the loci, and annotate the flags back against those loci....


def add_top_level_clinvar_flag(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    bop bop
    Parameters
    ----------
    mt :

    Returns
    -------
    same MT with two fields = clinvar yes/no, and a keep==no
    """
    return mt.annotate_rows(
        clinvar=hl.if_else(
            (mt.info.clinvar_sig.lower().contains(PATHOGENIC))
            & ~(mt.info.clinvar_sig.lower().contains(CONFLICTING)),
            ONE_INT,
            MISSING_INT,
        ),
        keep=MISSING_INT,
    )


def read_and_filter_mt(
    mt_path: str, filter_af: float, green_genes: hl.SetExpression
) -> hl.MatrixTable:
    """
    import the MT from file, and run some standard filters
    recycling as much as possible, inc splitting CSQ by gene

    Parameters
    ----------
    mt_path : path to the MT to use for input
    filter_af : AF filter to use
    green_genes : ENSG values to filter to

    Returns
    -------
    the filtered MT

    TODO add a field audit call...
    TODO work out what would go in that audit...
    """

    mt = hl.read_matrix_table(mt_path)
    mt = filter_on_quality_flags(mt)
    mt = filter_to_population_rare(mt, filter_af)
    return split_rows_by_gene_and_filter_to_green(mt, green_genes=green_genes)


def find_standard_interesting(mt: hl.MatrixTable, acmg_data: dict[str, dict]):
    """
    finds all variants which are:
        - in genes flagged by the ACMG specification
        - have Clinvar annotations of interest
        - are not restricted to only specific variant sub-types

    Parameters
    ----------
    mt : The MatrixTable of all variants
    acmg_data : Parsed ACMG data

    Returns
    -------
    same MT, with 'simple' clinvar P/LP flagged as 'keep'
    """

    plain_ids = hl.set(
        [
            key
            for key, value in acmg_data.items()
            if not any([value.get('specific_type'), value.get('specific_variant')])
        ]
    )

    return mt.annotate_rows(
        keep=hl.if_else(
            (plain_ids.contains(mt.geneIds)) & (mt.clinvar == ONE_INT),
            ONE_INT,
            mt.keep,
        )
    )


def find_specific_types(
    mt: hl.MatrixTable, acmg_data: dict[str, dict]
) -> hl.MatrixTable:
    """
    Finds all variants in the ACMG specification where the variant:
        - in genes flagged by the ACMG specification
        - have Clinvar annotations of interest
        - are limited to specific transcript consequences

    Note, the test is run for each different consequence category
    but the 'keep' flag is used in every case

    Parameters
    ----------
    mt :
    acmg_data : ACMG data

    Returns
    -------
    same MT, with categories flagged where the variant types are appropriate
    """

    specific_types = {
        key: value for key, value in acmg_data.items() if value.get('specific_type')
    }
    types = {gene['specific_type'] for gene in specific_types.values()}
    for var_type in types:
        assert (
            var_type in VARIANT_CONSEQUENCES
        ), f'this script doesn\'t tolerate {var_type}'  # so fix it
        consequences = VARIANT_CONSEQUENCES[var_type]
        genes_for_this_consequence = hl.set(
            [key for key, value in specific_types if var_type in value['specific_type']]
        )
        mt = mt.annotate_rows(
            keep=hl.if_else(
                (genes_for_this_consequence.contains(mt.geneIds))
                & (
                    hl.len(
                        mt.vep.transcript_consequences.filter(
                            lambda x: hl.len(
                                consequences.intersection(hl.set(x.consequence_terms))
                            )
                            > 0
                        )
                    )
                    > 0
                )
                & (mt.clinvar == ONE_INT),
                ONE_INT,
                mt.keep,
            )
        )
    return mt


def find_exact_variants(
    mt: hl.MatrixTable, acmg_data: dict[str, dict]
) -> hl.MatrixTable:
    """
    flags variants where:
        - gene is flagged by the ACMG specification
        - have Clinvar annotations of interest
        - there is the exact protein consequence flagged in ACMG data

    Parameters
    ----------
    mt : the whole matrix table
    acmg_data : ACMG data

    Returns
    -------
    same MT, with categories flagged where the exact variants were found
    """

    specific_variants = {
        key: value['specific_variant']
        for key, value in acmg_data.items()
        if value.get('specific_variant')
    }

    # filthy double layered lambda
    # standard: require a gene match
    # lambda logic 1, per 'p.' change:
    #     find all CSQ HGVSp containing the exact change: True/False
    # lambda logic 2: per variant:
    #     filter for HGVSp test==True, and retain any variants with > 0 CSQs
    # The reason for this is that the VEP HGVSp isn't trivial to parse:
    # the entries are "ENSG####.n;p.A123B", inc. versioned ENSG ID
    # so constructing an exact match would be brittle/impossible
    # Instead we collect all p. for a gene, and use lambda + hl.any
    # This means in one round of evaluation we can find any loose matches
    for gene, variants in specific_variants.items():
        variant_set = hl.set(variants)
        # pylint: disable=W0108
        mt = mt.annotate_rows(
            keep=hl.if_else(
                (mt.geneIds == gene)
                & (
                    hl.len(
                        mt.vep.transcript_consequences.filter(
                            lambda tc_con: variant_set.any(
                                lambda change: tc_con.hgvsp.contains(change)
                            )
                        )
                    )
                    > 0
                ),
                ONE_INT,
                mt.keep,
            )
        )

    return mt


def filter_mt_to_keep(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    reduces the MT to relevant rows only
    Parameters
    ----------
    mt : all variants in the MT

    Returns
    -------
    the same MT, filtered to all rows where 'keep' is set to 1
    keep is set to 1 when any individual category is approved
    """
    return mt.filter_rows(mt.keep == ONE_INT)


def main(mt_path: str, acmg_path: str, output: str):
    """
    main method
    """

    # get the relevant acmg data in full
    acmg_data = read_json_from_path(acmg_path)

    # read the MT & do filtering
    mt = read_and_filter_mt(
        mt_path=mt_path,
        filter_af=0.01,
        green_genes=hl.set(acmg_data.keys()),
    )

    # assign a clinvar boolean flag once to simplify logic
    mt = add_top_level_clinvar_flag(mt)
    mt = checkpoint_and_repartition(
        mt,
        checkpoint_root=output_path('additional_findings.mt', category='tmp'),
        extra_logging='after applying filters and splitting genes',
    )

    # NOW FOR THE FUN PART!
    mt = find_standard_interesting(mt=mt, acmg_data=acmg_data)
    mt = find_specific_types(mt=mt, acmg_data=acmg_data)
    mt = find_exact_variants(mt=mt, acmg_data=acmg_data)
    mt = filter_mt_to_keep(mt)

    # all remaining rows get an AIP-compatible label
    mt = mt.annotate_rows(info=mt.info.annotate(categorybooleanaf=ONE_INT))

    # write the resulting content to a file/folder
    write_matrix_to_vcf(mt, file_name=output)


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--mt', help='the MT to use as input')
    parser.add_argument('--acmg', help='the parsed ACMG data')
    parser.add_argument('-o', help='VCF output name')
    args = parser.parse_args()
    main(mt_path=args.mt, acmg_path=args.acmg, output=args.o)
