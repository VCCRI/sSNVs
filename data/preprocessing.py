import hail as hl

from misc import *

def preprocessing(
    exomes_ht,
    context_ht,
    mutation_ht,
    sex_split,
    with_gerp=False,
    with_mutability=True,
    standard_annotation=True,
):
    """Filtering steps for all synonymous variants.

    exomes_ht -- WES variants in ht format
    context_ht -- context Hail table
    mutation_ht -- mutation rates
    sex_split -- how many males, how many females
    with_gerp -- whether to annotate with GERP scores
    with_mutability -- whether to annotate with mutability
    standard_annotation -- whether to apply additional standard annotation
    """

    exomes = hl.read_table(exomes_ht)

    # Allele number (AN) adjustment.
    exomes = exomes.filter(get_an_adj_criteria(exomes, sex_split))

    # Filter the table so that only synonymous variants
    # with AF>0 and filter PASS are retained.
    # The first condition is necessary because in gnomAD variants
    # that were excluded from the analysis through QC have AF=0.
    # The 'most_severe_consequence' condition removes all
    # pLoF variants, including variants in canonical splice sites.
    synonymous = exomes.filter(
        (exomes.freq[0].AF > 0)
        & (exomes.filters.length() == 0)
        & (exomes.vep.most_severe_consequence == "synonymous_variant")
    )

    # For each variant there is a list of
    # 'transcript_consequences'. The code below extracts information
    # about the first transcript where the variant is synonymous and
    # .biotype is 'protein_coding' (there may be several such
    # transcripts). Those mutations with 'synonymous_variant' as their
    # most severe consequence that don't have such a transcript (a
    # very small number) will become NAs.
    syn_coding = synonymous.annotate(
        transcript_consequences=synonymous.vep.transcript_consequences.find(
            lambda x: (x.consequence_terms == ["synonymous_variant"])
            & (x.biotype == "protein_coding")
        )
    )

    if standard_annotation:
        syn_coding = syn_coding.annotate(
            ref_codon=syn_coding.transcript_consequences.codons.split("/")[0],
            alt_codon=syn_coding.transcript_consequences.codons.split("/")[1],
        )

    # Methylation and other context data.
    context = hl.read_table(context_ht)
    context = context[syn_coding.key]
    # The 2020 version of MAPS uses methylation.
    # Function 'prepare_ht' annotates the input table with methylation level,
    # coverage (optional), CpG/Non-CpG info, context for mutability
    # (ref allele in the middle plus two bases to the left and to the right)
    # and other, less important information.
    # For example, a variant that has changed '|...|..t|Cta|...|'
    # to '|...|..t|Tta|...|' will have 'ref_codon'="Cta",
    # 'alt_codon'="Tta", 'ref'="C", 'alt'="T" and 'context'="TCT".
    if with_gerp:
        syn_coding = prepare_ht(
            syn_coding.annotate(
                context=context.context,
                methylation=context.methylation,
                gerp=context.gerp,
            ),
            trimer=True,
            annotate_coverage=False,
        )
    else:
        syn_coding = prepare_ht(
            syn_coding.annotate(
                context=context.context, methylation=context.methylation
            ),
            trimer=True,
            annotate_coverage=False,
        )

    # Annotate with mutability (i.e. mutation rate per base pair per generation).
    # The 'mutation_ht' table contains 104 entries: 32 context sequences
    # (4^3 divided by 2 because of complementarity), times
    # 3 possible mutations of the middle nucleotide at methylation level 0,
    # plus for context sequences ACG, CCG, GCG and TCG
    # (which all have C adjacent to G) there are 8 additional possible
    # CpG transitions at methylation levels 1 and 2.
    # However, there will be only 100 different mutability values
    # because some context transitions (e.g. AAT->ACT) are always
    # nonsynonymous mutations.
    if with_mutability:
        mutation_rates = hl.read_table(mutation_ht)
        syn_coding = syn_coding.annotate(
            mu=mutation_rates[
                hl.struct(
                    context=syn_coding.context,
                    ref=syn_coding.ref,
                    alt=syn_coding.alt,
                    methylation_level=syn_coding.methylation_level,
                )
            ].mu_snp
        )

    return syn_coding


def annotate_quartiles(ht, variable):
    quartiles = ht.aggregate(
        hl.agg.approx_quantiles(ht[variable], [0, 0.25, 0.5, 0.75, 1])
    )
    ht = ht.annotate(
        quartile=hl.case()
        .when(ht[variable] <= quartiles[1], 1)
        .when(
            (ht[variable] > quartiles[1]) & (ht[variable] <= quartiles[2]),
            2,
        )
        .when(
            (ht[variable] > quartiles[2]) & (ht[variable] <= quartiles[3]),
            3,
        )
        .when(ht[variable] > quartiles[3], 4)
        .or_missing()
    )
    return ht
