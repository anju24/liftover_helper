import collections
import logging
import sys
from copy import copy
from io import StringIO

import vcf
from vcf.parser import _Info as VcfInfo

contig_spec = collections.namedtuple('Contig', 'id,length')
calldata_spec = collections.namedtuple('CallData', 'GT')
CHROMS = [str(i) for i in range(1, 23)] + ['X', 'Y']

# the following is an minimum vcf needed to generate new ones
# If there are no records in a vcf file, we use this template to add in mismatch sites
EXAMPLE_VCF = [
    "##fileformat=VCFv4.2",
    "##reference=file:///resources/b37/GRCh38.p12.fa",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##contig=<ID=chr2>",
    "##contig=<ID=chr6>",
    "##contig=<ID=chr15>",
    "##contig=<ID=chr19>",
    "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1",
    "chr2  21012603     .      C     T     100    PASS  .    GT  1/1",
]

def convert_hg19_vcf_to_grch37_vcf(input_vcf_file, output_vcf_file):
    """
    Liftover lifts from GRCh38 -> hg19. hg19 contigs have 'chr' prefix,
    that needs to be removed for GRCH37. Also filtering to only acceptable chrs
    autosomes, X and Y. No alternate contigs are considered
    :param input_vcf_file: hg19 vcf file liftedover from GRCh38
    :param output_vcf_file: GRCh37 contigs vcf file
    """
    chroms_to_keep = ['chr' + chrom for chrom in CHROMS]
    grch37_contigs = collections.OrderedDict()
    with open(input_vcf_file) as in_fp, open(output_vcf_file, 'w') as out_fp:
        reader = vcf.Reader(in_fp)
        # update contigs to remove 'chr' prefix
        for contig in reader.contigs:
            if contig in chroms_to_keep:
                value = contig_spec(reader.contigs[contig].id.replace('chr', ''), reader.contigs[contig].length)
                grch37_contigs.update({contig.replace('chr', ''): value})

        reader.contigs = grch37_contigs
        writer = vcf.Writer(out_fp, reader, lineterminator='\n')
        for record in reader:
            if record.CHROM in chroms_to_keep:
                # update chrom in the records
                record.CHROM = record.CHROM.replace('chr', '')
                writer.write_record(record)



MISMATCH_SITES = {
    'chr2:21012603': {
        '38_coordinates': {
            'chrom': 'chr2',
            'start': 21012603,
            'end': 21012604,
            'base': 'C',
        },
        '37_coordinates': {
            'chrom': 'chr2',
            'start': 21235475,
            'end': 21235476,
            'base': 'T',
        },
    },
    'chr6:7563750': {
        '38_coordinates': {
            'chrom': 'chr6',
            'start': 7563750,
            'end': 7563751,
            'base': 'G',
        },
        '37_coordinates': {
            'chrom': 'chr6',
            'start': 7563983,
            'end': 7563984,
            'base': 'T',
        },
    },
    'chr15:48515440': {
        '38_coordinates': {
            'chrom': 'chr15',
            'start': 48515440,
            'end': 48515441,
            'base': 'T',
        },
        '37_coordinates': {
            'chrom': 'chr15',
            'start': 48807637,
            'end': 48807638,
            'base': 'C',
        },
    },
    'chr19:55154216': {
        '38_coordinates': {
            'chrom': 'chr19',
            'start': 55154216,
            'end': 55154217,
            'base': 'C',
        },
        '37_coordinates': {
            'chrom': 'chr19',
            'start': 55665584,
            'end': 55665585,
            'base': 'A',
        },
    },
}


def find_overlapping_mismatch_site(record):
    for key, site in MISMATCH_SITES.items():
        if record.CHROM == site['38_coordinates']['chrom']:
            if record.POS <= site['38_coordinates']['start'] <= record.end:
                return key


def update_grch38_ref_to_grch37_for_record_if_needed(record, mismatched_site_key=None):
    """
    If record overlaps mismatched sites, update ref, alt and genotype accordingly
    Some assumptions:
    - records have only one sample (true for color data, one file per sample)
    - will not try to fix record, if genotype is malformed or missing
    - will not try to fix record, if more than 2 alleles
    """
    if not mismatched_site_key:
        mismatched_site_key = find_overlapping_mismatch_site(record)
    if not mismatched_site_key:
        return record
    mismatched_site = MISMATCH_SITES[mismatched_site_key]

    # get genotype
    # code only for one sample per vcf file. This is how color files are
    assert (len(record.samples) == 1)
    sample_call = record.samples[0]
    if 'GT' not in record.FORMAT:
        # missing genotype issue not handled
        raise ValueError(f'GT not in record {record}')
    gt_indices = sample_call.gt_alleles
    expected_gts = set(str(i) for i in range(len(record.ALT) + 1))
    if not set(gt_indices).issubset(expected_gts):
        # unknown genotype issue not handled
        raise ValueError(f'Unknown genotypes for record {record}, genotype {gt_indices}')

    ref = list(record.REF)
    pos = mismatched_site['38_coordinates']['start'] - record.POS

    if mismatched_site['38_coordinates']['base'] == ref[pos]:
        ref[pos] = mismatched_site['37_coordinates']['base']
    updated_ref = ''.join(ref)
    orig_allele_seqs = []
    # alleles from original record
    for allele in sample_call.gt_alleles:
        # only add reference, if it is in genotype, het variant
        if allele == '0':
            orig_allele_seqs.append(record.REF)
        else:
            orig_allele_seqs.append(record.ALT[int(allele) - 1])
    updated_alts = []
    updated_gt = []
    for allele_seq in orig_allele_seqs:
        # original allele matches updated ref, genotype needs '0' added
        if allele_seq == updated_ref:
            updated_gt.append(0)
        else:
            if allele_seq not in updated_alts:
                updated_alts.append(allele_seq)
            updated_gt.append(updated_alts.index(allele_seq) + 1)
    updated_gt = '/'.join(map(str, sorted(updated_gt)))

    if len(updated_alts) > 2:
        # case we will not try to handle, write original record as is
        raise ValueError(f'Updated record has more than 2 alts {updated_alts}, original record {record}')

    if len(updated_alts) == 0:
        # return None if the variant was not a variant in 37
        record = None
    else:
        # add genotype data to record. Loses all other fields
        record.FORMAT = 'GT'
        record.samples[0].data = calldata_spec(updated_gt)

        record.REF = updated_ref
        record.ALT = []
        record.add_info('preprocessed')
        for alt in updated_alts:
            record.ALT.append(vcf.model._Substitution(alt))
    return record


def _create_fake_record():
    """
    If the input vcf file is empty, use a fake record to use as template to create mismatch records
    """
    with StringIO('\n'.join(EXAMPLE_VCF)) as fh:
        records = list(vcf.Reader(fh))
        return records[0]


def convert_grch38_ref_mismatch_sites_to_grch37(input_vcf_file, output_vcf_basename):
    """
    For ACMG59 reportable range there are 4 sites that have
    reference mismatch between GRCh37 and GRCh38
    All ref and alts in variants overlapping these sites
    will need to be updated to 37 reference
    output file will contain variants overlapping mismatch sites
    and all other variants with original record
    """
    logger = logging.getLogger(__name__)
    output_vcf_file = f'{output_vcf_basename}.vcf'
    reader = vcf.Reader(filename=input_vcf_file)
    records = list(reader)
    if len(records) == 0:
        record = _create_fake_record()
        logger.warning(f'VCF file {input_vcf_file} is empty. We will be creating the mismatch sites using a template')

    mismatched_site_overlap = {}
    for record in records:
        mismatched_site_key = find_overlapping_mismatch_site(record)
        if mismatched_site_key:
            mismatched_site_overlap[mismatched_site_key] = True
            try:
                update_grch38_ref_to_grch37_for_record_if_needed(record, mismatched_site_key)
            except ValueError as e:
                logger.info(f'Record {record.CHROM}:{record.POS} with mismatch site {mismatched_site_key} encountered error {e}')

    reader.infos['PREPROCESSED'] = VcfInfo(
        'PREPROCESSED',
        0,
        'Flag',
        'The record was pre-processed. Added when a record needed to be changed for liftover',
        '',
        '',
    )
    # if there are no overlapping variants in mismatched sites,
    # create a homozygous variant matching 37 as ref and 38 as alt
    for key, site in MISMATCH_SITES.items():
        if key not in mismatched_site_overlap.keys():
            # TODO: separate out creation of a record
            mismatch_record = copy(record)
            mismatch_record.ID = '.'
            mismatch_record.QUAL = 100
            mismatch_record.FILTER = []
            mismatch_record.FORMAT = 'GT'
            mismatch_record.samples[0].data = calldata_spec('1/1')
            mismatch_record.add_info('preprocessed')
            mismatch_record.CHROM = site['38_coordinates']['chrom']
            mismatch_record.POS = site['38_coordinates']['start']
            mismatch_record.REF = site['37_coordinates']['base']
            mismatch_record.ALT = [
                vcf.model._Substitution(site['38_coordinates']['base'])
            ]
            records.append(mismatch_record)

    contig_order = {c: i for i, c in enumerate(reader.contigs)}

    def sort_key(record):
        """
        Sorts records by (CHROM,POS,REF).
        If contigs are specified in the VCF file and record CHROM matches a contig,
        contig order is maintained.
        Any unmatched CHROMs will throw an error
        """
        if record.CHROM not in contig_order:
            raise ValueError(
                f'Unexpected chrom {record.CHROM} found. Expected one of {contig_order.keys()}'
            )
        return (contig_order[record.CHROM], record.POS, record.REF)

    records.sort(key=sort_key)

    with open(output_vcf_file, 'w') as out_fp:
        writer = vcf.Writer(out_fp, reader, lineterminator='\n')
        for record in records:
            writer.write_record(record)



if __name__ == '__main__':
    function = getattr(sys.modules[__name__], sys.argv[1])
    function(sys.argv[2], sys.argv[3])
