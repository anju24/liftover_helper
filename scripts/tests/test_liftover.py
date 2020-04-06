import collections
import os
import tempfile

import vcf

from datasets.vcf import liftover
from util import testutil
from util.files import open
import unittest


class RecordTest(unittest.TestCase):
    @testutil.file(
        'mem:/hg19.vcf', """\
##fileformat=VCFv4.2
##INFO=<ID=ReverseComplementedAlleles,Number=0,Type=Flag,Description="The REF \
and the ALT alleles have been reverse complemented in liftover since the \
mapping from the previous reference to the current one was on the negative strand.">
##INFO=<ID=SwappedAlleles,Number=0,Type=Flag,Description="The REF and the ALT alleles have \
been swapped in liftover due to changes in the reference. It is possible that not all INFO \
annotations reflect this swap, and in the genotypes, only the GT, PL, and AD fields have \
been modified. You should check the TAGS_TO_REVERSE parameter that was used during the \
LiftOver to be sure.">
##contig=<ID=chr1,length=249250621>
##contig=<ID=chr2,length=243199373>
##contig=<ID=chr3,length=198022430>
##contig=<ID=chr6_ssto_hap7,length=4928567>
##contig=<ID=chr6_mcf_hap5,length=4833398>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1
chr1    97915604        .       A       G       100     PASS    .       GT      0/1
chr2    97915605        .       T       A       100     PASS    .       GT      0/1
"""
    )
    def test_convert_hg19_vcf_to_grch37_vcf(self):
        temp_dir = tempfile.gettempdir()
        output_file = os.path.join(temp_dir, 'output_grch37.vcf')
        liftover.convert_hg19_vcf_to_grch37_vcf('mem:/hg19.vcf', output_file)

        with open(output_file) as fh:
            records = vcf.Reader(fh)
            expected_contigs = collections.OrderedDict([('1', liftover.contig_spec('1', 249250621)),
                                                        ('2', liftover.contig_spec('2', 243199373)),
                                                        ('3', liftover.contig_spec('3', 198022430))])
            self.assertEqual(records.contigs, expected_contigs)
            record = next(records)
            self.assertEqual(record.CHROM, '1')
            self.assertEqual(record.POS, 97915604)
            record = next(records)
            self.assertEqual(record.CHROM, '2')
            self.assertEqual(record.POS, 97915605)

    @testutil.file(
        'mem:/grch38.vcf', """\
##fileformat=VCFv4.2
##reference=file:///resources/b37/GRCh38.p12.fa
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1
chr2    21012602        .       ACATG   A       100     PASS    .       GT      0/1
chr2    21012603        .       C       A       100     PASS    .       GT      0/1
chr2    21012603        .       CATG   C       100     PASS    .       GT      0/1
chr2    21012603        .       C       CAAT    100     PASS    .       GT      0/1
chr2    21012603        .       C       T,CAAT  100     PASS    .       GT      1/2
chr2    21012604        .       ATG   A       100     PASS    .       GT      0/1
chr2    21012622        .       C       T       100     PASS    .       GT      0/1
"""
    )
    def test_record_overlaps_mismatch(self):
        records = list(vcf.Reader(open('mem:/grch38.vcf')))
        expected_result = {
            '38_coordinates': {
                'chrom': 'chr2',
                'start': 21012603,
                'end': 21012604,
                'base': 'C'
            },
            '37_coordinates': {
                'chrom': 'chr2',
                'start': 21235475,
                'end': 21235476,
                'base': 'T'
            }
        }
        for i in range(0, 5):
            self.assertEqual(liftover.record_overlaps_mismatch_sites(records[i]), expected_result)
        for i in range(5, 7):
            self.assertFalse(liftover.record_overlaps_mismatch_sites(records[i]))

    @testutil.file(
        'mem:/grch38.vcf', """\
##fileformat=VCFv4.2
##reference=file:///resources/b37/GRCh38.p12.fa
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1
chr2    21012602        .       ACATG   A       100     PASS    .       GT      0/1
chr2    21012602        .       ACATG   A       100     PASS    .       GT      1/1
chr2    21012602        .       AC   A       100     PASS    .       GT      0/1
chr2    21012602        .       AC   A       100     PASS    .       GT      1/1
chr2    21012603        .       C       A       100     PASS    .       GT      0/1
chr2    21012603        .       C       A       100     PASS    .       GT      1/1
chr2    21012603        .       C       T       100     PASS    .       GT      0/1
chr2    21012603        .       C       T       100     PASS    .       GT      1/1 # most probable
chr2    21012603        .       C       T       100     PASS    .       GT      ./.
chr2    21012603        .       C       T       100     PASS    .       GT      1/2 # shouldnt happen
chr2    21012603        .       CATG   C       100     PASS    .       GT      0/1
chr2    21012603        .       CATG   C       100     PASS    .       GT      1/1
chr2    21012603        .       C       CAAT    100     PASS    .       GT      0/1
chr2    21012603        .       C       T,CAAT  100     PASS    .       GT      1/2
chr2    21012604        .       ATG   A       100     PASS    .       GT      0/1
chr2    21012622        .       C       T       100     PASS    .       GT      0/1
"""
    )
    def test_update_record(self):
        records = list(vcf.Reader(open('mem:/grch38.vcf')))
        expected_results = [
            ['ATATG', 'ACATG,A', '1/2'],
            ['ATATG', 'A', '1/1'],
            ['AT', 'AC,A', '1/2'],
            ['AT', 'A', '1/1'],
            ['T', 'C,A', '1/2'],
            ['T', 'A', '1/1'],
            ['T', 'C', '0/1'],
            None,
            ['C', 'T', './.', True],
            ['C', 'T', '1/2', True],
            ['TATG', 'CATG,C', '1/2'],
            ['TATG', 'C', '1/1'],
            ['T', 'C,CAAT', '1/2'],
            ['T', 'CAAT', '0/1'],
            ['ATG', 'A', '0/1'],
            ['C', 'T', '0/1'],
        ]
        for i in range(len(records)):
            expected_record = expected_results[i]
            if expected_record and len(expected_record) == 4:
                self.assertRaises(ValueError, liftover.update_grch38_ref_to_grch37_for_record_if_needed, records[i])
                continue

            observed_record = liftover.update_grch38_ref_to_grch37_for_record_if_needed(records[i])
            if expected_record is None:
                self.assertIsNone(observed_record)
            else:
                self.assertEqual(observed_record.REF, expected_record[0])
                observed_alts = ','.join(map(str, observed_record.ALT))
                self.assertEqual(observed_alts, expected_record[1])
                self.assertEqual(observed_record.samples[0].data.GT, expected_record[2])

        # TODO: tests for is_anchor_base

if __name__ == '__main__':
    unittest.main()
