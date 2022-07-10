"""Testing facilities for crops.elements.intervals"""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

from crops.elements import sequences as ces

import unittest

_SEQUENCE_1 = """MLRIPVTRALIGLSKSPKGCVRTTATAASNLIEVFVDGQSVMVEPGTTVLQACEKVGMQIPRFCYHERLSVAGNCRMCL
VEIEKAPKVVAACAMPVMKGWNILTNSEKSKKAREGVMEFLLANHPLDCPICDQGGECDLQDQSMMFGSDRSRFLEGKR
AVEDKNIGPLVKTIMTRCIQCTRCIRFASEIAGVDDLGTTGRGNDMQVGTYIEKMFMSELSGNIIDICPVGALTSKPYA
FTARPWETRKTESIDVMDAVGSNIVVSTRTGEVMRILPRMHEDINEEWISDKTRFAYDGLKRQRLTQPMIRNEKGLLTY
TTWEDALSRVAGMLQSFQGNDVAAIAGGLVDAEALVALKDLLNRVDSDSLCTEEVFPTAGAGTDLRSNYLLNTTIAGVE
EADVILLVGTNPRFEAPLFNARIRKSWLHNDLKVALIGSPVDLTYRYDHLGDSPKILQDIASGNHPFSQILKEAKKPMV
VLGSSALQRSDGTAILAAVSNIAQNIRLSSGVTGDWKVMNILHRIASQVAALDLGYKPGVEAIRKNPPKVLFLLGADGG
CITRQDLPKDCFIIYQGHHGDVGAPMADVILPGAAYTEKSATYVNTEGRAQQTKVAVTPPGLAREDWKIIRALSEIAGM
TLPYDTLDQVRSRLEEVSPNLVRYDDVEGANYFQQANELSKLVNQQLLADPLVPPQLTIKDFYMTDSISRASQTMAKCV
KAVTEGIQAVEEPSIC
"""
_SEQUENCE_2 = """PKLVLVRHGQSEWNEKNLFTGWVDVKLSAKGQQEAARAGELLKEKKVYPDVLYTSKLSRAIQTANIALEKADRLWIPVN
RSWRLNERHYGDLQGKDKAETLKKFGEEKFNTYRRSFDVPPPPIDASSPFSQKGDERYKYVDPNVLPETESLALVIDRL
LPYWQDVIAKDLLSGKTVMIAAHGNSLRGLVKHLEGISDADIAKLNIPTGIPLVFELDENLKPSKPSYYLDPEAAAAGA
AAVANQGKK
"""
_SEQUENCE_3 = "GATACTCAGATAG"
_SEQUENCE_4 = "CUAUCUGAGUAUC"
_SEQUENCE_5 = "GATACTNAGATAG"
_SEQUENCE_6 = "*AT*CTNA*ATAG"
_SEQUENCE_7 = "-AT-CTNA-ATAG"
_SEQUENCE_8 = "G-T-CTN-GATAG"

_CROPMAP_1 = {1: None, 2: 1, 3: 2, 4: None, 5: 3, 6: 4, 7: 5,
              8: 6, 9: None, 10: 7, 11: 8, 12: 9, 13: 10}
_CROPMAP_2 = {1: 2, 2: 3, 3: 5, 4: 6, 5: 7,
              6: 8, 7: 10, 8: 11, 9: 12, 10: 13}

_HEADER_1 = ">1IXY_1|Chains A[auth C], C[auth D]|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|"

class TestCropsSequences(unittest.TestCase):
    def test_guess_type_1(self):
        expected_type = "Protein"

        obtained_type = ces.guess_type(_SEQUENCE_1)

        self.assertEqual(obtained_type, expected_type)

    def test_guess_type_2(self):
        expected_type = "Protein"

        obtained_type = ces.guess_type(_SEQUENCE_2)

        self.assertEqual(obtained_type, expected_type)

    def test_guess_type_3(self):
        expected_type = "DNA"

        obtained_type = ces.guess_type(_SEQUENCE_3)

        self.assertEqual(obtained_type, expected_type)

    def test_guess_type_4(self):
        expected_type = "RNA"

        obtained_type = ces.guess_type(_SEQUENCE_4)

        self.assertEqual(obtained_type, expected_type)

    def test_sequence_addseq_1(self):
        expected_seq = _SEQUENCE_3

        seq = ces.sequence(seqid=1, oligomer_id='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.addseq(newid="alternative", newseq=_SEQUENCE_3)

        self.assertTrue("alternative" in seq.seqs)
        self.assertEqual(seq.seqs["alternative"], expected_seq)

    def test_sequence_delseq_1(self):
        expected_seq = _SEQUENCE_5

        seq = ces.sequence(seqid=1, oligomer_id='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.addseq(newid="alternative", newseq=_SEQUENCE_3)

        seq.delseq(delid="alternative")

        self.assertFalse("alternative" in seq.seqs)
        self.assertTrue("mainseq" in seq.seqs)
        self.assertEqual(seq.seqs["mainseq"], expected_seq)

    def test_sequence_delseq_2(self):
        expected_seq = ""

        seq = ces.sequence(seqid=1, oligomer_id='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.addseq(newid="alternative", newseq=_SEQUENCE_3)

        seq.delseq(wipeall=True)

        self.assertFalse("alternative" in seq.seqs)
        self.assertTrue("mainseq" in seq.seqs)
        self.assertEqual(seq.seqs["mainseq"], expected_seq)

    def test_sequence_mainseq_1(self):
        expected_seq = _SEQUENCE_5

        seq = ces.sequence(seqid=1, oligomer_id='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        returned_seq = seq.mainseq()

        self.assertEqual(expected_seq, returned_seq)

    def test_sequence_mainseq_2(self):
        expected_seq = _SEQUENCE_3

        seq = ces.sequence(seqid=1, oligomer_id='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.mainseq(add=_SEQUENCE_3)
        returned_seq = seq.mainseq()

        self.assertEqual(expected_seq, returned_seq)

    def test_sequence_guess_biotype_1(self):
        expected_type = ces.guess_type(_SEQUENCE_3)

        seq = ces.sequence(seqid=1, oligomer_id='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_3)

        returned_type = seq.guess_biotype()

        self.assertEqual(expected_type, returned_type)

    def test_sequence_dump_1(self):
        expected_output = """>crops|1IXY_1|Chains C,D|Source: RCSB PDB|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|
        GATACTNAGATAG
        """

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        returned_output = seq.dump(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_dump_2(self):
        expected_output = """>crops|2IXY_2|Chains A,C|Source: MADEUP|No info
        GATACTNAGATAG
        """

        seq = ces.sequence(seqid=2, oligomer_id='2IXY', chains={'A', 'C'},
                           seq=_SEQUENCE_5, header=_HEADER_1, extrainfo='No info',
                           source='MADEUP')

        returned_output = seq.dump(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_dumpmap_1(self):
        expected_output = """>crops|1IXY_1|Chains C,D|Source: RCSB PDB|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|
        1  0
        2  1
        3  2
        4  0
        5  3
        6  4
        7  5
        8  6
        9  0
        10  7
        11  8
        12  9
        13  10
        """

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].strip('*')

        returned_output = seq.dumpmap(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_dumpmap_2(self):
        expected_output = """>crops|2IXY_2|Chains A,C|Source: MADEUP|No info
        1  0
        2  1
        3  2
        4  0
        5  3
        6  4
        7  5
        8  6
        9  0
        10  7
        11  8
        12  9
        13  10
        """

        seq = ces.sequence(seqid=2, oligomer_id='2IXY', chains={'A', 'C'},
                           seq=_SEQUENCE_5, header=_HEADER_1, extrainfo='No info',
                           source='MADEUP')

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].strip('*')

        returned_output = seq.dumpmap(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_length_1(self):
        expected_length = 249

        seq = ces.sequence(seq=_SEQUENCE_2)

        obtained_length = seq.length()

        self.assertEqual(expected_length, obtained_length)

    def test_sequence_length_2(self):
        expected_length = 13

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        obtained_length = seq.length()

        self.assertEqual(expected_length, obtained_length)

    def test_sequence_full_length_1(self):
        expected_length = 10

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].strip('*')

        obtained_length = seq.full_length()

        self.assertEqual(expected_length, obtained_length)

    def test_sequence_ncrops_1(self):
        expected_ncrops = 3

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].strip('*')

        obtained_ncrops = seq.ncrops()

        self.assertEqual(expected_ncrops, obtained_ncrops)

    def test_sequence_ngaps_1(self):
        expected_ngaps = 3

        seq = ces.sequence(seq=_SEQUENCE_2)

        seq.seqs['gapseq'] = [_SEQUENCE_7]

        obtained_ngaps = seq.ngaps()

        self.assertEqual(expected_ngaps, obtained_ngaps)

    def test_sequence_ngaps_2(self):
        expected_ngaps = [3, 3]

        seq = ces.sequence(seq=_SEQUENCE_2)

        seq.seqs['gapseq'] = [_SEQUENCE_7, _SEQUENCE_8]

        obtained_ngaps = seq.ngaps()

        self.assertEqual(expected_ngaps, obtained_ngaps)

    def test_sequence_update_cropsheader_1(self):
        expected_output = ">crops|1IXY_2|Chains C,D|Source: RCSB PDB|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|"

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq.name = '2'

        seq.update_cropsheader()

        returned_output = seq.crops_header

        self.assertEqual(expected_output, returned_output)

    def test_sequence_cropinfo_1(self):
        expected_info = '#Residues cropped: 3 (2 not from terminals) ; % cropped: 23.08 (15.38 not from terminals)'

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].strip('*')

        obtained_info = seq.ncrops()

        self.assertEqual(expected_info, obtained_info)
