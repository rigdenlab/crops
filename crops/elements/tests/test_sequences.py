"""Testing facilities for crops.elements.intervals"""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

from crops.elements import sequences as ces

import unittest
import os

_SEQUENCE_1 = ("MLRIPVTRALIGLSKSPKGCVRTTATAASNLIEVFVDGQSVMVEPGTTVLQACEKVGMQIP" +
               "RFCYHERLSVAGNCRMCLVEIEKAPKVVAACAMPVMKGWNILTNSEKSKKAREGVMEFLLA" +
               "NHPLDCPICDQGGECDLQDQSMMFGSDRSRFLEGKRAVEDKNIGPLVKTIMTRCIQCTRCI" +
               "RFASEIAGVDDLGTTGRGNDMQVGTYIEKMFMSELSGNIIDICPVGALTSKPYAFTARPWE" +
               "TRKTESIDVMDAVGSNIVVSTRTGEVMRILPRMHEDINEEWISDKTRFAYDGLKRQRLTQP" +
               "MIRNEKGLLTYTTWEDALSRVAGMLQSFQGNDVAAIAGGLVDAEALVALKDLLNRVDSDSL" +
               "CTEEVFPTAGAGTDLRSNYLLNTTIAGVEEADVILLVGTNPRFEAPLFNARIRKSWLHNDL" +
               "KVALIGSPVDLTYRYDHLGDSPKILQDIASGNHPFSQILKEAKKPMVVLGSSALQRSDGTA" +
               "ILAAVSNIAQNIRLSSGVTGDWKVMNILHRIASQVAALDLGYKPGVEAIRKNPPKVLFLLG" +
               "ADGGCITRQDLPKDCFIIYQGHHGDVGAPMADVILPGAAYTEKSATYVNTEGRAQQTKVAV" +
               "TPPGLAREDWKIIRALSEIAGMTLPYDTLDQVRSRLEEVSPNLVRYDDVEGANYFQQANEL" +
               "SKLVNQQLLADPLVPPQLTIKDFYMTDSISRASQTMAKCVKAVTEGIQAVEEPSIC")

_SEQUENCE_2 = ("PKLVLVRHGQSEWNEKNLFTGWVDVKLSAKGQQEAARAGELLKEKKVYPDVLYTSKLSRAI" +
               "QTANIALEKADRLWIPVNRSWRLNERHYGDLQGKDKAETLKKFGEEKFNTYRRSFDVPPPP" +
               "IDASSPFSQKGDERYKYVDPNVLPETESLALVIDRLLPYWQDVIAKDLLSGKTVMIAAHGN" +
               "SLRGLVKHLEGISDADIAKLNIPTGIPLVFELDENLKPSKPSYYLDPEAAAAGAAAVANQGKK")

_SEQUENCE_3 = "GATACTCAGATAG"
_SEQUENCE_4 = "CUAUCUGAGUAUC"
_SEQUENCE_5 = "GATACTNAGATAG"
_SEQUENCE_5_2 = "CTATCTGAGTATC"
_SEQUENCE_6 = "+AT+CTNA+ATAG"
_SEQUENCE_6_2 = "+TA+CTGA+TATC"
_SEQUENCE_7 = "-AT-CTNA-ATAG"
_SEQUENCE_8 = "G-T-CTN-GATAG"

_CROPMAP_1 = {1: None, 2: 1, 3: 2, 4: None, 5: 3, 6: 4, 7: 5,
              8: 6, 9: None, 10: 7, 11: 8, 12: 9, 13: 10}
_CROPMAP_2 = {1: 2, 2: 3, 3: 5, 4: 6, 5: 7,
              6: 8, 7: 10, 8: 11, 9: 12, 10: 13}

_HEADER_1 = ">1IXY_1|Chains A[auth C], C[auth D]|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|"
_HEADER_2 = ">1IXY_2|Chains B[auth E], D[auth F]|5'-D(*CP*TP*AP*TP*CP*TP*GP*AP*GP*TP*AP*TP*C)-3'|"


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

        seq = ces.sequence(seqid=1, oligomer='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.addseq(newid="alternative", newseq=_SEQUENCE_3)

        self.assertTrue("alternative" in seq.seqs)
        self.assertEqual(seq.seqs["alternative"], expected_seq)

    def test_sequence_delseq_1(self):
        expected_seq = _SEQUENCE_5

        seq = ces.sequence(seqid=1, oligomer='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.addseq(newid="alternative", newseq=_SEQUENCE_3)

        seq.delseq(delid="alternative")

        self.assertFalse("alternative" in seq.seqs)
        self.assertTrue("mainseq" in seq.seqs)
        self.assertEqual(seq.seqs["mainseq"], expected_seq)

    def test_sequence_delseq_2(self):
        expected_seq = ""

        seq = ces.sequence(seqid=1, oligomer='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.addseq(newid="alternative", newseq=_SEQUENCE_3)

        seq.delseq(wipeall=True)

        self.assertFalse("alternative" in seq.seqs)
        self.assertTrue("mainseq" in seq.seqs)
        self.assertEqual(seq.seqs["mainseq"], expected_seq)

    def test_sequence_mainseq_1(self):
        expected_seq = _SEQUENCE_5

        seq = ces.sequence(seqid=1, oligomer='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        returned_seq = seq.mainseq()

        self.assertEqual(expected_seq, returned_seq)

    def test_sequence_mainseq_2(self):
        expected_seq = _SEQUENCE_3

        seq = ces.sequence(seqid=1, oligomer='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_5)

        seq.mainseq(add=_SEQUENCE_3)
        returned_seq = seq.mainseq()

        self.assertEqual(expected_seq, returned_seq)

    def test_sequence_guess_biotype_1(self):
        expected_type = ces.guess_type(_SEQUENCE_3)

        seq = ces.sequence(seqid=1, oligomer='1IXY', chains={'C', 'D'},
                           seq=_SEQUENCE_3)

        returned_type = seq.guess_biotype()

        self.assertEqual(expected_type, returned_type)

    def test_sequence_dump_1(self):
        expected_output = (">crops|1IXY_1|Chains C,D|Source: RCSB PDB|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|" +
                           os.linesep + "GATACTNAGATAG" + os.linesep)
        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        returned_output = seq.dump(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_dump_2(self):
        expected_output = (">crops|2IXY_2|Chains A,C|Source: MADEUP|No info" +
                           os.linesep + "GATACTNAGATAG" + os.linesep)

        seq = ces.sequence(seqid=2, oligomer='2IXY', chains={'A', 'C'},
                           seq=_SEQUENCE_5, header=_HEADER_1, extrainfo='No info',
                           source='MADEUP')

        returned_output = seq.dump(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_dumpmap_1(self):
        expected_output = (">crops|1IXY_1|Chains C,D|Source: RCSB PDB|5'-D(*GP*AP*TP*AP*CP*TP*3DRP*AP*GP*AP*TP*AP*G)-3'|" +
                           os.linesep + "1  0" + os.linesep + "2  1" + os.linesep +
                           "3  2" + os.linesep + "4  0" + os.linesep + "5  3" +
                           os.linesep + "6  4" + os.linesep + "7  5" + os.linesep +
                           "8  6" + os.linesep + "9  0" + os.linesep + "10  7" +
                           os.linesep + "11  8" + os.linesep + "12  9" + os.linesep +
                           "13  10" + os.linesep)

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].replace("+", "")

        returned_output = seq.dumpmap(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_dumpmap_2(self):
        expected_output = (">crops|2IXY_2|Chains A,C|Source: MADEUP|No info" +
                           os.linesep + "1  0" + os.linesep + "2  1" + os.linesep +
                           "3  2" + os.linesep + "4  0" + os.linesep + "5  3" +
                           os.linesep + "6  4" + os.linesep + "7  5" + os.linesep +
                           "8  6" + os.linesep + "9  0" + os.linesep + "10  7" +
                           os.linesep + "11  8" + os.linesep + "12  9" + os.linesep +
                           "13  10" + os.linesep)
        seq = ces.sequence(seqid=2, oligomer='2IXY', chains={'A', 'C'},
                           seq=_SEQUENCE_5, header=_HEADER_1, extrainfo='No info',
                           source='MADEUP')

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].replace("+", "")

        returned_output = seq.dumpmap(out='string')

        self.assertEqual(expected_output, returned_output)

    def test_sequence_length_1(self):
        expected_length = 246

        seq = ces.sequence(seq=_SEQUENCE_2)

        obtained_length = seq.length()

        self.assertEqual(expected_length, obtained_length)

    def test_sequence_length_2(self):
        expected_length = 13

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        obtained_length = seq.length()

        self.assertEqual(expected_length, obtained_length)

    def test_sequence_full_length_1(self):
        expected_length = 13

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].replace("+", "")

        obtained_length = seq.full_length()

        self.assertEqual(expected_length, obtained_length)

    def test_sequence_ncrops_1(self):
        expected_ncrops = 3

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)

        seq.cropmap = _CROPMAP_1
        seq.cropbackmap = _CROPMAP_2
        seq.seqs['cropseq'] = _SEQUENCE_6
        seq.seqs['fullseq'] = seq.seqs['mainseq']
        seq.seqs['mainseq'] = seq.seqs['cropseq'].replace("+", "")

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
        seq.seqs['mainseq'] = seq.seqs['cropseq'].replace("+", "")

        obtained_info = seq.cropinfo()

        self.assertEqual(expected_info, obtained_info)

    def test_oligoseq_add_sequence_1(self):
        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_5_2, header=_HEADER_2)

        expected_keys = {'1', '2'}
        expected_seqs = {seq, seq_2}

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        obtained_keys = set(oseq.imer.keys())
        obtained_seqs = set(oseq.imer.values())

        self.assertEqual(expected_keys, obtained_keys)
        self.assertEqual(expected_seqs, obtained_seqs)

    def test_oligoseq_del_sequence_1(self):
        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_5_2, header=_HEADER_2)

        expected_keys = {'2'}
        expected_seqs = {seq_2}

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        oseq.del_sequence(1)

        obtained_keys = set(oseq.imer.keys())
        obtained_seqs = set(oseq.imer.values())

        self.assertEqual(expected_keys, obtained_keys)
        self.assertEqual(expected_seqs, obtained_seqs)

    def test_oligoseq_purge_1(self):
        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_5_2, header=_HEADER_2)

        expected_keys = set()
        expected_seqs = set()

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        oseq.purge()

        obtained_keys = set(oseq.imer.keys())
        obtained_seqs = set(oseq.imer.values())

        self.assertEqual(expected_keys, obtained_keys)
        self.assertEqual(expected_seqs, obtained_seqs)

    def test_oligoseq_set_cropmaps_1(self):
        expected_seqs = []
        expected_seqs.append(_SEQUENCE_5)
        expected_seqs.append(_SEQUENCE_6)
        expected_seqs.append(_SEQUENCE_6.replace("+", ""))
        expected_seqs.append(_SEQUENCE_5_2)
        expected_seqs.append(_SEQUENCE_6_2)
        expected_seqs.append(_SEQUENCE_6_2.replace("+", ""))

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_5_2, header=_HEADER_2)

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        cropmaps = {}
        for key in oseq.imer:
            cropmaps[key] = {'cropmap': _CROPMAP_1,
                             'cropbackmap': _CROPMAP_2}

        oseq.set_cropmaps(cropmaps)

        obtained_seqs = []
        obtained_seqs.append(oseq.imer['1'].seqs['fullseq'])
        obtained_seqs.append(oseq.imer['1'].seqs['cropseq'])
        obtained_seqs.append(oseq.imer['1'].seqs['mainseq'])
        obtained_seqs.append(oseq.imer['2'].seqs['fullseq'])
        obtained_seqs.append(oseq.imer['2'].seqs['cropseq'])
        obtained_seqs.append(oseq.imer['2'].seqs['mainseq'])

        self.assertDictEqual(oseq.imer['1'].cropmap, _CROPMAP_1)
        self.assertDictEqual(oseq.imer['1'].cropbackmap, _CROPMAP_2)
        self.assertDictEqual(oseq.imer['2'].cropmap, _CROPMAP_1)
        self.assertDictEqual(oseq.imer['2'].cropbackmap, _CROPMAP_2)
        for n in range(len(expected_seqs)):
            self.assertEqual(expected_seqs[n], obtained_seqs[n])

    def test_oligoseq_length_1(self):
        expected_length_1 = 13
        expected_length_2 = 10

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_6_2.replace("+", ""), header=_HEADER_2)

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        obtained_length_1 = oseq.length(1)
        obtained_length_2 = oseq.length(2)

        self.assertEqual(expected_length_1, obtained_length_1)
        self.assertEqual(expected_length_2, obtained_length_2)

    def test_oligoseq_nchains_1(self):
        expected_value = 4

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_6_2.replace("+", ""), header=_HEADER_2)

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        obtained_value = oseq.nchains()

        self.assertEqual(expected_value, obtained_value)

    def test_oligoseq_nseqs_1(self):
        expected_value = 2

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_6_2.replace("+", ""), header=_HEADER_2)

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        obtained_value = oseq.nseqs()

        self.assertEqual(expected_value, obtained_value)

    def test_oligoseq_chainlist_1(self):
        expected_chains = {'C', 'D', 'E', 'F'}

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_6_2.replace("+", ""), header=_HEADER_2)

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        obtained_chains = oseq.chainlist()

        self.assertEqual(expected_chains, obtained_chains)

    def test_oligoseq_whatseq_1(self):
        chains = [{'C', 'D'}, {'E', 'F'}]

        expected_seqnum = ['1', '1', '2', '2']

        seq = ces.sequence(seq=_SEQUENCE_5, header=_HEADER_1)
        seq_2 = ces.sequence(seq=_SEQUENCE_6_2.replace("+", ""), header=_HEADER_2)

        oseq = ces.oligoseq(oligomer_id=seq.oligomer_id, imer={seq.name: seq})
        oseq.add_sequence(seq_2)

        obtained_seqnum = []
        for aset in chains:
            for ch in aset:
                obtained_seqnum.append(oseq.whatseq(ch))

        self.assertListEqual(expected_seqnum, obtained_seqnum)
