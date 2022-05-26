"""Testing facilities for crops.io.taggers"""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

from crops.io import taggers as cit

import unittest

_HEADERS_LIST=""">sp|Q6GZX4|001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1
>tr|Q3SA23|Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1
>UniRef50_Q9K794 Putative AgrB-like protein n=2 Tax=Bacillus TaxID=1386 RepID=AGRB_BACHD
>UPI0000000005 status=active
>sp|P05067 archived from Release 18.0 01-MAY-1991 SV=3
>tr|Q55167 archived from Release 17.0 01-JUN-2001 SV=1
>3J4F:A|PDBID|CHAIN|SEQUENCE
>3j4f_D resolution: 8.60 experiment: EMIC release_date: 2013-07-24 [ 3 : ALL ]
>6ln2_A resolution: 3.20 experiment: XRAY release_date: 2020-03-18 [ 476678 : ALL ] ['29-1054'] <SEQSE>29,1054<SEQSE> <100>1<100> <95>1<95> <90>1<90> <70>0<70> <50>0<50>
>6AVG_2|Chains B,C|T-cell receptor alpha variable 4,TCR alpha chain|Homo sapiens
>abcde|Chains j,k,l
>pdb|6AVG_2|B C
>6DWU_1|Chains AA[auth AC], AB[auth AE], AC[auth AG], A[auth AA], CA[auth AK], CB[auth AM], CC[auth AO], C[auth AI], EA[auth AS], EB[auth AU], EC[auth BA], E[auth AQ], GA[auth BE], GB[auth BG], GC[auth BI], G[auth BC], IA[auth BM], IB[auth BO], IC[auth BQ], I[auth BK], KA[auth BU], KB[auth CA], K[auth BS], MA[auth CE], MB[auth CG], M[auth CC], OA[auth CK], OB[auth CM], O[auth CI], QA[auth CQ], QB[auth CS], Q[auth CO], SA[auth DA], SB[auth DC], S[auth CU], UA[auth DG], UB[auth DI], U[auth DE], WA[auth DM], WB[auth DO], W[auth DK], YA[auth DS], YB[auth DU], Y[auth DQ]|Cationic trypsin|Bos taurus (9913)
>crops|HJ41_1|Chains A,B,C|these are some comments
>uc30-1808-30|Representative=A0A094PJX3 n=1 Descriptions=[Uncharacterized protein] Members=A0A094PJX3
>uc90-1808-39|Representative=U1GGI9 n=4 Descriptions=[ATP-dependent Clp protease adapter protein ClpS|ATP-dependent Clp protease adapter ClpS] Members=U7JT03,A0A239WXC9,U1GGI9,A0A2W5CTJ1
>uc90-1808-14|Representative=N1NU19 n=76 Descriptions=[Conjugal transfer protein TraV|Type IV conjugative transfer system protein TraV|Sex pilus assembly protein|Conjugative transfer protein TraV|Sex pilus assembly|Type IV conjugative transfer system lipoprotein (TraV)] Members=A0A126JTC8,A0A0E0SXV4,A0A0F1AU30,A0A167SUV6,A0PB04,U9YS98,A0A144XRS0,D2Y9Z2,A0A142BNX4,C4NVF2,A0A0J0IP13,A0A0K3RSM2,A0A097J0G2,A4IUB3,W1JB02,A0A0H3YD74,A0A1W6ARV1,T0QDZ9,A0A0V9GPN4,A0A0V9J4A8,E7DBH5,A0A2S4RFM4,C4NUX0,A0A2R4DIJ7,A0A023VGG0,A0A0A0CK90,U3PHU8,A0A096Y6F3,A0A0K2S451,B6VP02,A0A2I8R9F4,A0A2R4KLJ4,K0H5V6,K4Y465,D3VLT8,F5BPS4,A0A0B6VLZ4,T0P0P2,A0A0D7LBA5,E1ITU2,A0A0N9NKB1,A0A1C0D4B4,H9TK43,A0A0F0SXM7,A0A0A8NVJ2,A0A220UU09,K4W3W9,F2Z7X9,A0A0M1V7U1,A0A0W1W990,C6GA28,A0A2I8XNZ0,A0A2K9V063,A0A0R4CY03,A0A2S1XWH3,A0A1Q8YU31,V9SJY4,A0A0C4XWE1,A0A2L0TKR8,A0A2N1EL23,A0A0R9Q1C9,K4N9S7,N1NU19,A0A1B3B6M7,A0A0D5XR52,A0A0F0T558,A0A2I7QEW9,A0A1Y0F4W8,A0A221KL77,E7DBG1,A0A064D5Y7,A0A1U7FUF6,A0A0H3AEJ2,D3VGJ5,A0A1W6AR35,A0A1J0E5N4
"""

class TestCropsTaggers(unittest.TestCase):
    def test_retrieve_id_1(self):
        expected_tags = {'mainid': 'Q6GZX4', 'chains': {'Q6GZX4'},
                         'seqid': '1', 'source': 'UniProtKB/SwissProt',
                         'comments': '001R_FRG3G Putative transcription factor 001R OS=Frog virus 3 (isolate Goorha) OX=654924 GN=FV3-001R PE=4 SV=1'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[0])

        self.assertDictEqual(retrieved_tags, expected_tags)


    def test_retrieve_id_2(self):
        expected_tags = {'mainid': 'Q3SA23', 'chains': {'Q3SA23'},
                         'seqid': '1', 'source': 'UniProtKB/TrEMBL',
                         'comments': 'Q3SA23_9HIV1 Protein Nef (Fragment) OS=Human immunodeficiency virus 1  OX=11676 GN=nef PE=3 SV=1'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[1])

        self.assertDictEqual(retrieved_tags, expected_tags)


    def test_retrieve_id_3(self):
        expected_tags = {'mainid': 'Q9K794', 'chains': {'Q9K794'},
                         'seqid': '1', 'source': 'UniRef50',
                         'comments': 'Putative AgrB-like protein n=2 Tax=Bacillus TaxID=1386 RepID=AGRB_BACHD'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[2])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_4(self):
        expected_tags = {'mainid': 'UPI0000000005', 'chains': {'UPI0000000005'},
                         'seqid': '1', 'source': 'UniParc',
                         'comments': 'status=active'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[3])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_5(self):
        expected_tags = {'mainid': 'P05067', 'chains': {'P05067'},
                         'seqid': '1', 'source': 'UniProtKB/SwissProt (archived)',
                         'comments': 'archived from Release 18.0 01-MAY-1991 SV=3'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[4])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_6(self):
        expected_tags = {'mainid': 'Q55167', 'chains': {'Q55167'},
                         'seqid': '1', 'source': 'UniProtKB/TrEMBL (archived)',
                         'comments': 'archived from Release 17.0 01-JUN-2001 SV=1'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[5])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_7(self):
        expected_tags = {'mainid': '3j4f', 'chains': {'A'},
                         'seqid': None, 'source': 'RCSB PDB',
                         'comments': 'PDBID|CHAIN|SEQUENCE'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[6])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_8(self):
        expected_tags = {'mainid': '3j4f', 'chains': {'D'},
                         'seqid': None, 'source': 'MrBUMP',
                         'comments': 'resolution: 8.60 experiment: EMIC release_date: 2013-07-24 [ 3 : ALL ]'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[7])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_9(self):
        expected_tags = {'mainid': '6ln2', 'chains': {'A'},
                         'seqid': None, 'source': 'MrBUMP',
                         'comments': "resolution: 3.20 experiment: XRAY release_date: 2020-03-18 [ 476678 : ALL ] ['29-1054'] <SEQSE>29,1054<SEQSE> <100>1<100> <95>1<95> <90>1<90> <70>0<70> <50>0<50>"}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[8])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_10(self):
        expected_tags = {'mainid': '6avg', 'chains': {'B', 'C'},
                         'seqid': '2', 'source': 'RCSB PDB',
                         'comments': 'T-cell receptor alpha variable 4,TCR alpha chain|Homo sapiens'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[9])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_11(self):
        expected_tags = {'mainid': 'abcde', 'chains': {'j', 'k', 'l'},
                         'seqid': None, 'source': 'RCSB PDB', 'comments': None}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[10])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_12(self):
        expected_tags = {'mainid': '6avg', 'chains': {'B', 'C'},
                         'seqid': '2', 'source': 'PDBe', 'comments': None}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[11])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_13(self):
        expected_tags = {'mainid': '6dwu',
                         'chains': {'AA', 'AC', 'AE', 'AG', 'AI', 'AK', 'AM',
                                    'AO', 'AQ', 'AS', 'AU', 'BA', 'BC', 'BE',
                                    'BG', 'BI', 'BK', 'BM', 'BO', 'BQ', 'BS',
                                    'BU', 'CA', 'CC', 'CE', 'CG', 'CI', 'CK',
                                    'CM', 'CO', 'CQ', 'CS', 'CU', 'DA', 'DC',
                                    'DE', 'DG', 'DI', 'DK', 'DM', 'DO', 'DQ',
                                    'DS', 'DU'},
                         'seqid': '1',
                         'source': 'RCSB PDB',
                         'comments': 'Cationic trypsin|Bos taurus (9913)'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[12])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_14(self):
        expected_tags = {'mainid': 'hj41', 'chains': {'A', 'B', 'C'},
                         'seqid': '1', 'source': 'CROPS',
                         'comments': 'these are some comments'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[13])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_15(self):
        expected_tags = {'mainid': 'A0A094PJX3', 'chains': {'A0A094PJX'},
                         'seqid': '30', 'source': 'UniClust30_2018_08',
                         'comments': 'Representative=A0A094PJX3 n=1 Descriptions=[Uncharacterized protein] Members=A0A094PJX3'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[14])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_16(self):
        expected_tags = {'mainid': 'U1GGI9',
                         'chains': {'A0A239WXC9', 'A0A2W5CTJ', 'U1GGI9', 'U7JT03'},
                         'seqid': '39',
                         'source': 'UniClust90_2018_08',
                         'comments': 'Representative=U1GGI9 n=4 Descriptions=[ATP-dependent Clp protease adapter protein ClpS|ATP-dependent Clp protease adapter ClpS] Members=U7JT03,A0A239WXC9,U1GGI9,A0A2W5CTJ1'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[15])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_retrieve_id_17(self):
        expected_tags = {'mainid': 'N1NU19',
                         'chains': {'A0A023VGG0', 'A0A064D5Y7', 'A0A096Y6F3', 'A0A097J0G2',
                                    'A0A0A0CK90', 'A0A0A8NVJ2', 'A0A0B6VLZ4', 'A0A0C4XWE1',
                                    'A0A0D5XR52', 'A0A0D7LBA5', 'A0A0E0SXV4', 'A0A0F0SXM7',
                                    'A0A0F0T558', 'A0A0F1AU30', 'A0A0H3AEJ2', 'A0A0H3YD74',
                                    'A0A0J0IP13', 'A0A0K2S451', 'A0A0K3RSM2', 'A0A0M1V7U1',
                                    'A0A0N9NKB1', 'A0A0R4CY03', 'A0A0R9Q1C9', 'A0A0V9GPN4',
                                    'A0A0V9J4A8', 'A0A0W1W990', 'A0A126JTC8', 'A0A142BNX4',
                                    'A0A144XRS0', 'A0A167SUV6', 'A0A1B3B6M7', 'A0A1C0D4B4',
                                    'A0A1J0E5N', 'A0A1Q8YU31', 'A0A1U7FUF6', 'A0A1W6AR35',
                                    'A0A1W6ARV1', 'A0A1Y0F4W8', 'A0A220UU09', 'A0A221KL77',
                                    'A0A2I7QEW9', 'A0A2I8R9F4', 'A0A2I8XNZ0', 'A0A2K9V063',
                                    'A0A2L0TKR8', 'A0A2N1EL23', 'A0A2R4DIJ7', 'A0A2R4KLJ4',
                                    'A0A2S1XWH3', 'A0A2S4RFM4', 'A0PB04', 'A4IUB3',
                                    'B6VP02', 'C4NUX0', 'C4NVF2', 'C6GA28', 'D2Y9Z2',
                                    'D3VGJ5', 'D3VLT8', 'E1ITU2', 'E7DBG1', 'E7DBH5',
                                    'F2Z7X9', 'F5BPS4', 'H9TK43', 'K0H5V6', 'K4N9S7',
                                    'K4W3W9', 'K4Y465', 'N1NU19', 'T0P0P2', 'T0QDZ9',
                                    'U3PHU8', 'U9YS98', 'V9SJY4', 'W1JB02'},
                         'seqid': '14',
                         'source': 'UniClust90_2018_08',
                         'comments': 'Representative=N1NU19 n=76 Descriptions=[Conjugal transfer protein TraV|Type IV conjugative transfer system protein TraV|Sex pilus assembly protein|Conjugative transfer protein TraV|Sex pilus assembly|Type IV conjugative transfer system lipoprotein (TraV)] Members=A0A126JTC8,A0A0E0SXV4,A0A0F1AU30,A0A167SUV6,A0PB04,U9YS98,A0A144XRS0,D2Y9Z2,A0A142BNX4,C4NVF2,A0A0J0IP13,A0A0K3RSM2,A0A097J0G2,A4IUB3,W1JB02,A0A0H3YD74,A0A1W6ARV1,T0QDZ9,A0A0V9GPN4,A0A0V9J4A8,E7DBH5,A0A2S4RFM4,C4NUX0,A0A2R4DIJ7,A0A023VGG0,A0A0A0CK90,U3PHU8,A0A096Y6F3,A0A0K2S451,B6VP02,A0A2I8R9F4,A0A2R4KLJ4,K0H5V6,K4Y465,D3VLT8,F5BPS4,A0A0B6VLZ4,T0P0P2,A0A0D7LBA5,E1ITU2,A0A0N9NKB1,A0A1C0D4B4,H9TK43,A0A0F0SXM7,A0A0A8NVJ2,A0A220UU09,K4W3W9,F2Z7X9,A0A0M1V7U1,A0A0W1W990,C6GA28,A0A2I8XNZ0,A0A2K9V063,A0A0R4CY03,A0A2S1XWH3,A0A1Q8YU31,V9SJY4,A0A0C4XWE1,A0A2L0TKR8,A0A2N1EL23,A0A0R9Q1C9,K4N9S7,N1NU19,A0A1B3B6M7,A0A0D5XR52,A0A0F0T558,A0A2I7QEW9,A0A1Y0F4W8,A0A221KL77,E7DBG1,A0A064D5Y7,A0A1U7FUF6,A0A0H3AEJ2,D3VGJ5,A0A1W6AR35,A0A1J0E5N4'}

        retrieved_tags = cit.retrieve_id(_HEADERS_LIST[16])

        self.assertDictEqual(retrieved_tags, expected_tags)

    def test_makeheader_1(self):
        tags = ['MAINID', '1', {'A', 'B', 'C'}, 'MADEUP', 'A made up protein information']

        expected_header = '>crops|MAINID_1|Chains A,B,C|Source: MADEUP|A made up protein information'

        new_header = cit.makeheader(mainid=tags[0], seqid=tags[1], chains=tags[2],
                                   source=tags[3], extrainfo=tags[4])

        self.assertEqual(new_header, expected_header)

    def test_makeheader_2(self):
        tags = ['MAINID', '1', {'A', 'B', 'C'}, 'MADEUP', 'A made up protein information']

        expected_header = '>MAINID_1|Chains A,B,C'

        new_header = cit.makeheader(mainid=tags[0], seqid=tags[1], chains=tags[2],
                                    source=tags[3], extrainfo=tags[4], short=True)

        self.assertEqual(new_header, expected_header)
