from typing import Sequence
from src.objects import Spectrum
from unit_tests.test_alignment import database_preprocessing
from unit_tests import database_preprocessing
import unittest
from src.scoring import mass_comparisons, scoring
from src import database, gen_spectra, utils
from src.preprocessing import merge_search, preprocessing_utils

class test_scoring(unittest.TestCase):
    def setUp(self):
        self.tosortkey = [
            {'a': 'b', 'b': 2}, 
            {'a': 'e', 'b': 18},
            {'a': 'd', 'b': 11},
            {'a': 'a', 'b': -1}, 
            {'a': 'c', 'b': 3}    
        ]
        self.tosortindex = [
            ('b', 2), ('e', 18), ('d', 11), ('a', -1), ('c', 3)
        ]
    # To use once more rigorous framework is developed

    spectra_folder = "/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/Hybrid_inputs"
    db_file = "/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/raw_inputs/NOD2_E3/filtered_mouse_database.fasta"
    #import the database and observed spectrum
    global db, spectrum

    db = database_preprocessing.database_and_spectra_preprocessing(spectra_folder, db_file)
    spectrum = Spectrum(spectrum=[70.06481170654297, 72.07975006103516, 86.09618377685547, 120.08119201660156, 129.10186767578125, 132.10205078125, 183.14903259277344, 211.14425659179688, 231.06002807617188, 302.0948181152344, 316.1502380371094, 342.23724365234375, 343.2393493652344, 344.1456298828125, 443.20904541015625, 617.3176879882812, 635.3029174804688, 665.8477783203125, 732.345458984375, 732.85546875, 1134.5634765625, 1135.568603515625, 1233.637939453125, 1316.615234375, 1463.6756591796875], abundance=[5611.4921875, 759.39697265625, 13583.19140625, 1904.96875, 377.9483947753906, 4277.7255859375, 2619.271484375, 3084.39013671875, 1280.637939453125, 413.639892578125, 1418.86376953125, 6746.6708984375, 1324.3463134765625, 2561.01708984375, 949.918212890625, 388.83099365234375, 522.100830078125, 384.86090087890625, 1058.61181640625, 387.35174560546875, 804.1146240234375, 738.4801025390625, 975.6796875, 496.6000061035156, 672.422119140625], total_intensity=53342.53186035156, ms_level=2, scan_number=256, precursor_mass=640.007114, precursor_charge=3, file_name='/home/naco3124/raw_inputs/NOD2_E3/mzml/NOD2_E3.mzML', id='NOD2_E3.20155.20196.3.pkl', other_metadata={})

    def test_optimized_compare_masses(self):
        #Test the scoring.optimized_compare_masses function
        #easy test case:
        self.assertEqual(mass_comparisons.optimized_compare_masses([1, 2, 4], [1, 3, 4], 1, False), 2)
        #testing tolerance:
        self.assertEqual(mass_comparisons.optimized_compare_masses([1, 2, 4], [1, 3, 4], 1000000, False), 3)
        #test the optimized compare masses function with some input spectra and a database of hybrids

    def test_score_sequence(self):
        #Test the scoring.score_sequence function
        #easy test case:
        self.assertEqual(scoring.score_sequence([1, 2, 4], [1, 3, 4], 1, False), 2)
        #testing tolerance:
        self.assertEqual(scoring.score_sequence([1, 2, 4], [1, 3, 4], 1000000, False), 3)

    def test_hybrid_score(self):
        #Test the scoring.hybrid_score function
        hybrid_seq = 'ABC-DEF'
        lesser_point = .5
        greater_point = 1.0
        sample_spectrum = Spectrum(gen_spectra.gen_spectrum('ABCDEF'))
        # say our b ions found are A, C, E
        # and y ions found are D, A
        # our scoring then works like
        # .5(bA) + .5(bC) + 1(bE) + .5 (yD) + 1(yA) 
        self.assertEqual(scoring.hybrid_score(sample_spectrum, hybrid_seq, 20, lesser_point, greater_point), 3.5)

    def test_precursor_distance(self):
        #Test the scoring.precursor_distance function
        observed_precursor = 10
        reference_precursor = 15
        expected_precursor_dist = 5
        self.assertEqual(scoring.precursor_distance(observed_precursor, reference_precursor), expected_precursor_dist)
        #test the other way
        self.assertEqual(scoring.precursor_distance(reference_precursor, observed_precursor), expected_precursor_dist)

    def test_total_mass_error(self):
        #Test the scoring.total_mass_error function
        #Case 1, both spectra are the exact same
        same_spec = Spectrum(gen_spectra.gen_spectrum('DLTQLAL'))
        self.assertEqual(0, scoring.total_mass_error(same_spec, 'DLTQLAL', 0))
        #Case 2, spectra and ideal are different
        self.assertEqual(scoring.total_mass_error(spectrum, 'DDVALYNFSKYFIPLL', 20), 0.011992253671792241)

    def test_digest_score(self):
        #Test the scoring.digest_score function
        #Testing with a non-hybrid
        self.assertEqual(scoring.digest_score('DAP', db, 'trypsin'), 2)
        #Testing with a cis-spliced hybrid
        self.assertEqual(scoring.digest_score('DAP-GSTMYPGIADR', db, 'trypsin'), 2)
