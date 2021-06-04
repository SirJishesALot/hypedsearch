from typing import Tuple
import unittest
import os
from unit_tests import database_preprocessing
from src import database
from src.preprocessing import merge_search, preprocessing_utils


class test_preprocessing(unittest.TestCase):
    def setUp(self) -> None:
        spectra_folder = "/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/Hybrid_inputs"
        database_file = '/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/raw_inputs/NOD2_E3/filtered_mouse_database.fasta'
        global spectra_files
        # get all the spectra file names
        spectra_files = []
        for (root, _, filenames) in os.walk(spectra_folder):
            for fname in filenames:
                spectra_files.append(os.path.join(root, fname))

        # build/load the database
        db = database.build(database_file)

    def test_load_spectra(self):
        #Test the preprocessing_utils.load_spectra function
        spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(spectra_files, 20, 0, 0.0)
        self.assertEqual((spectra, boundaries, mz_mapping) , (target_spectra, target_boundaries, target_mz_mapping))

    # def test_merge(self):
    #     #Testing the merge_search.merge function
    #     self.assertEqual(merge_search.merge(spectrum[0], indices, kmers, boundaries))

    # def test_make_database_set(self):
    #     # Testing the merge_search.make_database_set function

    # def test_match_masses(self):
    #     # Testing the merge_search.match_masses function
    #     self.assertEqual(merge_search.match_masses(boundaries, db), Tuple(dict, dict, db))