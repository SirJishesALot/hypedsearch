from typing import Sequence
from src.objects import Alignments, Spectrum, Database
import unittest
import os

from src import gen_spectra
from src.alignment import alignment, alignment_utils
from src import database
from src.preprocessing import merge_search, preprocessing_utils
from src import utils

class test_utils(unittest.TestCase):
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

    def test_gen_extensions(self):
        # test b extensions
        seq = 'DLQTLAL'
        ion = 'b'
        # Generate spectra for some peptide. Note that we are only testing this function's ability to generate extensions
        # from the database so the input spectrum is not relevant. 
        spectrum = Spectrum(gen_spectra.gen_spectrum('DLQTLALWSRM'), [], 0, 0, -1, gen_spectra.get_precursor('DLQTLALWSRM'), 1, '', '', {})
        # Import database (mouse_filtered.fasta) in NOD2_E3 and do all usual preprocessing
        db = database_preprocessing("/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/Hybrid_inputs/hybrid_nod2e3.mzML", "/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/raw_inputs/NOD2_E3/filtered_mouse_database.fasta", True, 5, 20, 0, 0.0, 20, 10, '', 1, 5, False, '', '')
        # Generate kmer set. This is going to be every kmer of size len(seq)

        # Generate all extensions
        b_ext = []
        
        b_ext += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'b', db)]
        
        # Check there are no missing extensions
        # These extensions were found manually
        b_test_ext = ['DLQTLALEVAQQK', 'DLQTLALEVARQK']
        self.assertEqual(sorted(b_ext), sorted(b_test_ext))
        #calculate extension length
        extension_len = utils.predicted_len_precursor(spectrum, seq) - len(seq)

        # test y extensions
        ion = 'y'
        # Generate spectra for some peptide. Note that we are only testing this function's ability to generate extensions
        # from the database so the input spectrum is not relevant. 
        spectrum = Spectrum(gen_spectra.gen_spectrum('DLQTLALWSRM'), [], 0, 0, -1, gen_spectra.get_precursor('DLQTLALWSRM'), 1, '', '', {})
        # Import database (mouse_filtered.fasta) in NOD2_E3 and do all usual preprocessing
        db = database_preprocessing("/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/Hybrid_inputs/hybrid_nod2e3.mzML", "/mnt/c/Users/Maxim/Documents/Layer_Lab/Database/raw_inputs/NOD2_E3/filtered_mouse_database.fasta", True, 5, 20, 0, 0.0, 20, 10, '', 1, 5, False, '', '')
        # Generate kmer set. This is going to be every kmer of size len(seq)

        # Generate all extensions
        y_ext = []
        
        y_ext += [x for x in alignment_utils.extend_non_hybrid(seq, spectrum, 'y', db)]
        
        # Check there are no missing extensions
        # These extensions were found manually
        y_test_ext = ['GGPGAGDLQTLAL', 'LGGSPGDLQTLAL']
        self.assertEqual(sorted(y_ext), sorted(y_test_ext))
        #calculate extension length
        extension_len = utils.predicted_len_precursor(spectrum, seq) - len(seq)


def database_preprocessing(
        spectra_files,
        database_file, 
        verbose: bool = True, 
        min_peptide_len: int = 5, 
        max_peptide_len: int = 20, 
        peak_filter: int = 0, 
        relative_abundance_filter: float = 0.0,
        ppm_tolerance: int = 20, 
        precursor_tolerance: int = 10, 
        digest: str = '',
        cores: int = 1,
        n: int = 5,
        DEBUG: bool = False, 
        truth_set: str = '', 
        output_dir: str = ''
        ):

        DEV = False
        truth = None

        fall_off = None

        spectra_files = [spectra_files]
        # build/load the database
        verbose and print('Loading database...')
        db = database.build(database_file)
        verbose and print('Done')

        # load all of the spectra
        verbose and print('Loading spectra...')
        spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(
            spectra_files, 
            ppm_tolerance,
            peak_filter=peak_filter, 
            relative_abundance_filter=relative_abundance_filter
        )
        verbose and print('Done')

        # get the boundary -> kmer mappings for b and y ions
        matched_masses_b, matched_masses_y, db = merge_search.match_masses(boundaries, db, max_peptide_len)

        # keep track of the alingment made for every spectrum
        results = {}
        return db