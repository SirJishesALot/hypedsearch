# This is an exact clone of identification.py with functions renamed for clarity and all code relating to creating an 
# alignment removed

from typing import Tuple
from src.objects import Database, Spectrum, MPSpectrumID, DEVFallOffEntry
from src.cppModules import gen_spectra
from src.utils import ppm_to_da, to_percent, overlap_intervals, hashable_boundaries, is_json, is_file
from src import utils
from src.scoring import scoring, mass_comparisons
from src.preprocessing import merge_search, preprocessing_utils
from src import database
from src.file_io import JSON

import time
import os
import copy
import json

# top results to keep for creating an alignment
TOP_X = 50

def database_and_spectra_preprocessing(
    spectra_folder: str, 
    database_file: str, 
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
) -> Tuple[dict, Spectrum]:
    '''Load in all the spectra and try to create an alignment for every spectrum

    :param spectra_files: file names of input spectra
    :type spectra_files: list
    :param database_file: file name of the fasta database
    :type database_file: str
    :param verbose: print progress to the console. 
        (default is True)
    :type verbose: bool
    :param min_peptide_len: the minimum length alignment to create
        (default is 5)
    :type min_peptide_len: int
    :param max_peptide_len: the maximum length alignment to create
        (default is 20)
    :type max_peptide_len: int
    :param peak_filter: If set to a number, this metric is used over the relative abundance filter. 
        The most abundanct X peaks to use in the alignment. 
        (default is 0)
    :type peak_filter: int
    :param relative_abundance_filter: If peak_filter is set, this parameter is ignored. The 
        relative abundance threshold (in percent as a decimal) a peak must be of the total 
        intensity to be used in the alignment. 
        (default is 0.0)
    :type relative_abundance_filter: float
    :param ppm_tolerance: the parts per million error allowed when trying to match masses
        (default is 20)
    :type ppm_tolerance: int
    :param precursor_tolerance: the parts per million error allowed when trying to match
        a calculated precursor mass to the observed precursor mass
        (default is 10)
    :type precurosor_tolerance: int
    :param digest: the type of digest used in the sample preparation. If left blank, 
        a digest-free search is performed. 
        (default is '')
    :type digest: str
    :param cores: the number of cores allowed to use in running the program. If a number 
        provided is greater than the number of cores available, the maximum number of 
        cores is used. 
        (default is 1)
    :type cores: int
    :param n: the number of aligments to keep per spectrum. 
        (default is 5)
    :type n: int
    :param DEBUG: DEVELOPMENT USE ONLY. Used only for timing of modules. 
        (default is False)
    :type DEBUG: bool
    :param truth_set: the path to a json file of the desired alignments to make for each spectrum. 
        The format of the file is {spectrum_id: {'sequence': str, 'hybrid': bool, 'parent': str}}. 
        If left an empty string, the program proceeds as normal. Otherwise results of the analysis
        will be saved in the file 'fall_off.json' saved in the output directory specified.
        (default is '')
    :type truth_set: str
    :param output_dir: the full path to the output directory to save all output files.
        (default is '')
    :type output_dir: str

    :returns: alignments for all spectra save in the form {spectrum.id: Alignments}
    :rtype: dict
    '''
    # get all the spectra file names
    spectra_files = []
    for (root, _, filenames) in os.walk(spectra_folder):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))

    # build/load the database
    verbose and print('Loading database...')
    db = database.build(database_file)
    verbose and print('Done')
    
    # # load all of the spectra
    # verbose and print('Loading spectra...')
    # spectra, boundaries, mz_mapping = preprocessing_utils.load_spectra(
    #     spectra_files, 
    #     ppm_tolerance,
    #     peak_filter=peak_filter, 
    #     relative_abundance_filter=relative_abundance_filter
    # )
    # verbose and print('Done')

    # # get the boundary -> kmer mappings for b and y ions
    # matched_masses_b, matched_masses_y, db = merge_search.match_masses(boundaries, db, max_peptide_len)

    # # if we only get 1 core, don't do the multiprocessing bit
    # if cores == 1:
    #     # go through and id all spectra
    #     for i, spectrum in enumerate(spectra):

    #         print(f'Creating alignment for spectrum {i+1}/{len(spectra)} [{to_percent(i+1, len(spectra))}%]', end='\r')

    #         # get b and y hits
    #         b_hits, y_hits = [], []
    #         for mz in spectrum.spectrum:

    #             # get the correct boundary
    #             mapped = mz_mapping[mz]
    #             b = boundaries[mapped]
    #             b = hashable_boundaries(b)

    #             if b in matched_masses_b:
    #                 b_hits += matched_masses_b[b]

    #             if b in matched_masses_y:
    #                 y_hits += matched_masses_y[b]

    return db, spectrum