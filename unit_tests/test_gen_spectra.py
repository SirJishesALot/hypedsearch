from typing import Sequence
from src.constants import AMINO_ACIDS, DOUBLY_CHARGED_B_BASE, DOUBLY_CHARGED_Y_BASE, PROTON_MASS, SINGLY_CHARGED_B_BASE, SINGLY_CHARGED_Y_BASE, WATER_MASS
import unittest
import os
from src import gen_spectra

class test_gen_spectra(unittest.TestCase): 
    def setUp(self): 
        self.sequence = 'MALWAR'
    
    def test_b_ions(self):
        #test the b_ions function with the sequence, "MALWAR".
        target_sequence = [AMINO_ACIDS['M'], AMINO_ACIDS['M'] + AMINO_ACIDS['A'], 
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'],
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'],
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'] + AMINO_ACIDS['A'],
        AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'] + AMINO_ACIDS['A'] + AMINO_ACIDS['R']]
        #adding singly charged base
        singly_b_seq = [n + SINGLY_CHARGED_B_BASE for n in target_sequence]
        self.assertEqual(gen_spectra.b_ions(self.sequence, 1), singly_b_seq)
        #adding doubly charged base
        doubly_b_seq = [(n + DOUBLY_CHARGED_B_BASE) / 2 for n in target_sequence]
        #testing doubly charged
        self.assertEqual(gen_spectra.b_ions(self.sequence, 2), doubly_b_seq)
    
    def test_y_ions(self):
        #Test the y_ions function with the sequence, 'MALWAR".
        target_sequence = [AMINO_ACIDS['R'], AMINO_ACIDS['A'] + AMINO_ACIDS['R'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'] + AMINO_ACIDS['L'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'] + AMINO_ACIDS['L'] + AMINO_ACIDS['A'],
        AMINO_ACIDS['R'] + AMINO_ACIDS['A'] + AMINO_ACIDS['W'] + AMINO_ACIDS['L'] + AMINO_ACIDS['A'] + AMINO_ACIDS['M']]
        #add in single/doubly charged
        singly_charged_seq = [n + SINGLY_CHARGED_Y_BASE for n in target_sequence]
        self.assertEqual(gen_spectra.y_ions(self.sequence, 1), singly_charged_seq)
        #testing doubly charged
        doubly_charged_seq = [(n + DOUBLY_CHARGED_Y_BASE) / 2 for n in target_sequence]
        self.assertEqual(gen_spectra.y_ions(self.sequence, 2), doubly_charged_seq)

    def test_calc_masses(self):
        #Test the calc masses function with the sequence, "MALWAR".
        #The masses should be all of the b and y ion masses in one list.
        self.assertEqual(sorted(gen_spectra.calc_masses(self.sequence)), (sorted((gen_spectra.b_ions(self.sequence)) 
            + (gen_spectra.y_ions(self.sequence)))))
    
    def test_max_mass(self):
        #Test the max mass function with the sequence, "MALWAR"
        Expected_max = 728.379201
        charge = 1
        charge2 = 2
        self.assertEqual(gen_spectra.max_mass(self.sequence, 'b', charge), Expected_max + SINGLY_CHARGED_B_BASE)
        self.assertEqual(gen_spectra.max_mass(self.sequence, 'y', charge2), (Expected_max + DOUBLY_CHARGED_Y_BASE) / 2)
    
    def test_get_precursor(self):
        #Test the get_precursory function with the sequence, "MALWAR"
        Expected_prec = AMINO_ACIDS['M'] + AMINO_ACIDS['A'] + AMINO_ACIDS['L'] + AMINO_ACIDS['W'] + AMINO_ACIDS['A'] + AMINO_ACIDS['R'] + WATER_MASS
        self.assertEqual(gen_spectra.get_precursor(self.sequence), Expected_prec + PROTON_MASS)
    
    def test_gen_spectrum(self):
        #Test the gen_spectrum function with the sequence, "MALWAR"
        Expected_spectrum = gen_spectra.b_ions(self.sequence) + gen_spectra.y_ions(self.sequence)
        self.assertEqual(gen_spectra.gen_spectrum(self.sequence), Expected_spectrum)