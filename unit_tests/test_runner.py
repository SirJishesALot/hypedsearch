import unittest
import os
from src import utils
from src import gen_spectra
from src.objects import Spectrum

boundaries1 = [[2,10], [3,7]]
print(utils.hashable_boundaries(boundaries1))
