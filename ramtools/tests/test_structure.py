import pytest
import numpy as np

from ramtools.tests.base_test import BaseTest
from ramtools.structure import calc_number_density
from ramtools.utils.io import get_fn


class TestStructure(BaseTest):
    def test_number_density(self, water_trj):
        area = water_trj.unitcell_lengths[0][0] * water_trj.unitcell_lengths[0][1]
        dim = 2
        box_range = [0, 1]
        n_bins = 5

        calc_number_density.calc_number_density(
            water_trj, area=area, dim=dim, box_range=box_range, n_bins=n_bins
        )
