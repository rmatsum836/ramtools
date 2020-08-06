import pytest
import numpy as np
import tempfile

from ramtools.tests.base_test import BaseTest
from ramtools.structure import calc_number_density, calc_angle_distribution
from ramtools.utils.io import get_fn


class TestStructure(BaseTest):
    @pytest.mark.parametrize("n_bins, shift", [(100, True), (100, False), (101, True), (101, False)])
    def test_number_density(self, gph_pore, n_bins, shift):
        area = gph_pore.unitcell_lengths[0][0] * gph_pore.unitcell_lengths[0][1]
        dim = 1
        box_range = [0.837, 2.837]
        n_bins = n_bins

        calc_number_density.calc_number_density(
            gph_pore, 
            area=area,
            dim=dim,
            shift=shift,
            box_range=box_range,
            n_bins=n_bins
        )

    @pytest.mark.parametrize("shift", [True, False])
    def test_s_order_parameter(self, gph_pore_mda, shift):
        box_range = [8.37+3, 28.37+3]
        calc_angle_distribution.calc_water_order_parameter(gph_pore_mda[0],
            gph_pore_mda[1], bin_width=1, cutoffs=box_range, dim=1, filepath=tempfile.gettempdir(),
            shift=shift)
