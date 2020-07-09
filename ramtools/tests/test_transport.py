import pytest
import numpy as np

from ramtools.tests.base_test import BaseTest
from ramtools.transport import calc_transport
from ramtools.utils.io import get_fn

class TestTransport(BaseTest):

    def test_conductivity(self):
        n_mol = 240
        volume = 141
        D_cat = 2e-10
        D_an = 2e-10

        calc_transport.calc_conductivity(n_mol, volume, D_cat, D_an)

    def test_hfshear(self, il_trj):
        edr = get_fn('il.edr')

        calc_transport.calc_hfshear(edr, il_trj, temperature=300)
