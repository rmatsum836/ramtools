import pytest

from ramtools.utils.io import get_fn
from ramtools.utils import gromacs_tools
from ramtools.utils.read_files import read_xvg
from ramtools.tests.base_test import BaseTest

class TestUtils(BaseTest):

    def test_read_xvg(self):
        xvg = read_xvg(get_fn('energy.xvg'))

    def test_create_itp(self):
        gromacs_tools.create_itp('init.itp', n_atoms=100, fx_const=100)

    def test_make_ndx(self):
        gromacs_tools.make_ndx(get_fn('water_spce.gro'), 'resname SOL', 'resname SOL', 'index.ndx')
