import mdtraj as md
import pytest

from ramtools.utils.io import get_fn

class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def water_trj(self):
        trj = md.load(get_fn('spce.trr'), top=get_fn('water_spce.gro'))

        return trj

    @pytest.fixture
    def gph_pore(self):
        trj = md.load(get_fn('nvt_small.trr'), top=get_fn('nvt.gro'))
  
        return trj
