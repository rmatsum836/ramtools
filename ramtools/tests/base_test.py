import mdtraj as md
import pytest

from ramtools.utils.io import get_fn

class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def water_trj(self):
        trj = md.load(get_fn('spce.xtc'), top=get_fn('spce.gro'))[:100]

        return trj

    @pytest.fixture
    def il_trj(self):
        trj = md.load(get_fn('il.trr'), top=get_fn('il.gro'))[:100]

        return trj
