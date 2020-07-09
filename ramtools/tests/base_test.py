import mdtraj as md
import pytest

from ramtools.utils.io import get_fn

class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def water_trj(self):
        trj = md.load(get_fn('tip3p.xtc'), top=get_fn('tip3p.gro'))

        return trj
