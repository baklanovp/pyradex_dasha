import sys
from os.path import dirname, abspath, join
import numpy as np
import astropy.units as u

ROOT_DIRECTORY = dirname(dirname(dirname(abspath(__file__))))
sys.path.append(join(ROOT_DIRECTORY))
print(ROOT_DIRECTORY)

from pyradex.core import Radex
from pyradex import synthspec

def test_synthspec():
    R = Radex(column=1e15, temperature=20, density=1e3)
    R.run_radex()
    wcs = np.linspace(103.0, 103.1, 1000)*u.GHz
    S = synthspec.SyntheticSpectrum.from_RADEX(wcs, R, linewidth=10*u.km/u.s)
