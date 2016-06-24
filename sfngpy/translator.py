from spectral_cube import SpectralCube
import astropy.units as u
from radio_beam import Beam
from scipy.ndimage.interpolation import map_coordinates
import numpy as np
# Default behaviours
defaultBeam = Beam(major=15 * u.arcsec, minor=15 * u.arcsec)


def SampleWithConvolution(file, positions, beam=defaultBeam,
                          order=1, **kwargs):
    s = SpectralCube.read(file)
    spaxis = s.spectral_axis.value
    spaxis.shape += (1,)
    spaxis_ones = np.ones_like(spaxis)
    s2 = s.convolve_to(beam)
    ravals = spaxis_ones * positions.ra.value
    decvals = spaxis_ones * positions.dec.value
    vvals = spaxis * np.ones_like(positions.ra.value)
    x, y, v = s.wcs.all_world2pix(ravals, decvals, vvals, 0)
    output = map_coordinates(s2.filled_data[:], [v, y, x],
                             order=order, **kwargs)
    # import pdb; pdb.set_trace()
    return output
