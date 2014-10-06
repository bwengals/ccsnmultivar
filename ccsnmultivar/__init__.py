__version__ = '0.0.2'

__all__ = ['basis','designmatrix','multivar']

# want to have users just do import ccsnmultivar as cm and have access to the following:

# imports of multivar classes
from ccsnmultivar.multivar import Multivar
# imports of designmatrix classes
from ccsnmultivar.designmatrix import DesignMatrix
# imports of basis classes
from ccsnmultivar.basis import PCA, KPCA
from ccsnmultivar.basis import load_waveforms


