import os

#from amuse.lab import *

from modules import __SCRIPT_PATH__

#general output paths
__DATA_DIR__ = __SCRIPT_PATH__ + '/data/'

__ANIMATION_DIR__ = __DATA_DIR__ + '/animations/'
__FRAME_DIR__ = __ANIMATION_DIR__ + '/_frames/'

__MODEL_DIR__ = __DATA_DIR__ + '/models/'
__TEST_MODEL_DIR__ = __MODEL_DIR__ + '/test/'
__FULL_MODEL_DIR__ = __MODEL_DIR__ + '/full/'

__MERGER_DIR__ = __DATA_DIR__ + '/merger/'

__global_dirs__ = [__DATA_DIR__ ,
                   __ANIMATION_DIR__,
                   __FRAME_DIR__,
                   __MODEL_DIR__,
                   __TEST_MODEL_DIR__,
                   __FULL_MODEL_DIR__,
                   __MERGER_DIR__]
                   
for direct in __global_dirs__:
    if not os.path.exists(direct):
        os.makedirs(direct)
        
def galaxy_structures(galaxy, n_disk, n_bulge, n_halo=None):
    halo = galaxy[n_disk+n_bulge:]
    disk = galaxy[:n_disk]
    bulge = galaxy[n_disk:n_disk+n_bulge]
    return halo, disk, bulge