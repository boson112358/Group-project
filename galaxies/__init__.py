import os
import inspect
import sys

#finds script path
filename = inspect.getframeinfo(inspect.currentframe()).filename
__SCRIPT_PATH__ = os.path.dirname(os.path.abspath(os.path.abspath(os.path.join(filename, os.pardir))))

from .galaxies import *