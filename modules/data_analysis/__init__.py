#import os
#import inspect
#import sys

#defines script path as the parent folder of the module (where the main script is located)
#filename = inspect.getframeinfo(inspect.currentframe()).filename
#__SCRIPT_PATH__ = os.path.dirname(os.path.abspath(os.path.abspath(os.path.join(filename, os.pardir))))

from .data_analysis import *