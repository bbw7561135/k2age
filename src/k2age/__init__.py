## @package k2age
#
import logging
from os.path import dirname
from binary import *
from star import *
from tracks import *

log_file = dirname(binary.__file__) + '/../k2age.log'
logging.basicConfig(filename=log_file, level=logging.DEBUG)
logging.info('NEW RUN')
