import sys
import os

import tools
from tools.Config import *

tools.do_and_move_parFile()

os.chdir(LG_dir)
os.system('make clean')
os.system('make')

run = ('./L-Galaxies '+LG_inParFile)
os.system(run)
os.chdir(here)
