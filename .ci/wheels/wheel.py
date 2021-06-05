import os
import subprocess
import sys

import setup_python

version = sys.argv[1]
interpreter = setup_python.setup_python(version)
process = subprocess.check_call([interpreter, ".ci/wheels/install.py"])
process = subprocess.check_call([interpreter, ".ci/wheels/build.py"])
