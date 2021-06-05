import os
import subprocess
import sys

import setup_python

version = sys.argv[1]
interpreter = setup_python.setup_python(version)
process = subprocess.check_call(
    [interpreter, "install.py"], cwd=os.path.dirname(__file__))
process = subprocess.check_call(
    [interpreter, "build.py"], cwd=os.path.dirname(__file__))
