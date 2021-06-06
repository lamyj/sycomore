import os
import subprocess
import sys

import setup_python

if sys.version_info < (2, 6):
    def check_call(*args, **kwargs):
        retcode = subprocess.call(*args, **kwargs)
        if retcode:
            raise Exception()
        return 0
    subprocess.check_call = check_call

version = sys.argv[1]
interpreter = setup_python.setup_python(version)
subprocess.check_call([interpreter, ".ci/wheels/install.py"])
subprocess.check_call([interpreter, ".ci/wheels/build.py"])
subprocess.check_call([interpreter, ".ci/wheels/post_build.py"])
