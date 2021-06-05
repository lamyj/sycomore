import os
import site
import subprocess
import sys

sys.path.append(site.USER_SITE)
os.environ["PATH"] = os.pathsep.join([
    os.environ["PATH"], 
    os.path.join(site.USER_BASE, "bin"),
    os.path.join(site.USER_SITE, "../Scripts")])

subprocess.check_call([sys.executable, "setup.py", "bdist_wheel"])
