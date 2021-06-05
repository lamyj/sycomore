import subprocess
import sys

os.environ["PATH"] = os.pathsep.join([
    os.environ["PATH"], 
    os.path.join(site.USER_BASE, "bin"),
    os.path.join(site.USER_SITE, "../Scripts")])

subprocess.check_call([sys.executable, "setup.py", "bdist_wheel"])
