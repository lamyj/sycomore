import glob
import os
import site
import subprocess
import sys

sys.path.append(site.USER_SITE)
os.environ["PATH"] = os.pathsep.join([
    os.environ["PATH"], 
    os.path.join(site.USER_BASE, "bin"),
    os.path.join(site.USER_SITE, "../Scripts")])
os.environ["CMAKE_PREFIX_PATH"] = os.pathsep.join([
    site.USER_BASE, os.path.join(site.USER_SITE, "pybind11")])

subprocess.check_call([sys.executable, "setup.py", "bdist_wheel"])
for wheel in glob.glob("dist/*whl"):
    if sys.platform.startswith("linux"):
        subprocess.check_call(["auditwheel", "repair", "-w", "dist", wheel])
        os.unlink(wheel)
    elif sys.platform == "darwin":
        subprocess.check_call(["delocate-listdeps", wheel])
        subprocess.check_call(["delocate-wheel", wheel])
