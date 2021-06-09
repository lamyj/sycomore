import io
import os
import site
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

subprocess.check_call([
    sys.executable, "-m", "pip", "install", "--upgrade", "pip"])

subprocess.check_call([
    sys.executable, "-m", "pip", "install", "--user", "--upgrade",
    "requests", "setuptools", "setuptools_scm", "wheel", 
    "cmake", "pybind11"])

if sys.platform.startswith("linux"):
    subprocess.check_call([
        sys.executable, "-m", "pip", "install", "--user", "--upgrade",
        "auditwheel"])
elif sys.platform == "darwin":
    subprocess.check_call([
        sys.executable, "-m", "pip", "install", "--user", "--upgrade",
        "delocate"])

sys.path.append(site.USER_SITE)
import requests

os.environ["PATH"] = os.pathsep.join([
    os.environ["PATH"], 
    os.path.join(site.USER_BASE, "bin"),
    os.path.join(site.USER_SITE, "../Scripts")])

with tempfile.TemporaryDirectory() as directory:
    response = requests.get("https://github.com/xtensor-stack/xsimd/archive/refs/tags/7.5.0.zip")
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        archive.extractall(directory)
    subprocess.check_call(
        ["cmake", "-DCMAKE_INSTALL_PREFIX={}".format(site.USER_BASE), "."], 
        cwd=os.path.join(directory, "xsimd-7.5.0"))
    subprocess.check_call(
        ["cmake", "--build", ".", "--config", "Release", "--target", "install"],
        cwd=os.path.join(directory, "xsimd-7.5.0"))
