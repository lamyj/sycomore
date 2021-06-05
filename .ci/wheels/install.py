import io
import os
import site
import subprocess
import sys
import tempfile
import urllib.request
import zipfile

subprocess.check_call([
    sys.executable, "-m", "pip", "install", "--force-reinstall", "pip"])

subprocess.check_call([
    sys.executable, "-m", "pip", "install", "--user", 
    "setuptools", "wheel", "setuptools_scm",
    "cmake", "pybind11"])

os.environ["PATH"] = os.pathsep.join([
    os.environ["PATH"], 
    os.path.join(site.USER_BASE, "bin"),
    os.path.join(site.USER_SITE, "../Scripts")])

with tempfile.TemporaryDirectory() as directory:
    with urllib.request.urlopen("https://github.com/xtensor-stack/xsimd/archive/refs/tags/7.5.0.zip") as fd:
        archive = zipfile.ZipFile(io.BytesIO(fd.read()))
        archive.extractall(directory)
    subprocess.check_call(
        [
            "cmake", 
            "-DCMAKE_INSTALL_PREFIX={}".format(
                os.path.join(site.USER_BASE, "include")), 
            "."], 
        cwd=os.path.join(directory, "xsimd-7.5.0"))
    subprocess.check_call(
        ["cmake", "--build", ".", "--config", "Release", "--target", "install"],
        cwd=os.path.join(directory, "xsimd-7.5.0"))
