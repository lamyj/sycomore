import io
import os
import re
import subprocess
import sys
import tempfile
import zipfile

data = subprocess.check_output(["apt-cache", "show", "pybind11-dev"])
match = re.search(r"^Version: ([\d.]+)", data.decode(), flags=re.M)
pybind11_version = tuple(int(x) for x in match.group(1).split("."))
if pybind11_version >= (2,2):
    pybind11 = ["pybind11-dev"]
else:
    pybind11 = []

os.environ["DEBIAN_FRONTEND"] = "noninteractive"
subprocess.check_call([
    "apt-get", "-y", "--no-install-recommends", "install",
    "cmake", "g++", "libboost-dev", "make", "ninja-build",
    *pybind11, "python3-dev", 
    "python3-requests", "libboost-test-dev", "python3-numpy"])

import requests

if not pybind11:
    with tempfile.TemporaryDirectory() as directory:
        response = requests.get("https://github.com/pybind/pybind11/archive/refs/tags/v2.6.2.zip")
        response.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            archive.extractall(directory)
        subprocess.check_call(
            [
                "cmake", 
                "-DPYBIND11_TEST=OFF", "-DPYTHON_EXECUTABLE={}".format(sys.executable),
                "-DCMAKE_INSTALL_PREFIX={}".format(os.path.expanduser("~/local")),
                "."], 
            cwd=os.path.join(directory, "pybind11-2.6.2"))
        subprocess.check_call(
            ["cmake", "--build", ".", "--config", "Release", "--target", "install"],
            cwd=os.path.join(directory, "pybind11-2.6.2"))

with tempfile.TemporaryDirectory() as directory:
    response = requests.get("https://github.com/xtensor-stack/xsimd/archive/refs/tags/7.5.0.zip")
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        archive.extractall(directory)
    subprocess.check_call(
        [
            "cmake", 
            "-DCMAKE_INSTALL_PREFIX={}".format(os.path.expanduser("~/local")),
            "."], 
        cwd=os.path.join(directory, "xsimd-7.5.0"))
    subprocess.check_call(
        ["cmake", "--build", ".", "--config", "Release", "--target", "install"],
        cwd=os.path.join(directory, "xsimd-7.5.0"))
