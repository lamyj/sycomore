import io
import os
import subprocess
import tempfile
import zipfile

subprocess.check_call([
    "apt-get", "-y", "--no-install-recommends", "install",
    "cmake", "g++", "libboost-dev", "make", "ninja-build",
    "pybind11-dev", "python3-pybind11", "python3-dev", 
    "python3-requests", "libboost-test-dev", "python3-numpy"])

import requests

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
