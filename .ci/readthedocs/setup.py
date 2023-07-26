import io
import os
import subprocess
import tempfile
import zipfile

import requests

xt = "https://github.com/xtensor-stack"
cmake_projects = {
    "xsimd": ("10.0.0", f"{xt}/xsimd/archive/refs/tags/{{}}.zip"),
    "xtl": ("0.7.5", f"{xt}/xtl/archive/refs/tags/{{}}.zip"),
    "xtensor": ("0.24.4", f"{xt}/xtensor/archive/refs/tags/{{}}.zip"),
    "xtensor-python": ("0.26.1", f"{xt}/xtensor-python/archive/refs/tags/{{}}.zip"),
}

local = os.path.expanduser("~/local")
cmake = ["cmake", f"-DCMAKE_INSTALL_PREFIX={local}", "."]
install = ["cmake", "--build", ".", "--config", "Release", "--target", "install"]

for name, (version, url) in cmake_projects.items():
    with tempfile.TemporaryDirectory() as directory:
        wd = os.path.join(directory, f"{name}-{version}")
        
        response = requests.get(url.format(version))
        response.raise_for_status()
        
        with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
            archive.extractall(directory)
        
        subprocess.check_call(cmake, cwd=wd)
        subprocess.check_call(install, cwd=wd)
