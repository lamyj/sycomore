import os
import subprocess
import sys

import cibuildwheel
import cibuildwheel.macos
configurations = cibuildwheel.macos.get_python_configurations(
    cibuildwheel.util.BuildSelector(
        build_config=sys.argv[1], skip_config=""), 
        [cibuildwheel.architecture.Architecture.x86_64])
for configuration in configurations:
    path = cibuildwheel.macos.install_cpython(
        configuration.version, configuration.url)
    python = os.path.join(path, "python3")
    pip = os.path.join(path, "pip3")
    subprocess.check_call([
        python, "-m", "pip", "install", "--force-reinstall", "pip"])
    subprocess.check_call([
        pip, "install", "setuptools", "wheel", "setuptools_scm"])
    subprocess.check_call([python, "setup.py", "bdist_wheel"])
