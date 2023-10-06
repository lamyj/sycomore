import subprocess
import sys

subprocess.check_call([
    "cmake", "-S", ".", f"-DPython_EXECUTABLE={sys.executable}",
    "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_INSTALL_PREFIX={sys.exec_prefix}",
    "-DBUILD_TESTING=OFF", "-B", "build"])
subprocess.check_call([
    "cmake", "--build", "build", "--target", "install", "--config", "Release"])

