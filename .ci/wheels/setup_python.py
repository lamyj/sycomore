import glob
import os
import subprocess
import sys
import tempfile
import urllib.request

version = sys.argv[1]

if sys.platform.startswith("linux"):
    # Assume we are running in manylinux
    roots = {
        "3.6": "/opt/python/cp36-cp36m/bin",
        "3.7": "/opt/python/cp37-cp37m/bin",
        "3.8": "/opt/python/cp38-cp38/bin",
        "3.9": "/opt/python/cp39-cp39/bin",
    }
    if version not in roots:
        raise NotImplementedError(
            "Unknown interpreter version: {0}".format(version))
    
    interpreter = glob.glob(os.path.join(roots[version], "python"))[0]
elif sys.platform == "darwin":
    # Assume x86_64
    urls = {
        "3.6": "https://www.python.org/ftp/python/3.6.8/python-3.6.8-macosx10.9.pkg",
        "3.7": "https://www.python.org/ftp/python/3.7.9/python-3.7.9-macosx10.9.pkg",
        "3.8": "https://www.python.org/ftp/python/3.8.10/python-3.8.10-macosx10.9.pkg",
        "3.9": "https://www.python.org/ftp/python/3.9.5/python-3.9.5-macos11.pkg",
    }
    if version not in urls:
        raise NotImplementedError("Unknown interpreter version: {0}".format(version))
    
    with tempfile.TemporaryDirectory() as directory:
        data = urllib.request.urlopen(urls[version]).read()
        path = os.path.join(directory, urls[version].split("/")[-1])
        with open(path, "wb") as fd:
            fd.write(data)
        subprocess.check_call([
            "sudo", "installer", "-pkg", path, "-target", "/"])
    
    frameworks = "/Library/Frameworks/Python.framework/Versions"
    interpreter = glob.glob(os.path.join(frameworks, version, "bin/python3"))[0]
elif sys.platform == "win32":
    versions = {
        "3.6": "3.6.8",
        "3.7": "3.7.9",
        "3.8": "3.8.10",
        "3.9": "3.9.5",
    }
    if version not in versions:
        raise NotImplementedError(
            "Unknown interpreter version: {0}".format(version))
    
    subprocess.check_call([
        "nuget", "install", "python", 
        "-Version", version, "-OutputDirectory", "C:\\wheel_python", 
        "-Verbosity", "quiet"])
    interpreter = glob.glob("C:\\wheel_python\\*\\tools\\python.exe")[0]
else:
    raise NotImplementedError(sys.platform)

if "GITHUB_ENV" in os.environ:
    with open(os.environ["GITHUB_ENV"], "a") as fd:
        fd.write("PYTHON={}\n".format(interpreter))
else:
    print(interpreter)
