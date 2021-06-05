import glob
import os
import shutil
import subprocess
import sys
import tempfile

if sys.version_info >= (3, 5):
    import urllib.request
    urlopen = urllib.request.urlopen
else:
    import urllib
    urlopen = urllib.urlopen

def setup_python(version):
    if sys.platform.startswith("linux"):
        return setup_linux(version)
    elif sys.platform == "darwin":
        return setup_macos(version)
    elif sys.platform == "win32":
        return setup_windows(version)
    else:
        raise NotImplementedError(sys.platform)

def setup_linux(version):
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
    
    return glob.glob(os.path.join(roots[version], "python"))[0]

def setup_macos(version):
    # Assume x86_64
    
    root = "https://www.python.org/ftp/python"
    urls = {
        "3.6": root+"/3.6.8/python-3.6.8-macosx10.9.pkg",
        "3.7": root+"/3.7.9/python-3.7.9-macosx10.9.pkg",
        "3.8": root+"/3.8.10/python-3.8.10-macosx10.9.pkg",
        "3.9": root+"/3.9.5/python-3.9.5-macos11.pkg",
    }
    if version not in urls:
        raise NotImplementedError(
            "Unknown interpreter version: {0}".format(version))
    
    directory = tempfile.mkdtemp()
    try:
        data = urlopen(urls[version]).read()
        path = os.path.join(directory, urls[version].split("/")[-1])
        # NOTE: manylinux1 has Python 2.4, causing a syntax error when using
        # "with open(...) as fd".
        fd = open(path, "wb")
        fd.write(data)
        fd.close()
        
        subprocess.check_call([
            "sudo", "installer", "-pkg", path, "-target", "/"])
    finally:
        shutil.rmtree(directory)
    
    frameworks = "/Library/Frameworks/Python.framework/Versions"
    return glob.glob(os.path.join(frameworks, version, "bin/python3"))[0]

def setup_windows(version):
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
    
    return glob.glob("C:\\wheel_python\\*\\tools\\python.exe")[0]
