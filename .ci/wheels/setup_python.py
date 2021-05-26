import sys

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
        raise NotImplementedError("Unknown interpreter version: {0}".format(version))
    
    print(roots[version])
elif sys.platform == "darwin":
    import pip
    pip.main(["install", "cibuildwheel"])
    import cibuildwheel.macos
    
    # Assume x86_64
    # From https://github.com/pypa/cibuildwheel/blob/v1.11.0/cibuildwheel/resources/build-platforms.toml#L37-L50
    urls = {
        "3.6": "https://www.python.org/ftp/python/3.6.8/python-3.6.8-macosx10.9.pkg",
        "3.7": "https://www.python.org/ftp/python/3.7.9/python-3.7.9-macosx10.9.pkg",
        "3.8": "https://www.python.org/ftp/python/3.8.10/python-3.8.10-macosx10.9.pkg",
        "3.9": "https://www.python.org/ftp/python/3.9.5/python-3.9.5-macos11.pkg",
    }
    if version not in urls:
        raise NotImplementedError("Unknown interpreter version: {0}".format(version))
    
    path = cibuildwheel.macos.install_cpython(version, urls[version])
    print(cibuildwheel.macos.SYMLINKS_DIR)
else:
    raise NotImplementedError(sys.platform)
