import subprocess

# NOTE: restrict xsimd version on macOS due to
# https://github.com/xtensor-stack/xsimd/issues/810
xsimd_version = "<9.0.0,>9.0.1" if sys.platform == "darwin" else ""
subprocess.check_call([
    "conda", "install", "--yes", "-c", "conda-forge",
    "boost", "cmake", "ninja", "numpy", "pybind11", "scipy",
    f"xsimd{xsimd_version}"])
