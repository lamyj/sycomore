# https://github.com/pypa/python-manylinux-demo/blob/master/setup.py

import os
import pathlib
import re
import shutil
import subprocess
import sys

import setuptools
import setuptools.command.build_ext

here = pathlib.Path(__file__).parent.absolute()

class build_ext(setuptools.command.build_ext.build_ext):
    user_options = setuptools.command.build_ext.build_ext.user_options + [
        ("cmake=", None, "Location of CMake. Defaults of CMake in PATH")
    ]

    def __init__(self, dist):
        super().__init__(dist)
        self.extra_cmake_options = [
            "-DCMAKE_BUILD_TYPE:STRING=Release",
            "-DBUILD_TESTING:BOOL=OFF",
            "-DBUILD_EXAMPLES:BOOL=OFF",
            "-DBUILD_STANDALONE_PYTHON_WRAPPERS:BOOL=ON",
            "-DUSE_OPENMP:BOOL=OFF",
            "-DPYTHON_EXECUTABLE:FILEPATH={}".format(sys.executable),
        ]
    
    def initialize_options(self):
        super().initialize_options()
        self.cmake = shutil.which("cmake")
    
    def finalize_options(self):
        super().finalize_options()

        if self.cmake is None:
            raise Exception("CMake not located in PATH")

        if not os.path.exists(self.cmake):
            raise Exception("CMake not found at {}".format(self.cmake))
    
    def run(self):
        for extension in self.extensions:
            self.build_extension(extension)
        
    def build_extension(self, ext):
        self._configure(ext)
        self._build(ext)
        self._install(ext)
    
    def _configure(self, extension):
        os.makedirs(self.build_temp, exist_ok=True)
        
        self.announce("Running CMake in {}".format(self.build_temp), 3)
        
        # WARNING: install_prefix needs to be absolute
        install_prefix = os.path.abspath(self.build_lib)
        # WARNING: we need to build in build_temp, self.spawn cannot change cwd
        subprocess.check_call(
            [self.cmake, "-DCMAKE_INSTALL_PREFIX={}".format(install_prefix)]
            + self.extra_cmake_options
            + [str(here)],
            cwd=self.build_temp)
    
    def _build(self, extension):
        command = [self.cmake, "--build", self.build_temp]
        if self.parallel:
            command.extend(["-j", str(self.parallel)])
        
        self.announce("Building in {}".format(self.build_temp), 3)
        
        self.spawn(command)
    
    def _install(self, extension):
        self.announce("Installing to {}".format(self.build_lib), 3)
        
        self.spawn([
            self.cmake, "--build", self.build_temp, "--target", "install"])

version = re.search(
    r"set\(sycomore_VERSION (.+?)\)", (here/"CMakeLists.txt").read_text())
if not version:
    raise Exception("Could not get version from CMakeLists.txt")
version = version.group(1)

long_description = (here/"README.md").read_text()

sources = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", "HEAD"]
    ).decode().splitlines()

setuptools.setup(
    name="sycomore",
    version=version,
    
    description="MRI simulation toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    url="https://iris.icube.unistra.fr/lamy/sycomore",
    
    author="Julien Lamy",
    author_email="lamy@unistra.fr",
    
    license="MIT",
    
    classifiers=[
        "Development Status :: 4 - Beta",
        
        "Environment :: Console",
        
        "Intended Audience :: Science/Research",
        
        "License :: OSI Approved :: MIT License",
        
        "Programming Language :: Python :: 3",
        
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    
    keywords = ["MRI", "simulation", "EPG", "Bloch"],
    
    packages=["sycomore"],
    package_dir={"sycomore": "src"},
    python_requires=">=3.5",
    ext_modules=[setuptools.Extension("_sycomore", sources)],
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy"],
)
