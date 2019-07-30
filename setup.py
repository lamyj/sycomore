import os
import re
import subprocess
import sys

import setuptools
import setuptools.command.build_ext

here = os.path.abspath(os.path.dirname(__file__))

class build_ext(setuptools.command.build_ext.build_ext):
    user_options = setuptools.command.build_ext.build_ext.user_options + [
        ("cmake=", None, "Location of CMake. Defaults of CMake in PATH")
    ]

    def __init__(self, dist):
        setuptools.command.build_ext.build_ext.__init__(self, dist)
        self.extra_cmake_options = [
            "-DCMAKE_BUILD_TYPE:STRING=Release",
            "-DBUILD_TESTING:BOOL=OFF",
            "-DBUILD_EXAMPLES:BOOL=OFF",
            "-DBUILD_STANDALONE_PYTHON_WRAPPERS:BOOL=ON",
            "-DUSE_OPENMP:BOOL=OFF",
            "-DPYTHON_EXECUTABLE:FILEPATH={}".format(sys.executable),
        ]
    
    def initialize_options(self):
        setuptools.command.build_ext.build_ext.initialize_options(self)
        self.cmake = subprocess.check_output(["which", "cmake"]).strip().decode()
    
    def finalize_options(self):
        setuptools.command.build_ext.build_ext.finalize_options(self)

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
        if not os.path.isdir(self.build_temp):
            os.makedirs(self.build_temp)
        
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
        
        self.announce("Building in {}".format(self.build_temp), 3)
        
        self.spawn(command)
    
    def _install(self, extension):
        self.announce("Installing to {}".format(self.build_lib), 3)
        
        self.spawn([
            self.cmake, "--build", self.build_temp, "--target", "install"])

with open(os.path.join(here, "CMakeLists.txt")) as fd:
    version = re.search(r"set\(sycomore_VERSION (.+?)\)", fd.read())
if not version:
    raise Exception("Could not get version from CMakeLists.txt")
version = version.group(1)

long_description = open(os.path.join(here, "README.md")).read()

if os.path.isfile(os.path.join(here, ".git")):
    sources = subprocess.check_output(
            ["git", "ls-tree", "-r", "--name-only", "HEAD"]
        ).decode().splitlines()
else:
    sources = []
    for dirpath, dirnames, filenames in os.walk(str(here)):
        sources.extend(os.path.join(dirpath, x) for x in filenames)

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
        
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    
    keywords = ["MRI", "simulation", "EPG", "Bloch"],
    
    packages=["sycomore"],
    package_dir={"sycomore": "src"},
    python_requires=">=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*,",
    ext_modules=[setuptools.Extension("_sycomore", sources)],
    cmdclass={"build_ext": build_ext},
    # WARNING: numpy >= 1.17 requires >= 3.5
    install_requires=["numpy" if sys.version >= "3.5" else "numpy<=1.16.4"],
)
