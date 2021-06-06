import multiprocessing
import os
import re
import shlex
import subprocess
import sys

import setuptools
import setuptools.command.build_ext

here = os.path.abspath(os.path.dirname(__file__))

class build_ext(setuptools.command.build_ext.build_ext):
    user_options = setuptools.command.build_ext.build_ext.user_options + [
        ("cmake=", None, "Location of CMake. Defaults to CMake in PATH"),
        ("pybind11=", None, "Location of pybind11. Defaults to autodetection"),
        ("cmake-options=", None, "Extra CMake options"),
    ]

    def __init__(self, dist):
        setuptools.command.build_ext.build_ext.__init__(self, dist)
        self.extra_cmake_options = [
            "-DCMAKE_BUILD_TYPE:STRING=Release",
            "-DBUILD_TESTING:BOOL=OFF",
            "-DBUILD_STANDALONE_PYTHON_WRAPPERS:BOOL=ON",
            "-DPYTHON_EXECUTABLE:FILEPATH={}".format(sys.executable),
        ]
    
    def initialize_options(self):
        setuptools.command.build_ext.build_ext.initialize_options(self)
        if sys.platform == "win32":
            line = subprocess.check_output(["where", "cmake"]).splitlines()[0]
            self.cmake = line.strip().decode()
        else:
            self.cmake = subprocess.check_output(["which", "cmake"]).strip().decode()
        self.pybind11 = None
        self.cmake_options = []
    
    def finalize_options(self):
        setuptools.command.build_ext.build_ext.finalize_options(self)
        
        if self.cmake is None:
            raise Exception("CMake not located in PATH")

        if not os.path.exists(self.cmake):
            raise Exception("CMake not found at {}".format(self.cmake))
        
        if self.cmake_options:
            self.cmake_options = shlex.split(self.cmake_options)
    
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
        
        command = (
            [self.cmake, "-DCMAKE_INSTALL_PREFIX={}".format(install_prefix)]
            + self.extra_cmake_options)
        if self.pybind11:
            command += ["-Dpybind11_DIR={}".format(self.pybind11)]
        command +=  self.cmake_options
        command += [str(here)]
        self.announce(" ".join(command), 3)
        
        # WARNING: we need to build in build_temp, self.spawn cannot change cwd
        subprocess.check_call(command, cwd=self.build_temp)
    
    def _build(self, extension):
        command = [
            self.cmake, "--build", self.build_temp, 
            "--parallel", str(multiprocessing.cpu_count())]
        
        self.announce("Building in {}".format(self.build_temp), 3)
        
        self.spawn(command)
    
    def _install(self, extension):
        self.announce("Installing to {}".format(self.build_lib), 3)
        
        self.spawn([
            self.cmake, "--build", self.build_temp, "--target", "install"])

with open(os.path.join(here, "CMakeLists.txt")) as fd:
    version = re.search(r"project\(\"sycomore\" VERSION (.+?)\)", fd.read())
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
    
    url="https://sycomore.readthedocs.io/",
    
    author="Julien Lamy",
    author_email="lamy@unistra.fr",
    
    license="MIT",
    
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        
        "Environment :: Console",
        
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
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
    
    setup_requires=["setuptools_scm"],
    use_scm_version=True,
    
    ext_modules=[setuptools.Extension("_sycomore", sources)],
    cmdclass={"build_ext": build_ext},
    install_requires=["numpy"],
)
