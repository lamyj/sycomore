import argparse
import glob
import os
import pickle
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree

import docutils.nodes

configuration = {}

def main():
    parser = argparse.ArgumentParser(
        description="Build and run the code samples in the documentation")
    parser.add_argument("build_directory")
    parser.add_argument("files", nargs="*", metavar="file")
    arguments = parser.parse_args()
    
    configuration = {}
    
    cache = open(os.path.join(arguments.build_directory, "CMakeCache.txt")).read()
    source = re.search(r"^sycomore_SOURCE_DIR[^=]*=(.*)$", cache, re.M).group(1)
    binary = re.search(r"^sycomore_BINARY_DIR[^=]*=(.*)$", cache, re.M).group(1)
    docs = os.path.join(source, "docs")
    configuration["include_path"] = os.path.join(source, "src")
    configuration["library_path"] = os.path.join(binary, "src")
    
    configuration["python"] = re.search(r"^PYTHON_EXECUTABLE[^=]*=(.*)$", cache, re.M).group(1)
    configuration["install"] = re.search(r"^CMAKE_INSTALL_PREFIX[^=]*=(.*)$", cache, re.M).group(1)
    
    if not arguments.files:
        arguments.files = []
        for dirpath, dirnames, filenames in os.walk(os.path.join(docs, "_build", "doctrees")):
            arguments.files.extend(
                os.path.join(dirpath, x) for x in filenames if x.endswith(".doctree"))
    
    for path in sorted(arguments.files):
        check_file(path, configuration)

class CodeVisitor(docutils.nodes.NodeVisitor):

    def visit_literal_block(self, node):
        """Called for "literal_block" nodes."""
        print(node, node.attributes.get("language"), node.astext())

    def unknown_visit(self, node):
        """Called for all other node types."""
        pass

def check_file(doctree, configuration):
    
    with open(doctree, "rb") as fd:
        document = pickle.load(fd)
    
    blocks = document.traverse(
        lambda x: isinstance(x, docutils.nodes.literal_block), descend=True)
    
    languages = {
        "cpp": (check_cpp, ".cpp"),
        "python": (check_python, ".py")
    }
    for block in blocks:
        if block["language"] not in languages:
            print("Skipping '{}' block in {}".format(block["language"], doctree))
            continue
        
        action, suffix = languages[block["language"]]
        
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.write(fd, block.astext().encode())
        os.close(fd)
        
        try:
            action(path, configuration)
        except Exception:
            print(doctree)
            print(block.astext())
            raise
        os.unlink(path)

def check_cpp(path, configuration):
    subprocess.check_call([
        "c++", "-std=c++11", 
        "-I", configuration["include_path"], 
        "-L", configuration["library_path"], "-l", "sycomore",
        "-o", os.path.splitext(path)[0],
        path])
    environment = os.environ.copy()
    environment["LD_LIBRARY_PATH"] = configuration["library_path"]
    subprocess.check_output([os.path.splitext(path)[0]], env=environment)

def check_python(path, configuration):
    python_prefix = subprocess.check_output([
        configuration["python"], "-c", 
        "from distutils.sysconfig import *; print(get_python_lib(True, prefix=''))"]).strip()
    
    environment = os.environ.copy()
    environment["LD_LIBRARY_PATH"] = os.path.join(configuration["install"], "lib")
    environment["PYTHONPATH"] = os.path.join(configuration["install"], python_prefix.decode())
    subprocess.check_call([configuration["python"], path], env=environment)

if __name__ == "__main__":
    sys.exit(main())
