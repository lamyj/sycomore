import argparse
import glob
import os
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree

import bs4
import markdown

configuration = {}

def main():
    parser = argparse.ArgumentParser(
        description="Build and run the code samples in the documentation")
    parser.add_argument("build_directory")
    arguments = parser.parse_args()
    
    configuration = {}
    
    cache = open(os.path.join(arguments.build_directory, "CMakeCache.txt")).read()
    source = re.search(r"^sycomore_SOURCE_DIR[^=]*=(.*)$", cache, re.M).group(1)
    binary = re.search(r"^sycomore_BINARY_DIR[^=]*=(.*)$", cache, re.M).group(1)
    docs = os.path.join(source, "docs")
    configuration["include_path"] = os.path.join(source, "src")
    configuration["library_path"] = os.path.join(binary, "src")
    
    for path in sorted(glob.glob(os.path.join(docs, "*.md"))):
        check_file(path, configuration)
    
def check_file(markdown_path, configuration):
    with open(markdown_path) as fd:
        html = markdown.markdown(fd.read(), extensions=["fenced_code"])
    document = bs4.BeautifulSoup(html, "html.parser")

    languages = {"cpp": (check_cpp, ".cpp")}
    for language, (action, suffix) in languages.items():
        for code_sample in document.find_all("code", class_=language):
            fd, path = tempfile.mkstemp(suffix=suffix)
            os.write(fd, code_sample.text.encode())
            os.close(fd)
            try:
                action(path, configuration)
            except Exception:
                print(markdown_path)
                print(code_sample.text)
                raise
            os.unlink(path)

def check_cpp(path, configuration):
    subprocess.check_call([
        "c++", "-std=c++11", 
        "-I", configuration["include_path"], 
        "-L", configuration["library_path"], "-l", "sycomore",
        "-o", os.path.splitext(path)[0],
        path])
    subprocess.check_output([os.path.splitext(path)[0]])

if __name__ == "__main__":
    sys.exit(main())
