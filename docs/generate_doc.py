import pathlib
import subprocess
import sys

import docutils.parsers.rst.states
import jinja2
import sphinx.ext.autodoc.directive

def main():
    here = pathlib.Path(__file__).parent
    docs = here
    
    documenter = Autodocumenter(docs)

    environment = jinja2.Environment(
        loader=jinja2.FileSystemLoader(documenter.app.srcdir))
    environment.globals["autodoc"] = documenter.autodoc
    
    for path in docs.rglob("*.rst.in"):
        template = environment.get_template(str(path.relative_to(docs)))
        contents = template.render()
        
        destination = path.parent/path.stem
        destination.write_text(contents)
    
    for path in docs.rglob("*.rstw"):
        subprocess.check_call(["pweave", "-F", ".", path.name], cwd=path.parent)

class Autodocumenter(object):
    def __init__(self, docs):
        
        self.app = sphinx.application.Sphinx(
            docs, docs, docs/"_build", docs/"_build"/"doctrees", None)

        self.environment = sphinx.environment.BuildEnvironment(self.app)
        self.reporter = sphinx.util.docutils.LoggingReporter("autodoc")

        settings = docutils.parsers.rst.states.Struct(tab_width=8)
        document = docutils.parsers.rst.states.Struct(settings=settings)
        self.state = docutils.parsers.rst.states.Struct(document=document)

    def autodoc(self, object_type, name):
        DocumenterClass = self.app.registry.documenters[object_type]
        
        options = (
            {"members": None, "undoc-members": None} 
            if object_type in ["class", "module"] else {})
        options = sphinx.ext.autodoc.directive.process_documenter_options(
            DocumenterClass, self.environment.config, options)
        
        params = sphinx.ext.autodoc.directive.DocumenterBridge(
            self.environment, self.reporter, options, 0, self.state)
        documenter = DocumenterClass(params, name)
        documenter.generate()
        return "\n".join(params.result)

if __name__ == "__main__":
    sys.exit(main())
