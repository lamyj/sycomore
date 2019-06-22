from ._sycomore import *

extensions = ["como", "epg", "units"]
import sys
for extension in extensions:
    module = globals()[extension]
    name = "sycomore.{}".format(extension)

    module.__name__ = name
    sys.modules[name] = module

from . import bloch
