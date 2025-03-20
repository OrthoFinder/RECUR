import os 
import sys 

try:
    from importlib.metadata import PackageNotFoundError, version
    __version__ = version(__name__)
except PackageNotFoundError:
    from ._version import __version__ as __version__

if getattr(sys, 'frozen', False):
    __location__ = os.path.split(sys.executable)[0]
else:
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
