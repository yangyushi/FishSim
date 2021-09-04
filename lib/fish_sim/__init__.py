from . import utility
from . import model
try:
    from . import cmodel
except ImportError:
    import warnings
    warnings.warn(
        "failed to import the cmodel module, if you don't need \
        the fast model in written C++ you can ignore this warning.",
        ImportWarning
    )
