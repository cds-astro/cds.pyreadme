""" VizieR tool for ReadMe generation """

__version__ = "1.3.4"

from .core import CDSTablesMaker, CDSException, \
    CDSTable, CDSAstropyTable, CDSNumpyTable, CDSFileTable, CDSMRTTable, CDSAsciiTable
from .CDSColumn import CDSColumn
