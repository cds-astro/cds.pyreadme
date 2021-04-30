"""
   G.Landais (CDS) 28 nov 2015
   Generate ReadMe and CDS standardized tables (in ASCII aligned columns)
"""

from astropy.io import ascii
from astropy.table import Table
from cdspyreadme.CDSColumn import CDSColumn
import numpy as np
import sys, os, re
from string import Template
from textwrap import wrap, fill
import datetime
import logging
import math

MAX_SIZE_README_LINE = 80
MAX_COL_INTLIMIT = 10000000


class CDSException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class CDSTable:
    """Manage table
    """

    def __init__(self, table):
        """ Constructor
        :param table: astropy table
        """
        self.logger = logging.getLogger('CDSTable')
        self.__bytebybytetemplate = None
        self.__linewidth = None
        self.__cds_columns = [] # array of CDSColumn
        self.notes = None
        self.table = table
        self.nlines = None
        self.nullvalue = None
        self.init_meta_columns()

    def get_column(self, name=None):
        """Get CDS meta columns
        :param name: column name (default is None)
        :return: CDSColumn
        """
        if name is None:
            return self.__cds_columns
        for column in self.__cds_columns:
            if column.name == name:
                return column
        return None

    def setByteByByteTemplate(self, name):
        """Set a User Byte-By-Byte template
        """
        self.__bytebybytetemplate = name

    def getByteByByteTemplate(self):
        """Get the Byte-By-Byte template
        """
        return self.__bytebybytetemplate

    def __writeTable(self, fo):
        col_length = len(self.__cds_columns)

        for nrec in range(len(self.table)):
            fo.write(self.__cds_columns[0].value(nrec))
            for i in range(1, col_length):
                fo.write(" " + self.__cds_columns[i].value(nrec))
            fo.write("\n")

    def makeCDSTable(self, fd=None):
        """Make the standardized table in ASCII aligned format.
        :param fd: file descriptor (by default, the methods creates a new file with CDSTable.name)
        """
        for col in self.__cds_columns:
            col.parse()
            if self.nullvalue:
                col.set_null_value(self.nullvalue)

        if fd is None:
            fd = open(self.name, "w") #, buffering=2048)
        self.__writeTable(fd)
        fd.close()

    def init_meta_columns(self):
        """Initialize list of CDSColumns  (self.__cds_columns)
        """
        if self.table is None:
            raise CDSException("table is empty")

        if isinstance(self.table, str): return;
        if self.__cds_columns: return self.__cds_columns

        self.__cds_columns = []
        for column in self.table.columns:
            col = CDSColumn(self.table[column])
            self.__cds_columns.append(col)

    def getlinewidth(self):
        """Get ASCII table line width
        :return: the line size
        """
        if self.__linewidth != None:
            return self.__linewidth

        self.__linewidth = 0
        for column in self.__cds_columns:
            self.__linewidth += column.size + 1
        self.__linewidth -= 1
        return self.__linewidth

    def gatherSexagesimalColumn(self):# not used
        """gather/detects sexagesimal columns if in different columns
        :return: True if found
        """
        if isinstance(self.table, Table): return

        array = ""
        for col in self.table.columns:
            if self.__column.dtype.name.startswith("i"):
                array +='i'
            elif self.__column.dtype.name.startswith("f"):
                array +='f'
            elif self.__column.dtype.name.startswith("s"):
                array += 's'
            else:
                array += ' '
        n = array.indexof("iifiif")
        if n < 9:
            return False
        return True



class CDSAstropyTable(CDSTable):
    """Manage astropy table
    """

    def __init__(self, table, name=None, description=None):
        """Constructor
        :param table: astropy table
        :param name: table name in output
        :param description: table description
        """
        if not isinstance(table, Table):
            raise CDSException("input is not Astropy Table")

        CDSTable.__init__(self, table)

        self.__filename = "astropytable"
        if name is None:
            self.name = self.__filename
        else:
            self.name = name
        self.description = description
        self.nlines = len(table)


class CDSNumpyTable(CDSTable):
    """Manage numpy table
       (TODO: problem with None values)
    """

    def __init__(self, table, name=None, description=None):
        """Constructor
        :param table: Numpy table
        :param name: table name in output
        :param description: table description
        """
        if not isinstance(table, np.ndarray):
            raise CDSException("input is not a Numpy Table")

        self.table = Table(table)
        CDSTable.__init__(self, self.table)

        self.__filename = "numpytable"
        if name is None:
            self.name = self.__filename
        else:
            self.name = name
        self.description = description
        self.nlines = len(table)


class CDSFileTable(CDSTable):
    """Manage Table in TSV/CSV or ASCII aligned column file
    """

    def __init__(self, table, name=None, description=None, data_start=None):
        """Constructor
        :param table: table file name
        :param name: table name in output
        :param description: table description
        """
        if not isinstance(table, str):
            raise CDSException("input is not a Table name")

        self.__filename = table
        if data_start:
            self.table = ascii.read(self.__filename, data_start=data_start)
        else:
            self.table = ascii.read(self.__filename) # automated search, capable to read csv header
        CDSTable.__init__(self, self.table)

        if name is None:
            self.name = self.__filename + ".cds"
        else:
            self.name = name
        self.description = description
        self.nlines = len(self.table)


class CDSAsciiTable(CDSFileTable):
    """CDS ASCII aligned file
       (long to execute because it creates a temporary csv file
        ant then it uses CDSFileTable)
    """
    def __init__(self, filename, name=None, description=None, data_start=None):
        self.__filename = filename
        self.to_sv(filename, filename+".tmp")
        CDSFileTable.__init__(self, filename+".tmp", name, description, data_start)
        #if self.gatherSexagesimalColumn():
        #    self.logger.warning("sexagesimal detected can't be well formatted by the program")

    def to_sv(self, filename, out_filename, separator=','):
        """Create  CSV table from ASCII aligned table
        :param filename: ASCII aligned filename in input
        :param out_filename: CSV filename in output
        :param separator: char separator in output (default CSV)
        """
        # memorize the full file
        with open(filename, "r") as fd:
            lines = fd.readlines()

        self.nlines = len(lines)
        self.size = 0

        # get an array filled of bool to indicate blanks
        for line in lines:
            n = len(line)
            if n > self.size: self.size = n
        filled = [False] * self.size
        for line in lines:
            for i in range(len(line)):
                if line[i] != ' ':
                    filled[i] = True

        # make the format
        col_length = []
        n = 0
        last = False
        for b in filled:
            if b is False:
                if last is False:
                    pass
                elif n > 0:
                    col_length.append(n)
                    n = 0
            n += 1
            last = b
        col_length.append(n)

        # simple test
        count = 0
        for n in col_length: count += n
        if count != self.size: raise Exception("error parsing")

        # create a new file
        with open(out_filename, "w") as out:
            for line in lines:
                cursor = 0
                buff = []
                for length in col_length:
                    buff.append(line[cursor:cursor+length])
                    cursor += length

                out.write(separator.join(buff))

    def makeCDSTable_2(self): # not used
        for col in self.get_column():
            col.parse()

        # copy table
        if self.__filename != self.name:
            with open(self.__filename, "r") as fd:
                with open(self.name, "w") as fout:
                    for line in fd:
                        fout.write(line)


class CDSMRTTable(CDSTable):
    """Manage MRT file to add information in the ReadMe
    """
    def __init__(self, tablein, tableout=None, description=None, set_limit=False):
        """Constructor
        :param tablein:  table in intput
        :param tableout: table in output
        :param description: description
        :param set_limit: add limits to the byte-by-byte description
        """
        CDSTable.__init__(self, tableout)

        self.name = tableout

        self.__bbb = None
        self.notes = []
        self.__begin_data = -1

        self.__tablename = tablein
        if tableout is None:
            self.__tablename_cds = tablein.replace("mrt", "").replace(".txt", "") + ".dat"
        self.__tablename_cds = tableout

        self.__linewidth = -1
        self.nlines = -1
        self.__columns = []

        self.description = description
        if self.description is None:
            self.description = ""

        self.__is_written = False
        self.__parseTable()

        if set_limit:
            try:
                self.__bytebybytes_statistics()
            except Exception as e:
                self.logger.error("byte-by-bytes parsing : {0}".format(str(e)))

    def __parseTable(self):
        num = 0
        linewidth = "0"
        regbbb = re.compile("^\s*[0-9]*[ -]+([0-9]*)\s+[A-Z][0-9.]+\s+.*$")
        regsep = re.compile("^[-]+$")
        regnot = re.compile("^\s*Note *[(][0-9]+[)].*")
        regtit = re.compile("^\s*Table\s*:\s(.*)\s*$")

        status = 0
        fd = open(self.__tablename, "r")
        for line in fd:
            num += 1
            if status == 5:
                continue

            if status == 0:
                mo = regtit.match(line)
                if mo:
                    if len(mo.group(1).strip()) > 1:
                        self.description = mo.group(1)
                    continue

                if line.find("Byte-by-byte Description ") == 0:
                    status = 1
                    self.__bbb = ""
                else:
                    if line.find("========") != 0:
                        self.description += " " + line.strip()

                continue

            elif status == 1:
                # continue until byte-by-bytes begining
                mo = regbbb.match(line)
                if mo:
                    status = 2
                else:
                    self.__bbb += line
                    continue

            if status == 2:
                # byte-by-bytes
                mo = regsep.match(line)
                if mo:
                    status = 3
                mo = regbbb.match(line)
                if mo:
                    linewidth = mo.group(1)
                self.__bbb += line

            elif status == 3:
                # notes -
                mo = regnot.match(line);
                if mo:
                    self.notes.append(line)
                    continue
                elif len(self.notes) == 0:
                    status = 5
                    self.__begin_data = num
                    continue

                mo = regsep.match(line)
                if mo:
                    status = 4
                    continue

                self.notes[-1] = self.notes[-1] + line

            elif status == 4:
                # data
                status = 5
                self.__begin_data = num
                continue

        fd.close()

        # print (self.__bbb)
        self.nlines = int(num) - int(self.__begin_data) + 1
        self.__linewidth = int(linewidth)

    def __bytebybytes_statistics(self):
        self.table = ascii.read(self.__tablename)
        self.logger.debug("initmetacolumn from bbb\n")

        self.init_meta_columns()
        columns = self.get_column()
        for col in columns:
            col.parse()
            
        regbbb = re.compile("^(\s*[\d \-]*\d\s+[A-Z][0-9.]*\s+[^\s]+\s+[^\s]+\s+)(.*$)$")

        ncol = 0
        bbb_out = []
        bbb_in = self.__bbb.split("\n")
        for line in bbb_in:
            line = line.replace('?=""', '?')

            mo = regbbb.match(line)
            if mo is None:
                bbb_out.append(line)
                continue


            if ncol >= len(columns):
                self.logger.warning("byte-by-byte limit can't be interpreted")
                return

            if columns[ncol].formatter.min is not None and columns[ncol].formatter.max is not None:
                s = mo.group(1) + "[" + str(columns[ncol].formatter.min) + "/" + str(columns[ncol].formatter.max) + "]"
                if mo.group(2)[0] != '?': s += " "
                s += mo.group(2).lstrip()
                bbb_out.append(s)
            else:
                bbb_out.append(line)
            ncol += 1

        if ncol == len(columns):
            self.__bbb = "\n".join(bbb_out) + "\n"

    def __writeTable(self):
        i = 0
        fd = open(self.__tablename, "r")
        fout = open(self.__tablename_cds, "w")

        for line in fd:
            i += 1
            if i < self.__begin_data:
                continue
            fout.write(line)

        fd.close()
        fout.close()
        self.__is_written = True

    def makeCDSTable(self):
        """Make the standardized table in
           ASCII format with aligned columns.
        """
        if self.__is_written: return
        self.__writeTable()

    def setByteByByteTemplate(self, name):
        """Set a User Byte-By-Byte template
        """
        raise Exception("not available in MRT")

    def getByteByByteTemplate(self):
        """Get the Byte-By-Byte template
        """
        return os.path.dirname(__file__) + "/bytebybyte.template"

    def getByteByByte(self):
        """Get byte-by-byte.
        :return: byte-by-byte string
        """
        return self.__bbb

    def getlinewidth(self):
        """get ASCII table line width
        :return: ASCII line size
        """
        return self.__linewidth

    def get_columns(self):
        """get CDS meta columns
        :return: CDSColumn array
        """
        return self.__columns


class CDSTablesMaker:
    """Generate standardized tables and ReadMe
    """

    def __init__(self, out=None, debug=False):
        """Constructor
        :param out: the output file (default: stdout)
        :param debug: True/False
        """
        self.__dir = os.path.dirname(os.path.realpath(__file__))
        self.__tables = []
        self.__ReadMetemplate = self.__dir + "/ReadMe.template"

        if out != None:
            sys.stdout = open(out, 'w')

        self.title = 'Title ?'
        self.author = '1st author ?'
        self.catalogue = ''
        self.date = 'Date ?'
        self.abstract = 'Abstract ?'
        self.authors = 'Authors ?'
        self.bibcode = 'ref ?'
        self.keywords = ''
        self.ref = None
        self.__templatevalue = None
        self.logger = logging.getLogger('CDSTablesMaker')

        if debug is True:
            self.logger.basicConfig(level=logging.DEBUG)

    def addTable(self, table, name=None, description=None, nullvalue=None):
        """Add a Table, memorize the meta-data and generate the standardized CDS table
        :param table: table (type accepted: astropy, numpy, filename, CDSTable)
        :param name: the name used in output
        :param description: table description
        :param nullvalue: set a nullvalue (aplyed for all columns)
        :return: CDStable created
        """
        self.logger.debug("add table")
        if isinstance(table, CDSTable):
            cdstable = table
            if name != None:
                cdstable.name = name
            if description != None:
                cdstable.description = description
        elif isinstance(table, Table):
            cdstable = CDSAstropyTable(table, name, description)
        elif isinstance(table, str):
            cdstable = CDSFileTable(table, name, description)
        elif isinstance(table, np.ndarray):
            cdstable = CDSNumpyTable(table, name, description)
        else:
            raise CDSException("type " + type(table) + " is not accepted (only String or astropy.Table)")

        self.__tables.append(cdstable)
        self.logger.debug("append table ok")

        cdstable.nullvalue = nullvalue
        return cdstable

    def writeCDSTables(self):
        """Write tables in ASCII format
        """
        for table in self.__tables:
            table.makeCDSTable()
            self.logger.debug("make CDS table " + table.name)

    def getTablesInfo(self):
        """get tables information.
        :return: info (dictionary)
        """
        info = []
        for table in self.__tables:
            info.append({'name': table.name, 'lines': len(table.table), 'width': table.getlinewidth()})
        return info

    def getTables(self):
        """Get the CDSTable list
        :return: list of CDSTable
        """
        return self.__tables

    def printTablesIndex(self, outBuffer=False):
        """Print the tables index
        :param outBuffer: true to get buffer, else write on output (default: False)
        """
        sz = [14, 0, 8]
        for tab in self.__tables:
            if len(tab.name) > sz[0]:
                sz[0] = len(tab.name)
            l = len(str(tab.getlinewidth()))
            if l > sz[1]:
                sz[1] = l

        fmtt = "{0:" + str(sz[0]) + "s} {1:" + str(sz[1]) + "d} {2:>" + str(sz[2]) + "s}  {3:s}"
        buff = fmtt.format("ReadMe", MAX_SIZE_README_LINE, ".", "this file") + "\n"
        for tab in self.__tables:
            buff += fmtt.format(self.__strFmt(tab.name),
                                tab.getlinewidth(),
                                str(tab.nlines),
                                self.__strFmt(tab.description))
            buff += "\n"

        if outBuffer: return buff
        sys.stdout.write(buff)

    def __strFmt(self, string):
        if string is None:
            return ""
        else:
            return string

    def __strFmtRa(self, column, fmtb, startb):
        ra = column.getSexaRA()
        n = startb
        nvalue = ""
        if column.hasNull: nvalue = "? "

        buff = fmtb.format(n, n + ra.RAh.size - 1, "",
                           ra.RAh.fortran_format, ra.RAh.unit, ra.RAh.name, nvalue + ra.RAh.description)
        n += ra.RAh.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + ra.RAm.size - 1, "",
                            ra.RAm.fortran_format, ra.RAm.unit, ra.RAm.name, nvalue + ra.RAm.description)
        n += ra.RAm.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + ra.RAs.size - 1, "",
                            ra.RAs.fortran_format, ra.RAs.unit, ra.RAs.name, nvalue + ra.RAs.description)
        return buff

    def __strFmtDe(self, column, fmtb, startb):
        de = column.getSexaDE()
        n = startb
        nvalue = ""
        if column.hasNull: nvalue = "? "
        buff = fmtb.format(n, n + de.DEsign.size - 1, "",
                           de.DEsign.fortran_format, de.DEsign.unit, de.DEsign.name, nvalue + de.DEsign.description)
        n += de.DEsign.size
        buff += '\n'
        buff += fmtb.format(n, n + de.DEd.size - 1, "",
                            de.DEd.fortran_format, de.DEd.unit, de.DEd.name, nvalue + de.DEd.description)
        n += de.DEd.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + de.DEm.size - 1, "",
                            de.DEm.fortran_format, de.DEm.unit, de.DEm.name, nvalue + de.DEm.description)
        n += de.DEm.size + 1
        buff += '\n'
        buff += fmtb.format(n, n + de.DEs.size - 1, "",
                            de.DEs.fortran_format, de.DEs.unit, de.DEs.name, nvalue + de.DEs.description)
        return buff

    def __splitLine(self, line, shift=0):
        """Split line 80 char
           :param shift: add left blank
        """
        if shift > MAX_SIZE_README_LINE:
            shift = 0
        return ("\n" + " " * shift).join(wrap(line, width=MAX_SIZE_README_LINE - shift))

    def __add_authors(self, line, shift=0):
        """Split the line containing the authors without separate given and surname
        :param line: authors list in a line
        :param shift: add left blank
        :return: authors formatted string
        """
        # Find all spaces followed by authors's initials
        space_letter_list = re.findall(" (?:[A-Z]\.-?)+", line)

        # Replace founded spaces by !
        for space_letter in space_letter_list:
            line = line.replace(space_letter, "!" + space_letter.strip())

        # Wrap the text by using spaces as breakpoint and then replace ! by spaces so given name and surname
        # are not separate
        new_line = fill(line, width=MAX_SIZE_README_LINE, subsequent_indent=shift * " ").replace("!", " ")

        if self.author:
            if len(new_line)>0:
                return self.author+", "+new_line
            else:
                return self.author
        return new_line

    def __add_keywords(self, line, shift=0):
        """Split the line containing the authors without separate given and surname
        :param line: keywords list in a line
        :param shift: add left blank
        :return: keywods formatted string
        """
        # Replace all spaces that are NOT precede by ; with !
        line = re.sub("(?<!;) ", "!", line)

        # Wrap the text by using spaces as breakpoint and then replace ! by spaces so there is no break line
        # between two ;
        new_line = fill(line, width=MAX_SIZE_README_LINE, subsequent_indent=shift * " ").replace("!", " ")

        return new_line

    def printByteByByte(self, table, outBuffer=False):
        """Print byte-by-byte
        :param table: the CDSTable
        :param outBuffer: true to get buffer, else write on output (default: False)
        """
        if isinstance(table, CDSMRTTable):
            buff = table.getByteByByte()
            if table.notes != None:
                buff += "-" * 80 + "\n"
                for line in table.notes:
                    try:
                        buff += line.strip().encode('ascii', 'replace').decode()
                    except Exception as err:
                        self.logger.error("error detected in Notes " + str(err))
                        buff += "?" * len(line) + "\n"
                    buff += "\n"
                buff += "-" * 80 + "\n"

            if outBuffer: return buff
            sys.stdout.write(buff)
            return

        columns = table.get_column()
        startb = 1
        sz = [0, 0, 1, 7]
        l = len(str(table.getlinewidth()))
        if l > sz[0]:
            sz[0] = l
            sz[1] = l
        self.logger.debug("size sz=" + str(sz))
        fmtb = "{0:" + str(sz[0]) + "d}-{1:" + str(sz[1]) + "d} {2:" + str(sz[2]) + "s}"
        for column in columns:
            if len(column.name) > sz[3]:
                sz[3] = len(column.name)
        fmtb += " {3:6s} {4:6s} {5:" + str(sz[3]) + "s} {6:s}"
        buff = ""
        nsplit = sz[0] + sz[1] + sz[2] + sz[3] + 16

        for column in columns:
            endb = column.size + startb - 1
            if column.formatter.fortran_format[0] == 'R':
                buff += self.__strFmtRa(column, fmtb, startb) + "\n"
            elif column.formatter.fortran_format[0] == 'D':
                buff += self.__strFmtDe(column, fmtb, startb) + "\n"
            else:
                description = column.description
                if column.hasNull:
                    nullflag = "?"
                else:
                    nullflag = ""

                borne = ""
                if column.min and column.max:
                    if column.formatter.fortran_format[0] == 'I':
                        if abs(column.min) < MAX_COL_INTLIMIT and abs(column.max) < MAX_COL_INTLIMIT:
                            if column.min == column.max:
                                borne = "[{0}]".format(column.min)
                            else:
                                borne = "[{0}/{1}]".format(column.min, column.max)
                    elif column.formatter.fortran_format[0] in ('E','F'):
                        borne = "[{0}/{1}]".format(math.floor(column.min*100)/100.,
                                                   math.ceil(column.max*100)/100.)

                description = "{0}{1} {2}".format(borne, nullflag, description)
                newline = fmtb.format(startb, endb, "",
                                      self.__strFmt(column.formatter.fortran_format),
                                      self.__strFmt(column.unit),
                                      self.__strFmt(column.name),
                                      description)

                if len(newline) > MAX_SIZE_README_LINE:
                    buff += ("\n").join(wrap(newline,
                                             subsequent_indent=" " * nsplit,
                                             width=MAX_SIZE_README_LINE))
                    buff += "\n"
                else:
                    buff += newline + "\n"
            startb = endb + 2

        if table.notes != None:
            buff += "-" * 80 + "\n"
            for line in table.notes: buff += line + "\n"
            buff += "-" * 80 + "\n"

        if outBuffer: return buff
        sys.stdout.write(buff)

    def __getByteByByteTemplate(self, table):
        templateValue = {'file': table.name,
                         'bytebybyte': self.printByteByByte(table, outBuffer=True)}

        templatename = table.getByteByByteTemplate()
        if templatename is None: templatename = self.__dir + "/bytebybyte.template"
        with open(templatename) as filein:
            src = Template(filein.read())
            return src.substitute(templateValue)

    def setReadmeTemplate(self, templatename, templateValue=None):
        """Set a user ReadMe template
        :param templateName: the template name
        :param templateValue: dictionary to fill added variable in templatename
        """
        self.__ReadMetemplate = templatename
        self.__templatevalue = templateValue

    def putRef(self, catname, title=""):
        """Put a reference.
        :param catname: catalgogue name (string)
        :param title: the title (string)
        """
        if self.ref is None: self.ref = []
        if catname is None: raise Exception("catname is required")
        self.ref.append((catname, title))

    def printRef(self, outBuffer):
        """The "See also" section in ReadMe
        :param outBuffer: true to get buffer, else write on output (default: False)
        """
        if self.ref is None or len(self.ref) == 0: return

        buf = ""
        # Find the highest string length in first references column
        max_len = len(max([i[0] for i in self.ref], key=len))

        # Write references and align the double dot symbols
        for ref in self.ref:
            buf += self.__splitLine(" {0:<{max}} : {1}".format(ref[0], ref[1], max=max_len)) + "\n"

        if outBuffer is True: return buf
        sys.stdout.write(buf)

    def makeReadMe(self, out=sys.stdout):
        """Print the ReadMe
        :param out: file descriptor (default sys.stdout)
        """
        templateValue = {'catalogue': self.catalogue,
                         'title': self.__splitLine(self.title),
                         'author': self.author,
                         'date': self.date,
                         'abstract': self.__splitLine(self.abstract, shift=2),
                         'authors': self.__add_authors(self.authors, shift=4),
                         'bibcode': "=" + self.bibcode,
                         'keywords': self.__add_keywords(self.keywords, shift=len("keywords: ")),
                         'tablesIndex': self.printTablesIndex(outBuffer=True),
                         'seealso': self.printRef(outBuffer=True),
                         'bytebybyte': '',
                         'today': datetime.datetime.now().strftime("%d-%b-%Y")}

        if self.__templatevalue is not None:
            for key in self.__templatevalue: templateValue[key] = self.__templatevalue[key]

        buff = ""
        for table in self.__tables:
            buff += self.__getByteByByteTemplate(table)
            buff += "\n"
        templateValue['bytebybyte'] = buff

        with open(self.__ReadMetemplate) as filein:
            src = Template(filein.read())
            result = src.substitute(templateValue)
            out.write(result)

    def toMRT(self):
        """Transform tables into MRT (ASCII aligned table with byte-by-byte header)
        """
        templateValue = {'title': self.__splitLine(self.title),
                         'authors': self.__add_authors(self.authors, shift=4)}
        mrt_template = self.__dir + "/MRT.template"

        for table in self.__tables:
            for col in table.get_column():
                col.parse()
            templateValue['bytebybyte'] = self.__getByteByByteTemplate(table)

            with open(table.name, "w") as fd:
                with open(mrt_template) as filein:
                    src = Template(filein.read())
                    result = src.substitute(templateValue)
                    fd.write(result)
                table.makeCDSTable(fd)
