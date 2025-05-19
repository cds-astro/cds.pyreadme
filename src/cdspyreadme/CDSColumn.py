"""
   G.Landais (CDS) 28 nov 2015
   CDSColumn, used to Generate ReadMe and CDS standardized tables (in ASCII aligned columns)
"""
import numpy
import re
import logging
from astropy.table import Column, MaskedColumn
import math
UNDEFINED_UNIT = "---"
DEFAULT_STRING_SIZE = 50


class SexaMeta:
    """Sexagesimal metadata"""

    def __init__(self, name: str, description: str, fortran_format: str, unit: str):
        """Constructor.
        :param name: column name
        :param description: column description
        :param fortran_format: Fortran format
        :param unit: column unit
        """
        self.name = name
        self.description = description
        self.fortran_format = fortran_format
        self.unit = unit
        self.size = int(fortran_format[1])


class SexaRA:
    def __init__(self, number_sec_digits: int = 0):
        """Constructor.
        :param number_sec_digits: number of digit in seconds
        """
        self.RAh = SexaMeta("RAh", "Right ascension (hour)", "I2", "h")
        self.RAm = SexaMeta("RAm", "Right ascension (minute)", "I2", "min")
        self.RAs = SexaMeta("RAs", "Right ascension (seconds)", "I2", "s")
        if number_sec_digits != 0:
            self.RAs.fortran_format = "F" + str(3 + number_sec_digits) + "." + str(number_sec_digits)
            self.RAs.size = (3 + number_sec_digits)


class SexaDE:
    def __init__(self, number_sec_digits: int):
        """Constructor.
        :param number_sec_digits: number of digit in seconds
        """
        self.DEsign = SexaMeta("DE-", "Declination (degree)", "A1", "---")
        self.DEd = SexaMeta("DEd", "Declination (degree)", "I2", "deg")
        self.DEm = SexaMeta("DEm", "Declination (minute)", "I2", "arcmin")
        self.DEs = SexaMeta("DEs", "Declination (seconds)", "I2", "arcsec")
        if number_sec_digits != 0:
            self.DEs.fortran_format = "F" + str(3 + number_sec_digits) + "." + str(number_sec_digits)
            self.DEs.size = (3 + number_sec_digits)


class CDSColumnFormatter:
    """CDS Column formatter interface
       with formats, size, min, max
    """
    def __init__(self):
        self.fortran_format = None
        self.format = None
        self.size = None
        self.min = None
        self.max = None
        self.out_format = None
        self.none_format = None

    def write(self, value: any):
        """write a value using the formatter
        :param value: value in input
        :return: the formatted value (for ASCII aligned serialization)
        """
        if isinstance(value, numpy.ma.core.MaskedConstant):
            return self.none_format.format("")
        return self.out_format.format(value)


class CDSColumnIntegerFormatter(CDSColumnFormatter):
    """ CDS Column  for integer
    """
    def __init__(self, column: Column, has_null: bool = False):
        """Constructor - build Format from an astropy column
        :param column: astropy column
        :param has_null: contain null values
        """
        CDSColumnFormatter.__init__(self)
        self.__set_slot(column, has_null)
        self.size = len(str(self.max))
        l = len(str(self.min))
        if self.size < l:
            self.size = l
        self.fortran_format = "I" + str(self.size)
        self.format = ">" + self.fortran_format[1:]

        self.out_format = "{0:" + self.format + "}"
        self.none_format = "{0:" + str(self.size) + "s}"

    def __set_slot(self, column: Column, has_null: bool):
        if has_null:  # optimized (2x more speed)
            mask_col = column  # MaskedColumn(column, mask=[col is None for col in column])
            mask_col.fill_value = -999999
            self.max = max(mask_col.filled())
            if self.max == -999999:
                self.max = None
            mask_col.fill_value = +999999
            self.min = min(mask_col.filled())
            if self.min == 999999:
                self.min = None
        else:
            self.max = max(column)
            self.min = min(column)


class CDSColumnStringFormatter(CDSColumnFormatter):
    """CDS column for String
    """
    def __init__(self, column: Column, has_null: bool):
        """Constructor - build Format from an astropy column
        :param column: astropy column
        :param has_null: contain null values
        """
        CDSColumnFormatter.__init__(self)

        try:
            if has_null:
                mask_col = column  # MaskedColumn(column, mask=[col is None for col in column])
                mask_col.fill_value = ""
                col_tmp = Column(mask_col.filled(), dtype=str)
                self.size = int(re.sub(r'^[^0-9]+(\d+)$', r'\1', col_tmp.dtype.str))
            else:
                self.size = int(re.sub(r'^[^0-9]+(\d+)$', r'\1', column.dtype.str))
        except Exception as e:
            logging.error(e)
            self.size = DEFAULT_STRING_SIZE

        self.fortran_format = "A" + str(self.size)
        self.format = str(self.size) + "s"
        self.out_format = "{0:" + str(self.size) + "s}"
        self.none_format = self.out_format


class FloatRepr:
    """Float representation
       extract precision, width, format....
       The class can be used for a single value or a list of values (method consume()) to get the format
       of a single or a list
    """
    SCI = 0
    DEC = 1
    MIXED = 2

    REG_FLOAT = re.compile("([+-]*)([^eE.]+)([.]*)([0-9]*)([eE]*-*)[0-9]*")

    def __init__(self):
        """Constructor.
        """
        self.type = -1
        self.sign = False  # contains sign

        self.e_width = -1  # width (type SCI)
        self.e_precision = -1  # precision (type SCI)
        self.f_coma = -1  # decimal, after coma (type DEC)
        self.f_dec = -1  # decimal, before coma (type DEC)

    def consume(self, value: str):
        """set format from a value
        :param value: (string) value to consume
        """
        mo = FloatRepr.REG_FLOAT.match(value)
        if mo is None: 
            raise Exception(f"{value} is not a float number")

        if self.sign is False:
            self.sign = True if mo.group(1) != "" else False

        if mo.group(5) != "":  # Scientific format
            e_precision = len(mo.group(2))+len(mo.group(4))
            if self.e_precision < e_precision:
                self.e_precision = e_precision

            w = len(value)
            if self.e_width < w:
                self.e_width = w

            if self.type in (FloatRepr.MIXED, FloatRepr.DEC):
                self.type = FloatRepr.MIXED
            else:
                self.type = FloatRepr.SCI

        else:  # Decimal format
            coma = len(mo.group(4))
            if self.f_coma < coma:
                self.f_coma = coma

            dec = len(mo.group(2))
            if self.f_dec < dec:
                self.f_dec = dec

            if self.type in (FloatRepr.MIXED, FloatRepr.SCI):
                self.type = FloatRepr.MIXED
            else:
                self.type = FloatRepr.DEC

    def width(self) -> int:
        """get size (when serialized in string)
        :return: the size
        """
        if self.type == FloatRepr.SCI:
            return self.e_width
        w = self.f_dec + self.f_coma + 1
        if self.sign:
            w += 1
        if self.type == FloatRepr.DEC:
            return w
        return w if w > self.e_width else self.e_width

    def precision(self) -> int:
        """get precision
        :return: precision
        """
        if self.type == FloatRepr.SCI:
            return self.e_precision
        w = self.f_dec + self.f_coma
        if self.type == FloatRepr.DEC:
            return w
        return w if w > self.e_precision else self.e_precision

    def fortran_format(self) -> str:
        """return FORTRAN format style
        """
        if self.type == FloatRepr.SCI:
            return f"E{self.width()}.{self.ePrecision}"
        elif self.type == FloatRepr.MIXED:
            return f"E{self.width()}.{self.precision()}"
        else:
            return f"F{self.width()}.{self.f_coma}"

    def c_format(self) -> str:
        """return C Format style
        """
        if self.type == FloatRepr.SCI:
            return f"{self.width()}.{self.ePrecision}e"
        elif self.type == FloatRepr.MIXED:
            return f"{self.width()}.{self.precision()+1}g"
        else:
            return f"{self.width()}.{self.f_coma}f"


class CDSColumnFloatFormatter(CDSColumnFormatter):
    """CDS column for float
    """
    def __init__(self, column: Column, has_null: bool):
        """Constructor - build Format from an astropy column
        :param column: astropy column
        :param has_null: contain null values

        Warning: compare float that can use (or not or sometimes)
                 scientific format (C-type %e or %f or %g)

        Ex: -12.12345     -> F10.5  (%10.5f)
            -1.212345e-1  -> E12.7
        """
        CDSColumnFormatter.__init__(self)
        self.__set_slot(column, has_null)
        self.__parse(column)

    def __parse(self, column: Column):
        fmt_saved = FloatRepr()

        for rec in column:
            # skip null values
            if rec is None:
                continue
            s = str(rec)
            fmt_saved.consume(s)

        self.size = fmt_saved.width()
        self.fortran_format = fmt_saved.fortran_format()
        self.format = fmt_saved.c_format()

        self.out_format = "{0:" + self.format + "}"
        self.none_format = "{0:" + str(self.size)+"s}"

    def __set_slot(self, column: Column, has_null: bool):
        if has_null:  # optimized (2x more speed)
            masked_col = column  # MaskedColumn(column, mask=[col is None for col in column])
            masked_col.fill_value = -999999
            self.max = max(masked_col.filled())
            if self.max == -999999:
                self.max = None
            masked_col.fill_value = +999999
            self.min = min(masked_col.filled())
            if self.min == 999999:
                self.min = None
        else:
            self.max = max(column)
            self.min = min(column)

    def write(self, value: any):
        """write a value using the formatter
        :param value: value in input
        :return: the formatted value (for ASCII aligned serialization)
        """
        if isinstance(value, numpy.ma.core.MaskedConstant):
            return self.none_format.format("")
        if math.isnan(value):
            return self.none_format.format("")
        return self.out_format.format(value)


#class CDSColumnSexaRAFormatter(CDSColumnFormatter):
#    """CDS column for RightAscension in sexagesimal
#    """
#    def __init__(self, formatter):
#        CDSColumnFormatter.__init__(self)
#        self.size = formatter.size
#        self.ffmt = formatter.fortran_format.replace('A', 'R')
#        self.format = "s"
#        self.out_format = "{0:" + self.format + "}"
#
#    def write(self, value):
#        if not isinstance(value, numpy.string_):
#            return self.out_format.format(" ")
#        return self.out_format.format(str(value).replace(":", " "))


#class CDSColumnSexaDEFormatter(CDSColumnFormatter):
#    """CDS column forDeclination in sexagesimal
#    """
#    def __init__(self, formatter):
#        CDSColumnFormatter.__init__(self)
#        self.size = formatter.size
#        self.fortran_format = formatter.fortran_format.replace('A', 'D')
#        self.format = "s"
#        self.out_format = "{0:" + self.format + "}"
#
#    def write(self, value):
#        if not isinstance(value, numpy.string_):
#            return self.out_format.format(" ")
#        return self.out_format.format(str(value).replace(":", " "))


class CDSColumn:
    """CDS Column decorator on astropy.table.Column
    """
    def __init__(self, column: Column):
        """Constructor.
        :param column: astropy Column
        """
        if not isinstance(column, Column):
            raise Exception("column must be astropy.table.Column")

        self.name = column.name
        self.size = None
        self.min = None
        self.max = None
        self.hasNull = None
        self.description = "Description of " + self.name
        self.unit = UNDEFINED_UNIT
        self.formatter = None

        self.__dbname = column.name
        self.__column = column
        self.__sexagesimal = [None, None]
        self.__force_format = None

    def set_format(self, fmt: str):
        """force a format (fortran format)
        :param fmt: the new format (ex:  %F10.6)
        """
        formatter = CDSColumnFormatter()
        formatter.fortran_format = fmt

        mo = re.match("^([EF])([0-9]+)[.]([0-9]+)", fmt)
        if mo:
            if mo.group(1) == 'E':
                f = "{0}.{1}e".format(mo.group(2), mo.group(3))
                formatter.size = int(mo.group(2))
                formatter.format = "%"+f
                formatter.out_format = "{{0:{}.{}e}}".format(mo.group(2), int(mo.group(3))-1)
            else:
                f = "{0}.{1}f".format(mo.group(2), mo.group(3))
                formatter.size = int(mo.group(2))
                formatter.format = "%"+f
                formatter.out_format = "{0:"+f+"}"

        else:
            mo = re.match("^([IA])([0-9]+)", fmt)
            if mo is None:
                raise Exception("format {} not accepted".format(fmt))

            if mo.group(1) == "I":
                f = "{}d".format(mo.group(2))
                formatter.size = int(mo.group(2))
                formatter.format = "%"+f
                formatter.out_format = "{0:>"+mo.group(2)+"}"
            else:
                f = "{}s".format(mo.group(2))
                formatter.size = int(mo.group(2))
                formatter.format = "%"+f
                formatter.out_format = "{0:"+f+"}"
 
        formatter.none_format = "{0:"+str(formatter.size)+"}"
        self.__force_format = formatter

    def set_null_value(self, null_value: any):
        """Assign null value to the Column
        (create an astropy  MaskedColumn)
        :param null_value: value
        """
        mask = []
        for value in self.__column:
            if isinstance(value, numpy.ma.core.MaskedConstant):
                mask.append(True)
            else:
                mask.append((value == null_value))
        # mask = [(value==null_value) for value in col]
        self.__column = MaskedColumn(self.__column, mask=mask)

    def parse(self):
        """the method parse the columns and set type, size format
        """
        if self.formatter is not None:
            return

        if self.__column.description:
            self.description = self.__column.description

        if self.unit is None or self.unit == UNDEFINED_UNIT:
            if self.__column.unit is not None:
                try:
                    self.unit = self.__column.unit.to_string("cds")
                except:
                    self.unit = self.__column.unit
            else:
                self.unit = self.__get_unit()

        if isinstance(self.__column, MaskedColumn) is False:
            for col in self.__column:
                if col is None:
                    logging.debug("detect None value "+str(col))
                    self.__column = MaskedColumn(self.__column, mask=[col is None for col in self.__column])
                    break

        self.hasNull = self.__has_null()
        type = self.__get_type()
        if type is int:
            self.formatter = CDSColumnIntegerFormatter(self.__column, self.hasNull)
        elif type is float:
            self.formatter = CDSColumnFloatFormatter(self.__column, self.hasNull)
        else:
            self.formatter = CDSColumnStringFormatter(self.__column, self.hasNull)

        self.size = self.formatter.size
        self.min = self.formatter.min
        self.max = self.formatter.max
        if self.__force_format:
            self.size = self.__force_format.size
            self.formatter.size = self.__force_format.size
            self.formatter.format = self.__force_format.format
            self.formatter.fortran_format = self.__force_format.fortran_format
            self.formatter.out_format = self.__force_format.out_format
            self.formatter.none_format = self.__force_format.none_format

    def __get_type(self):
        if self.__column.dtype.name.startswith("i"):
            return int
        elif self.__column.dtype.name.startswith("f"):
            return float
        elif self.__column.dtype.name.startswith("u"):
            return int
        elif self.__column.dtype.name.startswith("s"):
            return str

        # if contains null values dtype is 'object'
        #if self.hasNull is False:
        #    return str

        for value in self.__column:
            if value is None:
                continue
            if isinstance(value, numpy.ma.core.MaskedConstant):
                continue
            if isinstance(value, (int, numpy.integer)):
                return int
            elif numpy.isreal(value):
                return float
            else:
                return str
        return str

    def getName(self) -> str:
        """ get the name in input
        :return: the name (which can be different that original Column.name)
        """
        return self.__dbname

    def isSexa(self) -> bool:
        """Return True if the column is in sexagesimal
        :return: True is Sexagesimal
        """
        if self.__sexagesimal[0] is None and self.__sexagesimal[1] is None:
            return False
        return True

    def isSexaRA(self):
        """Return True if Right ascension
        :return: True if Sexagesimal RA
        """
        if self.__sexagesimal[0] is not None:
            return True
        return False

    def isSexaDE(self):
        """Return True if Declination
        :return: True if sexagesimal declination
        """
        if self.__sexagesimal[1] is not None:
            return True
        return False

    def setSexaRa(self, precision: int = 0):
        """set column as a Sexagesimal Right Ascension
        :param precision: number of seconds decimals
        """
        self.parse()

        if self.formatter.fortran_format[0] != 'A': raise Exception("bad sexagesimal format")
        if precision == 0:
            if self.size > 9:
                precision = self.size - 9  # HH:MM:SS.ss
        else:
            if self.size < 9 + precision or self.size > 10 + precision:
                raise Exception("bad sexagesimal format or bad precision (format: [+-]dd mm ss[.ss])")

        self.formatter.fortran_format = self.formatter.fortran_format.replace('A', 'R')
        self.__sexagesimal[0] = SexaRA(precision)

    def setSexaDe(self, precision: int = 0):
        """set column as a Sexagesimal Declination
        :param precision: number of seconds decimals
        """
        self.parse()
        if self.formatter.fortran_format[0] != 'A': raise Exception("bad sexagesimal format")
        if precision == 0:
            if self.size > 9:
                precision = self.size - 10  # [+-]HH:MM:SS.ss
        else:
            if self.size < 9 + precision or self.size > 10 + precision:
                raise Exception("bad sexagesimal format or bad precision (format: [+-]dd mm ss[.ss])")

        self.formatter.fortran_format = self.formatter.fortran_format.replace('A', 'D')
        self.__sexagesimal[1] = SexaDE(precision)

    def getSexaRA(self):
        """get Sexagesimal  Right ascension
        :return: SexaMeta information
        """
        return self.__sexagesimal[0]

    def getSexaDE(self):
        """get Sexagesimal DEclination
        :return: SexaMeta information
        """
        return self.__sexagesimal[1]

    def __has_null(self) -> bool:
        """Test null values in a numpy array
        """
        if isinstance(self.__column, MaskedColumn):
            return True
        #try:
        #    n = len(self.__column[self.__column.argmin(fill_value=0)])
        #except: #fail inf MaskedColumn
        #    return True
        return False

    def __get_unit(self) -> str:
        s = self.name.lower()
        if s.find("magnitude") > -1:
            return "mag"
        elif s.find("[ (]days[ )]") > -1:
            return "d"
        return UNDEFINED_UNIT

    def value(self, nrec: int) -> str:
        """get formatted value
        :return: the value formatted (ready fo ASCII aligned format)
        """
        return self.formatter.write(self.__column[nrec])

    def __str__(self) -> str:
        return "<CDSColumn name={0}>".format(self.name)
