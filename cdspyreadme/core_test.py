import unittest
import os, sys, time


from astropy.table import Table, MaskedColumn
import numpy as np
from cdspyreadme.core import CDSTablesMaker, CDSMRTTable, CDSAsciiTable

RESULT_DIR = "tests"
RESULT_PREFIX = "out_table"
CSV_TABLE = "table.csv"
ASCII_TABLE = "table.ascii"
MRT_TABLE = "table.mrt"


def __init_env():
    sys.stderr.write("change current dir to {}\n".format(RESULT_DIR))
    os.chdir(RESULT_DIR)

    for f in (CSV_TABLE, ASCII_TABLE, MRT_TABLE):
        if os.path.exists(f) is False:
            sys.stderr.write("error file {} not found for tests\n".format(f))
            sys.exit(1)

    # remove old files
    files = [f for f in os.listdir('.') if f.find(RESULT_PREFIX) == 0]
    for f in files:
        os.remove(f)

__init_env()


class ReadMeCase(unittest.TestCase):

    def __init__(self, argv):
        self.time = time.time()
        unittest.TestCase.__init__(self,argv)

    def __get_time(self, flag=""):
        newtime = time.time()
        t = newtime - self.time
        self.time = newtime
        sys.stderr.write("%elapse time {0} = {1}\n".format(flag, t))
        return t

    def __test_file(self, filename, min_line=0):
        if os.path.exists(filename) is False:
            return False
        if min_line > 0:
            with open(filename, "r") as fd:
                nline = 1
                for line in fd:
                    if nline >= min_line:
                        return True
                    nline +=1
            return False
        return True

    def __result_file(suffix):
        return "{}_{}".format(RESULT_PREFIX, suffix)

    def __readme_file(suffix):
        return "ReadMe_{}".format(suffix)

    def test_numpy(self):
        """test numpy table"""
        tablemaker = CDSTablesMaker()
        ntab = np.array([(1.1, 1, 'test'), (2.2, 2, 'test2'), (3.3, 3, None), (None, -1, '')],
                        dtype=[('mag', np.float), ('recno', np.int32), ('comment', np.str_, 10)])
        tablename = ReadMeCase.__result_file("numpy")
        tablemaker.addTable(ntab,
                            name=tablename,
                            description="test numpy table")
        tablemaker.writeCDSTables()
        self.assertTrue(self.__test_file(tablename, 3), "numpy file")

        readme = ReadMeCase.__readme_file("numpy")
        with open(readme, "w") as fd:
            tablemaker.makeReadMe(out=fd)

        self.assertTrue(self.__test_file(readme, 5), "Readme for numpy")
        sys.stderr.write("generate {0}/{1} {0}/{2}\n".format(RESULT_DIR, tablename, readme))

    def test_astropy(self):
        """test astropy table
           + test null_value
           + test masked column
        """
        tablemaker = CDSTablesMaker()
        ntab = Table([(1, 2, None, 4, 999),
                      (4.0, 1115.0, np.NaN, None, 999),
                      (1.1, 2., 999, 12.3, 12.3),
                      (-1.001, 2., 0, -99.12, np.NaN)],
                     names=['a', 'b', 'col3', 'col4'])
        ntab["col4"] = MaskedColumn(ntab["col4"], mask=[(val > 0) for val in ntab["col4"]])

        tablename = ReadMeCase.__result_file("astropy")
        table = tablemaker.addTable(ntab,
                            name=tablename,
                            description="test astropy table",
                            nullvalue=999)
        col = table.get_column("col4")
        col.set_null_value(-1.001)

        tablemaker.writeCDSTables()
        self.assertTrue(self.__test_file(tablename, 3), "astropy file")

        readme = ReadMeCase.__readme_file("astropy")
        with open(readme, "w") as fd:
            tablemaker.makeReadMe(out=fd)
        self.assertTrue(self.__test_file(readme, 5), "Readme for astropy")
        sys.stderr.write("generate {0}/{1} {0}/{2}\n".format(RESULT_DIR, tablename, readme))

    def test_file(self):
        """test csv file + ascii file
        + test sexa
        """
        tablemaker = CDSTablesMaker()

        self.__get_time("test_file")
        table1 = ReadMeCase.__result_file("file1")
        table = tablemaker.addTable(CSV_TABLE, table1, description="test")
        self.__get_time("test_file add_table")
        ra = table.get_column("col4")
        ra.setSexaRa()
        dec = table.get_column("col5")
        dec.setSexaDe()

        table2 = ReadMeCase.__result_file("file2")
        ascii_table = CDSAsciiTable(ASCII_TABLE, table2, description="test")
        tablemaker.addTable(ascii_table, table2, description="test")

        tablemaker.writeCDSTables()
        self.__get_time("test_file write")
        self.assertTrue(self.__test_file(table1, 3), "csv/ascii file")
        self.assertTrue(self.__test_file(table2, 3), "ascii file")

        readme = ReadMeCase.__readme_file("file")
        with open(readme, "w") as fd:
            tablemaker.makeReadMe(out=fd)

        self.assertTrue(self.__test_file(readme, 5), "csv/ascii for astropy")
        sys.stderr.write("generate {0}/{1} {0}/{2} {0}/{3}\n".format(RESULT_DIR, table1, table2, readme))
        #sys.stderr.write("generate {0}/{1} {0}/{2}\n".format(RESULT_DIR, table1, readme))

    def test_mrt(self):
        if self.__test_file("table.mrt", 3) is False:
            self.assertTrue(False, "MRT test not executes (table.mrt is not found)")
            return

        tablemaker = CDSTablesMaker()
        tablename = ReadMeCase.__result_file("mrt")
        mrt = CDSMRTTable(MRT_TABLE,
                          tablename,
                          description="test title MRT",
                          set_limit=True)
        tablemaker.addTable(mrt,
                            name=tablename,
                            description="test astropy table")
        tablemaker.writeCDSTables()
        self.assertTrue(self.__test_file(tablename, 3), "MRT file")

        readme = ReadMeCase.__readme_file("mrt")
        with open(readme, "w") as fd:
            tablemaker.makeReadMe(out=fd)

        self.assertTrue(self.__test_file(readme, 5), "Readme for MRT")
        sys.stderr.write("generate {0}/{1} {0}/{2}\n".format(RESULT_DIR, tablename, readme))

    def test_generate_mrt(self):
        tablemaker = CDSTablesMaker()
        tablemaker.title = "catalogue title!..."
        tablemaker.author = 'G.Landais'
        tab = Table([(1, 2, None, 4, 999),
                      (4.0, 1115.0, np.NaN, None, 999),
                      (1.1, 2., 999, 12.3, 12.3),
                      (-1.001, 2., 0, -99.12, np.NaN)],
                     names=['a', 'b', 'col3', 'col4'])
        tablename = ReadMeCase.__result_file("mrt_table")
        col = tab['a']
        col.name = 'ee'
        col.description = "test"
        col.unit = "deg"
        print(tab)
        tablemaker.addTable(tab,
                            name=tablename,
                            description="test MRT generation")
        tablemaker.toMRT()
        self.assertTrue(self.__test_file(tablename, 3), "MRT generation")