#from .GromacsParser import GromacsParserInterface as GROMACSParser
from .metainfo import m_env

from nomad.parsing.parser import FairdiParser
from gromacsparser.gromacs_parser import GromacsOutput


class GromacsParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/gromacs', code_name='Gromacs', code_homepage='http://www.gromacs.org/',
            domain='dft', mainfile_contents_re=r'gmx mdrun, VERSION')

    def parse(self, filepath, archive, logger=None):
        self._metainfo_env = m_env

        parser = GromacsOutput(filepath, archive, logger)

        parser.parse()
