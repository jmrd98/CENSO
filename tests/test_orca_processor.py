import unittest
import unittest.mock
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
from censo.procfact import ProcessorFactory
os.chdir(os.path.split(__file__)[0])
test_dir = os.getcwd()


class OrcaProcTest(unittest.TestCase):
    def setUp(self):
        self.processor = ProcessorFactory.create_processor("orca", test_dir)

    @unittest.mock.patch("censo.orca_processor.OrcaProc._make_call")
    def test_sp(self, mock_call):
        mock_call.return_value = 0
