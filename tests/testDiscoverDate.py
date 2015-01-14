"""Test for the plate 8131 measurement file
"""
import os
import unittest
from fitPlugPlateMeas.fitPlugPlateMeasWdg import readHeader

pathToThisFile = os.path.realpath(__file__)

# determine the path to the 8131 plate measurement file in the testData subdirectory
pathToPlate1413MeasFile = os.path.dirname(pathToThisFile)+"/testData/D1413_1"

class TestDiscoverDate(unittest.TestCase):
    """Every method beginning with "test" is run automatically by the
    unittest package

    all self.assert* methods are provided by the unittest package
    """
    def testDiscoverDate(self):
        """This runs before each test method.
        """
        # this file has no date written in the header
        # but it should still succeed.
        numHeaderLines, plateID, measDate = readHeader(pathToPlate1413MeasFile)
        self.assertTrue(measDate is None)

if __name__ == '__main__':
    # this is how to automatically run the unittests.
    from unittest import main
    main()


