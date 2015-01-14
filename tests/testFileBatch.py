#!/usr/bin/env python
"""Test a batch of files (every file in tests/testData), assure that they are parsed successfully
and have the correct type of information stored
"""
import glob
import os
import unittest
import datetime
from fitPlugPlateMeas.cmmPlateMeas import PlateMeas

pathToThisFile = os.path.realpath(__file__)

# each parsed file object must have these attributes, types

class TestFileBatch(unittest.TestCase):

    def testFileBatch(self):
        # find test file directory
        testFileDirectory = os.path.dirname(pathToThisFile)+"/testData"
        # get every filein the test file directory
        allFiles = glob.glob(testFileDirectory + "/D*") # get all files beginning with D
        for f in allFiles:
            # verify each file successfully parses!
            print "testing ", f
            parsedFileObj = PlateMeas(f)
            # now test all attributes
            # by running all methods that start with _verify
            # defined below
            for verifyAttr in dir(self):
                if verifyAttr.startswith("_verify"):
                    getattr(self, verifyAttr)(parsedFileObj)

    def _verifyPlateID(self, parsedFileObj):
        # check that the parsedFileObj has an attribute of plateID
        self.assertTrue(hasattr(parsedFileObj, "plateID"))
        # and that that value is an integer
        self.assertTrue(type(parsedFileObj.plateID)==int)

    def _verifyDate(self, parsedFileObj):
        # check that the parsedFileObj has an attribute of date
        self.assertTrue(hasattr(parsedFileObj, "date"))
        # and that that value is a datetime.date obj
        self.assertTrue(type(parsedFileObj.date)==datetime.date)

    def _verifyFitOffset(self, parsedFileObj):
        # check that the parsedFileObj has an attribute of fitOffset
        self.assertTrue(hasattr(parsedFileObj, "fitOffset"))
        # and that that value is a list of 2 floats
        self.assertTrue(len(parsedFileObj.fitOffset)==2)
        for offset in self.parserFileObj.fitOffset:
            self.assertTrue(type(offset)==float)

    def _verifyFitScale(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "fitScale"))
        self.assertTrue(type(parsedFileObj.fitScale) == float)

    def _verifyFitRotAngle(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "fitRotAngle"))
        self.assertTrue(type(parsedFileObj.fitRotAngle) == float)

    def _verifyQPMagnitude(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "qpMagnitude"))
        self.assertTrue(type(parsedFileObj.qpMagnitude) == float)

    def _verifyQPAngle(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "qpAngle"))
        self.assertTrue(type(parsedFileObj.qpAngle) == float)

    def _verifyResidRadErrRMS(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "residRadErrRMS"))
        self.assertTrue(type(parsedFileObj.residRadErrRMS) == float)

    def _verifyMangaResidRadErrRMS(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "mangaResidRadErrRMS"))
        # either None or a float
        try:
            self.assertTrue(parsedFileObj.mangaResidRadErrRMS is not None)
        except:
            self.assertTrue(type(parsedFileObj.mangaResidRadErrRMS) == float)

    def _verifyDiaErrRMS(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "diaErrRMS"))
        self.assertTrue(type(parsedFileObj.diaErrRMS) == float)

    def _verifyMangaDiaRadErrRMS(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "mangaDiaRadErrRMS"))
        # either None or a float
        try:
            self.assertTrue(parsedFileObj.mangaDiaRadErrRMS is not None)
        except:
            self.assertTrue(type(parsedFileObj.mangaDiaRadErrRMS) == float)

    def _verifyMaxDiaError(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "maxDiaError"))
        self.assertTrue(type(parsedFileObj.maxDiaError) == float)

    def _verifyQPResidRadErrorRMS(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj, "qpResidRadErrorRMS"))
        self.assertTrue(type(parsedFileObj.qpResidRadErrorRMS) == float)

    # these are dataTable tests

    def _verifyNomX(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "nomX"))
        for val in self.dataTable.nomX:
            self.assertTrue(type(val) == float)

    def _verifyNomY(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "nomY"))
        for val in self.dataTable.nomY:
            self.assertTrue(type(val) == float)

    def _verifyMeasX(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "measX"))
        for val in self.dataTable.measX:
            self.assertTrue(type(val) == float)

    def _verifyMeasY(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "measY"))
        for val in self.dataTable.measY:
            self.assertTrue(type(val) == float)

    def _verifyResidX(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "residX"))
        for val in self.dataTable.residX:
            self.assertTrue(type(val) == float)

    def _verifyResidY(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "residY"))
        for val in self.dataTable.residY:
            self.assertTrue(type(val) == float)

    def _verifyResidRad(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "residRad"))
        for val in self.dataTable.residRad:
            self.assertTrue(type(val) == float)

    def _verifyNomDia(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "nomDia"))
        for val in self.dataTable.nomDia:
            self.assertTrue(type(val) == float)

    def _verifyDiaError(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "diaErr"))
        for val in self.dataTable.diaErr:
            self.assertTrue(type(val) == float)

    def _verifyQPResidX(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "qpResidX"))
        for val in self.dataTable.qpResidX:
            self.assertTrue(type(val) == float)

    def _verifyQPResidY(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "qpResidY"))
        for val in self.dataTable.qpResidY:
            self.assertTrue(type(val) == float)

    def _verifyQPResidRad(self, parsedFileObj):
        self.assertTrue(hasattr(parsedFileObj.dataTable, "qpResidRad"))
        for val in self.dataTable.qpResidRad:
            self.assertTrue(type(val) == float)


if __name__ == '__main__':
    # this is how to automatically run the unittests.
    from unittest import main
    main()