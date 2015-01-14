"""Test for the plate 8131 measurement file
"""
import os
import unittest
import datetime
import itertools
from fitPlugPlateMeas.cmmPlateMeas import PlateMeas

pathToThisFile = os.path.realpath(__file__)

# determine the path to the 8131 plate measurement file in the testData subdirectory
pathToPlate8131MeasFile = os.path.dirname(pathToThisFile)+"/testData/D8131_1"

class TestPlate8131(unittest.TestCase):
    """Every method beginning with "test" is run automatically by the
    unittest package

    all self.assert* methods are provided by the unittest package
    """
    def setUp(self):
        """This runs before each test method.
        """
        # parse the file
        self.plateMeas = PlateMeas(pathToPlate8131MeasFile)

    def testPlateID(self):
        self.assertTrue(self.plateMeas.plateID == 8131)

    def testDate(self):
        # datetime docs: https://docs.python.org/2/library/datetime.html
        self.assertTrue(self.plateMeas.date == datetime.date(2014, 10, 14))

    def testFitOffset(self):
        fitOffset = [-.044, 0.079]
        # use almost equal for floats
        # (machine precision may bit you when comparing floats)
        self.assertAlmostEqual(self.plateMeas.fitOffset[0], fitOffset[0])
        self.assertAlmostEqual(self.plateMeas.fitOffset[1], fitOffset[1])
        self.assertTrue(len(self.plateMeas.fitOffset)==2)

    def testFitScale(self):
        fitScale = 1.000007
        self.assertAlmostEqual(fitScale, self.plateMeas.fitScale)

    def testFitRotAngle(self):
        fitRotAngle = 0.029
        self.assertAlmostEqual(fitRotAngle, self.plateMeas.fitRotAngle)

    def testQPMagnitude(self):
        qpMagnitude = 0.00001041
        self.assertAlmostEqual(qpMagnitude, self.plateMeas.qpMagnitude)

    def testQPAngle(self):
        qpAngle = -69.59
        self.assertAlmostEqual(qpAngle, self.plateMeas.qpAngle)

    def testResidRadErrRMS(self):
        residRadErrRMS = 0.0083
        self.assertAlmostEqual(residRadErrRMS, self.plateMeas.residRadErrRMS)

    def testMangaResidRadErrRMS(self):
        mangaResidRadErrRMS = 0.0079
        self.assertAlmostEqual(mangaResidRadErrRMS, self.plateMeas.mangaResidRadErrRMS)

    def testDiaErrRMS(self):
        diaErrRMS = 0.0078
        self.assertAlmostEqual(diaErrRMS, self.plateMeas.diaErrRMS)

    def testMangaDiaErrRMS(self):
        mangaDiaErrRMS = 0.0103
        self.assertAlmostEqual(mangaDiaErrRMS, self.plateMeas.mangaDiaErrRMS)

    def testMaxDiaErr(self):
        maxDiaErr = 0.015
        self.assertAlmostEqual(maxDiaErr, self.plateMeas.maxDiaErrRMS)

    def testQPResidRadErrRMS(self):
        qpResidRadErrRMS = 0.0079
        self.assertAlmostEqual(qpResidRadErrRMS, self.plateMeas.qpResidRadErrRMS)

    ###### following tests data table ###########

    def testNomX(self):
        nomX = [
            11.650, -0.079, 5.070, 210.257, 264.541, 316.010, -132.308, -102.618, -178.579,
            -235.190,-182.381,-110.285,-154.584,-189.156,-155.137,-198.382,-271.072,-235.999,
            -193.162,-147.968, -79.842, -64.297, 201.910, 223.395, 228.753, 135.628, 109.986,
            134.547,179.880,175.571,83.809,-30.895,-27.349,-86.281,-91.691, 119.757, 198.435,
            242.026,293.150,308.559,283.031,131.175, 90.484, 43.171, 51.761, 111.554, 93.949,
            63.219, 0.597, 12.398, 107.478, 84.021, -39.110, -56.805, -37.636, -17.194, -30.470,
            -255.748, -266.833, -63.506, -5.191, 16.687, -32.305, 19.434, 48.823, -67.416, -3.104,
            27.537, 11.262, -102.338, -108.337, -137.988, -135.251, 127.312, 169.627, 215.792, 133.809,
            116.670, 75.657, -181.334, -283.440, -293.342,
        ]

        for nomX_i, dataTableNomX_i in itertools.izip(nomX, self.plateMeas.nomX):
            # compare each element in list
            # izip iterates over multiple lists at once, eg:
            # izip([x1,x2,x3],[y1,y2,y3]):
            #   iter 1: x1, y1
            #   iter 2: x2, y2
            #   iter 3: x3, y3
            self.assertAlmostEqual(nomX_i, dataTableNomX_i)

    def testNomY(self):
        nomY = [
            -10.184, 162.760, 185.377, 72.308, 61.924, 71.483, -296.354, -284.962, -236.239,
            -186.346, -107.976, 198.064, 175.594, 173.351, 116.867, 105.085, 70.854, 8.442,
            11.148, 51.402, 117.334, 91.382, 214.878, 186.367, 137.729, 86.215, -131.470,
            -93.187, -122.949, -160.561, -231.492, -109.104, -63.378, -17.920, -10.636, -248.607,
            -227.725, -157.622, -105.030, -37.176, 127.990, 283.000, 263.175, 291.045, 254.908,
            -42.420, -53.546, -88.703, -113.476, -62.848, 163.114, 203.659, 263.602, 217.805,
            236.433, 217.420, 128.991, -75.264, -120.317, -315.439, -314.039, -320.223, -222.554,
            -245.030, -219.837, 40.756, 77.874, -11.879, -33.993, -113.114, -163.872, -243.535,
            -259.715, 12.078, 16.560, -6.787, -80.166, -52.234, 6.799, 244.927, 144.037, 88.194
        ]


        for nomY_i, dataTableNomY_i in itertools.izip(nomY, self.plateMeas.nomY):
            self.assertAlmostEqual(nomY_i, dataTableNomY_i)

    def testMeasX(self):
        measX = [
            11.610, -0.207, 4.921, 210.177, 264.461, 315.925, -132.201, -102.510, -178.504,
            -235.142, -182.378, -110.429, -154.720, -189.286, -155.243, -198.476, -271.154,
            -236.045, -193.213, -148.036, -79.939, -64.388, 201.753, 223.256, 228.638, 135.537,
            110.017, 134.556, 179.896, 175.610, 83.887, -30.884, -27.361, -86.316, -91.725,
            119.848, 198.521, 242.063, 293.167, 308.537, 282.918, 130.985, 90.301, 42.971, 51.582,
            111.533, 93.931, 63.221, 0.608, 12.386, 107.351, 83.877, -39.294, -56.959, -37.803,
            -17.354, -30.581, -255.751, -266.821, -63.384, -5.065, 16.813, -32.231, 19.523, 48.895,
            -67.473, -3.187, 27.497, 11.231, -102.328, -108.303, -137.906, -135.163, 127.263, 169.577,
            215.747, 133.804, 116.652, 75.610, -181.514, -283.556, -293.427
        ]


        for measX_i, dataTablemeasX_i in itertools.izip(measX, self.plateMeas.measX):
            self.assertAlmostEqual(measX_i, dataTablemeasX_i)

    def testMeasY(self):
        measY = [
            -10.090, 162.841, 185.449, 72.484, 62.131, 71.719, -296.356, -284.940, -236.258,
            -186.389, -107.983, 198.090, 175.598, 173.337, 116.873, 105.067, 70.810, 8.413,
            11.141, 51.412, 117.374, 91.432, 215.056, 186.552, 137.918, 86.353, -131.337,
            -93.040, -122.787, -160.398, -231.373, -109.035, -63.306, -17.876, -10.592, -248.470,
            -227.548, -157.425, -104.798, -36.937, 128.209, 283.146, 263.295, 291.148, 255.016,
            -42.284, -53.417, -88.581, -113.398, -62.750, 163.248, 203.788, 263.660 , 217.855,
            236.492, 217.487, 129.051, -75.304, -120.363, -315.403, -313.977, -320.141, -222.497,
            -244.947, -219.736, 40.810, 77.955, -11.780, -33.902, -113.086, -163.847, -243.534,
            -259.714, 12.218, 16.714, -6.607, -80.024, -52.096, 6.917, 244.913, 143.982, 88.135
        ]


        for measY_i, dataTablemeasY_i in itertools.izip(measY, self.plateMeas.measY):
            self.assertAlmostEqual(measY_i, dataTablemeasY_i)

    def testResidX(self):
        residX = [
            -0.002, -0.002, -0.012, -0.001, -0.007, -0.008,  0.001,  0.008,
            0.000, -0.001, -0.007,  0.001, -0.002,  0.003, -0.002,
            0.004, -0.001,  0.004, -0.000,  0.003,  0.007, -0.001, -0.006,
            -0.003, -0.003, -0.005,  0.007,  0.004, -0.004, -0.000,  0.004,
            -0.000, -0.000,  0.000,  0.005,  0.008,  0.013, -0.001,  0.005,
            0.000, -0.007, -0.004, -0.007, -0.009, -0.007,  0.000, -0.002,  0.000,
            -0.003, -0.000, -0.002,  0.002, -0.007,  0.000, -0.003, -0.006, -0.002,
            0.004, -0.003,  0.006,  0.011,  0.007,  0.005,  0.008,  0.004,  0.008,
            0.000, -0.003, -0.005, -0.003, -0.005,  0.003,  0.001, -0.000,  0.001,
            -0.006, -0.003, -0.002, -0.001, -0.011,  0.003,  0.005
        ]


        for residX_i, dataTableresidX_i in itertools.izip(residX, self.plateMeas.residX):
            self.assertAlmostEqual(residX_i, dataTableresidX_i)

    def testResidY(self):
        residY = [
            0.009,  0.001, -0.011, -0.010, -0.007, -0.004, -0.012, -0.003, -0.006,
            -0.002,  0.007,  0.001,  0.002,  0.001,  0.005,  0.003,  0.014,  0.011,
            0.012,  0.005,  0.000,  0.003, -0.005, -0.009, -0.007, -0.011, -0.001,
            0.000, -0.007, -0.004, -0.001,  0.006,  0.007,  0.009,  0.011, -0.001,
            -0.001, -0.004,  0.005,  0.004, -0.004, -0.002, -0.007, -0.000,  0.001,
            0.001,  0.003,  0.011, -0.001,  0.013, -0.001,  0.006, -0.003, -0.002,
            -0.003, -0.005, -0.005,  0.011,  0.011, -0.009, -0.012, -0.003, -0.004,
            -0.004, -0.001,  0.009,  0.003,  0.006,  0.006,  0.001,  0.002, -0.007,
            -0.008, -0.004, -0.011, -0.008, -0.004,  0.000,  0.000, -0.003,  0.008,  0.010
        ]


        for residY_i, dataTableresidY_i in itertools.izip(residY, self.plateMeas.residY):
            self.assertAlmostEqual(residY_i, dataTableresidY_i)

    def testResidRad(self):
        residRad = [
            0.009, 0.002, 0.016, 0.010, 0.010, 0.008, 0.012, 0.009, 0.006, 0.002,
            0.010, 0.001, 0.003, 0.003, 0.005, 0.005, 0.014, 0.012, 0.012, 0.006,
            0.007, 0.003, 0.008, 0.009, 0.008, 0.012, 0.007, 0.004, 0.008, 0.004,
            0.004, 0.006, 0.007, 0.009, 0.012, 0.008, 0.013, 0.004, 0.007, 0.004,
            0.008, 0.004, 0.010, 0.009, 0.007, 0.001, 0.003, 0.011, 0.003, 0.013,
            0.002, 0.006, 0.007, 0.002, 0.004, 0.008, 0.005, 0.012, 0.011, 0.011,
            0.016, 0.008, 0.007, 0.009, 0.004, 0.012, 0.003, 0.007, 0.008, 0.003,
            0.005, 0.007, 0.008, 0.004, 0.011, 0.011, 0.005, 0.002, 0.001, 0.012,
            0.009, 0.011
        ]


        for residRad_i, dataTableresidRad_i in itertools.izip(residRad, self.plateMeas.residRad):
            self.assertAlmostEqual(residRad_i, dataTableresidRad_i)

    def testNomDia(self):
        nomDia = [
            3.280, 2.167, 2.167, 2.167, 2.167, 2.167, 2.167, 2.167, 3.280, 3.280,
            2.167, 2.167, 2.167, 2.167, 2.167, 2.167, 3.280, 2.167, 2.167, 2.167,
            2.167, 3.280, 2.167, 2.167, 2.167, 2.167, 3.280, 3.280, 2.167, 2.167,
            2.167, 2.167, 3.280, 2.167, 3.280, 2.167, 3.280, 2.167, 3.280, 2.167,
            2.167, 2.167, 2.167, 3.280, 3.280, 3.280, 3.280, 2.167, 2.167, 2.167,
            3.280, 3.280, 2.167, 2.167, 3.280, 2.167, 2.167, 2.167, 2.167, 2.167,
            3.280, 2.167, 2.167, 3.280, 2.167, 2.167, 3.280, 3.280, 2.167, 2.167,
            2.167, 3.280, 3.280, 3.280, 3.280, 2.167, 2.167, 3.280, 3.280, 3.280,
            2.167, 2.167
        ]


        for nomDia_i, dataTablenomDia_i in itertools.izip(nomDia, self.plateMeas.nomDia):
            self.assertAlmostEqual(nomDia_i, dataTablenomDia_i)

    def testDiaErr(self):
        diaErr = [
            -0.015,  0.007,  0.011,  0.007,  0.009,  0.007,  0.004,  0.009, -0.010,
            -0.007,  0.008,  0.010,  0.009,  0.015,  0.011,  0.010, -0.009,  0.009,
            0.009,  0.009,  0.009, -0.012,  0.008,  0.008,  0.006,  0.005, -0.010,
            -0.012,  0.006,  0.006,  0.007,  0.006, -0.010,  0.008, -0.012,  0.007,
            -0.011,  0.006, -0.009,  0.005,  0.008,  0.010,  0.010, -0.007, -0.009,
            -0.013, -0.011,  0.005,  0.005,  0.005, -0.010, -0.009,  0.008,  0.006,
            -0.010,  0.008,  0.007,  0.010,  0.008,  0.007, -0.010,  0.006,  0.006,
            -0.007,  0.006,  0.005, -0.012, -0.013,  0.003,  0.004,  0.005, -0.009,
            -0.009, -0.010, -0.010,  0.005,  0.004, -0.011, -0.011, -0.006,  0.009,  0.012
        ]


        for diaErr_i, dataTablediaErr_i in itertools.izip(diaErr, self.plateMeas.diaErr):
            self.assertAlmostEqual(diaErr_i, dataTablediaErr_i)

    def testQPResidX(self):
        qpResidX = [
            -0.002, -0.002, -0.012,  0.001, -0.004, -0.004,  0.001,  0.008, -0.001,
            -0.003, -0.009,  0.002, -0.001,  0.004, -0.001,  0.004, -0.002,  0.002,
            -0.002,  0.002,  0.007, -0.000, -0.005, -0.001, -0.001, -0.003,  0.006,
            0.004, -0.005, -0.001,  0.003, -0.000, -0.000, -0.001,  0.004,  0.007,
            0.011, -0.002,  0.006,  0.002, -0.004, -0.004, -0.007, -0.010, -0.007,
            0.000, -0.002, -0.000, -0.003, -0.000, -0.001,  0.002, -0.006,  0.001,
            -0.003, -0.006, -0.002,  0.002, -0.006,  0.007,  0.011,  0.007,  0.005,
            0.008,  0.003,  0.008,  0.000, -0.003, -0.005, -0.004, -0.005,  0.003,
            0.001,  0.001,  0.002, -0.005, -0.003, -0.002,  0.000, -0.009,  0.003,  0.005
        ]


        for qpResidX_i, dataTableqpResidX_i in itertools.izip(qpResidX, self.plateMeas.qpResidX):
            self.assertAlmostEqual(qpResidX_i, dataTableqpResidX_i)

    def testQPResidY(self):
        qpResidY = [
            0.009, -0.001, -0.013, -0.009, -0.006, -0.003, -0.012, -0.003,
            -0.007, -0.003,  0.006, -0.001,  0.000,  0.000,  0.004,  0.002,
            0.014,  0.011,  0.012,  0.005, -0.001,  0.002, -0.004, -0.007,
            -0.006, -0.010,  0.000,  0.001, -0.007, -0.003,  0.001,  0.007,
            0.007,  0.008,  0.011,  0.002,  0.001, -0.003,  0.005,  0.004,
            -0.003, -0.002, -0.007, -0.002, -0.000,  0.001,  0.003,  0.012,
            0.000,  0.014, -0.000,  0.006, -0.006, -0.004, -0.005, -0.007,
            -0.006,  0.010,  0.010, -0.007, -0.010, -0.001, -0.003, -0.002,
            0.001,  0.009,  0.002,  0.006,  0.007,  0.001,  0.001, -0.007,
            -0.008, -0.004, -0.011, -0.008, -0.004,  0.000,  0.000, -0.005,
            0.008,  0.010
        ]


        for qpResidY_i, dataTableqpResidY_i in itertools.izip(qpResidY, self.plateMeas.qpResidY):
            self.assertAlmostEqual(qpResidY_i, dataTableqpResidY_i)

    def testQPResidRad(self):
        qpResidRad = [
            0.009, 0.002, 0.017, 0.009, 0.007, 0.005, 0.012, 0.009, 0.007, 0.005,
            0.010, 0.002, 0.001, 0.004, 0.004, 0.005, 0.014, 0.011, 0.012, 0.006,
            0.007, 0.002, 0.006, 0.007, 0.006, 0.010, 0.006, 0.004, 0.008, 0.003,
            0.003, 0.007, 0.007, 0.008, 0.012, 0.007, 0.011, 0.004, 0.008, 0.004,
            0.005, 0.004, 0.010, 0.010, 0.007, 0.001, 0.004, 0.012, 0.003, 0.014,
            0.001, 0.006, 0.009, 0.004, 0.006, 0.009, 0.006, 0.010, 0.011, 0.010,
            0.014, 0.007, 0.006, 0.008, 0.003, 0.012, 0.002, 0.007, 0.008, 0.004,
            0.005, 0.008, 0.008, 0.004, 0.011, 0.010, 0.005, 0.002, 0.001, 0.011,
            0.009, 0.011
        ]


        for qpResidRad_i, dataTableqpResidRad_i in itertools.izip(qpResidRad, self.plateMeas.qpResidRad):
            self.assertAlmostEqual(qpResidRad_i, dataTableqpResidRad_i)

if __name__ == '__main__':
    # this is how to automatically run the unittests.
    from unittest import main
    main()


