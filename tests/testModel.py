#!/usr/bin/env python
"""A simple test of whether the fitter is working well.

This should be turned into a unit test that reports pass/fail.
"""
import numpy
import matplotlib.pyplot as pyplot
from . import fitData

posArr = numpy.arange(-1.0, 1.01, 0.1, dtype=float)
nomPos = numpy.zeros((len(posArr)**2,2))
i = 0
for x in posArr:
    for y in posArr:
        nomPos[i] = (x, y)
        i += 1

pyplot.ion()

print "QuadrupoleModel"
qpModel = fitData.QuadrupoleModel()
qpModel.setMagnitudeAngle(0.1, 45.0)
print "desired coeffs=", qpModel.getCoeffs()
qpPos = qpModel.apply(nomPos)

for testCoeffs in ((1, 1), (10, 10), (0.5, -45)):
    print "init coeffs=", testCoeffs
    qpFit = fitData.ModelFit(
        model = fitData.QuadrupoleModel(testCoeffs),
        measPos = qpPos,
        nomPos = nomPos,
        doRaise = True,
    )
    print "fit coeffs=", qpFit.model.getCoeffs()


print "TransRotScaleModel"
trsModel = fitData.TransRotScaleModel((0.2, -0.3, 0.8, -0.9))
trsPos = trsModel.apply(nomPos)
print "desired coeffs=", trsModel.getCoeffs()

for testCoeffs in ((1, 1, 1, 1), (10, 10, 10, 10), (0.5, -0.6, 0.5, -0.5)):
    print "init coeffs=", testCoeffs
    trsFit = fitData.ModelFit(
        model = fitData.TransRotScaleModel(testCoeffs),
        measPos = trsPos,
        nomPos = nomPos,
        doRaise = True,
    )
    print "fit coeffs=", trsFit.model.getCoeffs()


# pyplot.plot(
#     nomPos[:,0], nomPos[:,1], "b.",
#     qpPos[:,0], qpPos[:,1], "r.",
# )
