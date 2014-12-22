#!/usr/bin/env python
from __future__ import with_statement
"""Fit plug plate measurements

History:
2010-06-17 ROwen    2.0: first version using Python instead of Igor
2011-01-13 ROwen    2.1: Show another digit of radial position error.
                    Display radial position error as Pos Err instead of Rad Err.
                    Fit and show quadrupole moment as additional data.
2011-01-14 ROwen    2.2: Graph residuals after removing quadrupole error.
                    Normalize reported quadrupole magnitude (always positive) and angle (range -180 to 180).
                    Fixed name dipole -> quadrupole.
                    Slightly improved quadrupole fitting by using an initial value for quadrupole magnitude.
                    Improve display of long filenames (which are sometimes used for debugging).
2011-02-25 ROwen    Moved data fitting to the fitData module.
                    Modified to use RO.Wdg.DropletApp.
2011-03-18 ROwen    Save fit data to a file.
2011-10-11 ROwen    Renamed to avoid conflicting with package name. Moved __main__ elsewhere.
                    Modified to only process files whose name matches D[0-9][0-9]*
2012-01-30 ROwen    Bug fix: the output file contained had X data instead of Y data for NomY and MeasY.
2014-01-16 CCS      Reporting RMS stats for Manga (larger size) holes separately.
"""
import os.path
import re
import Tkinter
import numpy
import scipy.optimize
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import RO.Constants
import RO.StringUtil
import RO.Wdg
from . import fitData
from .version import __version__

__all__ = ["FitPlugPlateMeasWdg"]

class GraphWdg(Tkinter.Frame):
    def __init__(self, master):
        Tkinter.Frame.__init__(self, master)

        # dictionary of fileName: posErrArr
        self.dataDict = dict()
        self.enablePlotData = True

        RO.Wdg.StrLabel(master = self, text="Position Errors for:").grid(row=0, column=0)
        self.fileMenu = RO.Wdg.OptionMenu(
            master = self,
            items = [],
            callFunc = self.doPlateMenu,
        )
        self.fileMenu.grid(row=0, column=1)

        self.plotFig = matplotlib.figure.Figure(figsize=(19, 8), frameon=True)
        self.figCanvas = FigureCanvasTkAgg(self.plotFig, self)
        figCnvWdg = self.figCanvas.get_tk_widget()
        figCnvWdg.grid(row=1, column=0, columnspan=5, sticky="news")
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(4, weight=1)
        self.plotAxis = None
        self.quadrupolePlotAxis = None
        self.resetPlotAxis()

    def resetPlotAxis(self):
        if self.plotAxis:
            self.plotFig.delaxes(self.plotAxis)
            self.plotFig.delaxes(self.quadrupolePlotAxis)
        self.plotAxis = self.plotFig.add_subplot(
            1, 2, 1,
            title = "Residual errors",
            xlabel = "X (mm)",
            ylabel = "Y (mm)",
            autoscale_on = False,
            xlim = (-350, 350),
            ylim = (-350, 350),
            aspect = "equal",
        )
        self.quadrupolePlotAxis = self.plotFig.add_subplot(
            1, 2, 2,
            title = "Residuals after removing quadrupole",
            xlabel = "X (mm)",
            ylabel = "Y (mm)",
            autoscale_on = False,
            xlim = (-350, 350),
            ylim = (-350, 350),
            aspect = "equal",
        )

    def addData(self, fileName, posErrArr):
        """Add data to be graphed

        Inputs:
        - fileName: name of file containing the data
        - posErrArr: a structured array containing fields:
            "nomPos", "residPosErr" and "quadrupoleResidPosErr"
        """
        self.dataDict[fileName] = posErrArr
        plateIDs = sorted(self.dataDict.keys())
        self.enablePlotData = False
        try:
            self.fileMenu.setItems(plateIDs)
        finally:
            self.enablePlotData = True
        self.fileMenu.set(fileName)

    def doPlateMenu(self, dumWdg=None):
        fileName = self.fileMenu.getString()
        if fileName:
            self.plotData(self.dataDict[fileName])

    def plotData(self, posErrArr):
        if not self.enablePlotData:
            return
        self.resetPlotAxis()

        q = self.plotAxis.quiver(
            posErrArr["nomPos"][:,0],
            posErrArr["nomPos"][:,1],
            posErrArr["residPosErr"][:,0],
            posErrArr["residPosErr"][:,1],
            units="dots",
            width=1,
            angles="xy",
            scale=0.001)
        self.plotAxis.quiverkey(q, 0.9, 0.92, 0.01, "0.01 mm")

        q = self.quadrupolePlotAxis.quiver(
            posErrArr["nomPos"][:,0],
            posErrArr["nomPos"][:,1],
            posErrArr["quadrupoleResidPosErr"][:,0],
            posErrArr["quadrupoleResidPosErr"][:,1],
            units="dots",
            width=1,
            angles="xy",
            scale=0.001)
        self.quadrupolePlotAxis.quiverkey(q, 0.9, 0.92, 0.01, "0.01 mm")

        self.figCanvas.draw()


class FitPlugPlateMeasWdg(RO.Wdg.DropletApp):
    """Fit plug plate CMM measurements

    Fit out translation, rotation and scale and compute the residual position error to determine
    if the plate is good enough.
    Also fit quadrupole to that residual error to look for systematic error in the drilling machine
    (when this gets bad enough we recalibrate the machine or preprocess the plFanuc files to take it out).
    Log and graph the results.
    """
    def __init__(self, master, filePathList=None):
        """Construct a FitPlugPlateMeasWdg

        Inputs:
        - master: master widget; should be root
        - filePathList: list of files to process
        """
        RO.Wdg.DropletApp.__init__(self,
            master = master,
            width = 135,
            height = 20,
            font = "Courier 12", # want a fixed width font
            printTraceback = True,
            recursionDepth = 1,
            patterns = "D[0-9][0-9]*",
        )

        self.logWdg.addOutput("""Plug Plate Fitter %s

File       Meas Date  Holes  Offset X  Offset Y   Scale     Rotation  Pos Err   Dia Err   Pos Err (M)   Dia Err (M)   Dia Err  QPole Mag  QPole Ang  QP Pos Err
                                mm        mm                  deg     RMS mm    RMS mm     RMS mm        RMS mm        Max mm      1e-6       deg       RMS mm
""" % (__version__,))

        self.graphTL = RO.Wdg.Toplevel(
            master = self,
            title = "Position Errors",
            geometry = "+10+350",
            closeMode = RO.Wdg.tl_CloseDisabled,
        )
        self.graphWdg = GraphWdg(master = self.graphTL)
        self.graphWdg.grid(row=0, column=0, sticky="news")
        self.graphTL.grid_columnconfigure(0, weight=1)
        self.graphTL.grid_rowconfigure(0, weight=1)

        desDType = [
            ("measPos", float, (2,)),
            ("fitPos", float, (2,)),
        ]
        self.savedDataArr = numpy.zeros(0, dtype=desDType)
        self.numSavedPlates = 0

        if filePathList:
            self.processFileList(filePathList)

    def processFileList(self, filePathList):
        RO.Wdg.DropletApp.processFileList(self, filePathList)

        if self.numSavedPlates < 1:
            return
        try:
            # fit quadrupole to saved data and report
            fitQuadrupole = fitData.ModelFit(
                model = fitData.QuadrupoleModel([0.01, 0.0]),
                measPos = self.savedDataArr["measPos"],
                nomPos = self.savedDataArr["fitPos"],
                doRaise = False,
            )
            quadrupoleMag, quadrupoleAng = fitQuadrupole.model.getMagnitudeAngle()
            quadrupoleResidRadRMS = fitData.arrayRMS(fitQuadrupole.getRadError())

            self.logWdg.addMsg("quadrupole fit for %5d plates: %65s %8.2f  %8.2f    %8.4f" % \
                (self.numSavedPlates, "", quadrupoleMag * 1.0e6, quadrupoleAng, quadrupoleResidRadRMS))
        except Exception, e:
            self.logWdg.addMsg("Failed to fit quadrupole for cumulative data: %s" % (RO.StringUtil.strFromException(e,)), severity=RO.Constants.sevError)


    def processFile(self, filePath):
        """Process one file of plug plate CMM measurements.
        """
        fileDir, fileName = os.path.split(filePath)

        dataArr, plateID, measDate = readFile(filePath)
        if plateID not in fileName:
            raise RuntimeError("File name = %s does not match plate ID = %s" % (fileName, plateID))

        # fit translation, rotation and scale; raise an exception if fit fails since we can't use the data
        fitTransRotScale = fitData.ModelFit(
            model = fitData.TransRotScaleModel(),
            measPos = dataArr["measPos"],
            nomPos = dataArr["nomPos"],
            doRaise=True,
        )
        xyOff, rotAngle, scale = fitTransRotScale.model.getTransRotScale()
        fitPos = fitTransRotScale.getFitPos()
        residPosErr = fitTransRotScale.getPosError()
        residRadErr = fitTransRotScale.getRadError()

        # Identify manga holes by their larger diameter.
        mangaNomDia = 3.28
        largeDia = dataArr["measDia"] > 3.
        mangaInds = numpy.flatnonzero(largeDia)
        nonMangaInds = numpy.flatnonzero(numpy.invert(largeDia))
        ## ad-hoc adjust nominal diameter for manga holes.
        dataArr["nomDia"][mangaInds] = mangaNomDia
        diaErr = dataArr["measDia"] - dataArr["nomDia"]


        newDataToSave = numpy.zeros(dataArr.shape, dtype=self.savedDataArr.dtype)
        newDataToSave["measPos"] = dataArr["measPos"]
        newDataToSave["fitPos"] = fitPos
        self.savedDataArr = numpy.concatenate((self.savedDataArr, newDataToSave))
        self.numSavedPlates += 1

        # handle quadrupole (if it cannot be fit then display what we already got)
        fitQuadrupole = fitData.ModelFit(
            model = fitData.QuadrupoleModel([0.01, 0.0]),
            measPos = dataArr["measPos"],
            nomPos = fitTransRotScale.getFitPos(),
            doRaise = False,
        )
        quadrupoleMag, quadrupoleAng = fitQuadrupole.model.getMagnitudeAngle()
        quadrupoleResidPosErr = fitQuadrupole.getPosError()
        quadrupleResidRadErr = fitQuadrupole.getRadError()

        if len(fileName) > 9:
            dispFileName = fileName + "\n         "
        else:
            dispFileName = fileName

        residRadErrRMS_nonManga = fitData.arrayRMS(residRadErr[nonMangaInds])
        diaErrRMS_nonManga = fitData.arrayRMS(diaErr[nonMangaInds])
        residRadErrRMS_manga = fitData.arrayRMS(residRadErr[mangaInds])
        diaErrRMS_manga = fitData.arrayRMS(diaErr[mangaInds])
        maxDiaErr = numpy.max(diaErr)
        quadrupleResidRadErrRMS = fitData.arrayRMS(quadrupleResidRadErr)
        self.logWdg.addMsg("%-9s %10s %5d  %8.3f  %8.3f   %8.6f  %8.3f %8.4f  %8.4f   %8.4f      %8.4f     %8.3f   %8.2f  %8.2f    %8.4f" % \
            (dispFileName, measDate, len(dataArr), xyOff[0], xyOff[1], scale, rotAngle, \
            residRadErrRMS_nonManga, diaErrRMS_nonManga, residRadErrRMS_manga, diaErrRMS_manga, maxDiaErr, \
            quadrupoleMag * 1.0e6, quadrupoleAng, quadrupleResidRadErrRMS))

        posErrDType = [
            ("nomPos", float, (2,)),
            ("residPosErr", float, (2,)),
            ("quadrupoleResidPosErr", float, (2,)),
        ]
        posErrArr = numpy.zeros(dataArr.shape, dtype=posErrDType)
        posErrArr["nomPos"] = dataArr["nomPos"]
        posErrArr["residPosErr"] = residPosErr
        posErrArr["quadrupoleResidPosErr"] = quadrupoleResidPosErr

        self.graphWdg.addData(fileName, posErrArr)

        # save fit data to a file
        inBriefName = fileName
        if inBriefName.startswith("D"):
            inBriefName = inBriefName[1:]
        outName = "fit_%s_%s.txt" % (inBriefName, measDate)
        outPath = os.path.join(fileDir, outName)
        with file(outPath, "w") as outFile:
            outFile.write("PlugPlate %s\nMeasDate %s\n" % (plateID, measDate))
            outFile.write("# Lengths in mm, angles in deg\n")
            outFile.write("FitOffset %8.3f %8.3f\nFitScale %8.6f\nFitRotAngle %8.3f\n" % \
                (xyOff[0], xyOff[1], scale, rotAngle))
            outFile.write("QPMagnitude %10.8f\nQPAngle %8.2f\n" % (quadrupoleMag, quadrupoleAng))
            outFile.write("ResidRadErrRMS (MANGA) %10.4f (%10.4f)\nDiaErrRMS (MANGA) %8.4f (%8.4f)\nMaxDiaErr %8.3f\nQPResidRadErrRMS %8.4f\n" % \
                (residRadErrRMS_nonManga, residRadErrRMS_manga, diaErrRMS_nonManga, diaErrRMS_manga, maxDiaErr, quadrupleResidRadErrRMS))
            outFile.write("DataTable\n")
            outFile.write("    NomX     NomY    MeasX    MeasY   ResidX   ResidY ResidRad   NomDia   DiaErr QPResidX QPResidY QPResidRad\n")
#                         |        |        |        |        |        |        |        |        |        |        |        |          |
            for ind, dataRow in enumerate(dataArr):
                outFile.write("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.3f\n" % (
                    dataRow[ "nomPos"][0], dataRow[ "nomPos"][1],
                    dataRow["measPos"][0], dataRow["measPos"][1],
                    residPosErr[ind][0], residPosErr[ind][1],
                    residRadErr[ind],
                    dataRow["nomDia"],
                    diaErr[ind],
                    quadrupoleResidPosErr[ind][0], quadrupoleResidPosErr[ind][1],
                    quadrupleResidRadErr[ind],
                ))

# match the Plug Plate: plateID
plateIDRE = re.compile(r"^Plug Plate: ([0-9a-zA-Z_]+) *(?:#.*)?$", re.IGNORECASE)
# matching Date: YYYY-MM-DD
measDateRE = re.compile(r"^Date: ([0-9-]+) *(?:#.*)?$", re.IGNORECASE)

def readHeader(filePath):
    """Read header information from a measurement file

    Inputs:
    - filePath: path to data file (see format below).

    Return:
    - nHeaderLines: number of lines of header
    - plateID (a string)
    - measDate (a string)

    The header includes the following (by example), plus possible comments starting with #
    Plug Plate: 2247
    Date: 2005-06-23

    Meas X    Meas Y     Nom X     Nom Y     Err X     Err Y  Meas D   Nom D   Round
    """
    plateID = None
    measDate = None
    lineNum = 0
    with file(filePath, "rU") as dataFile:
        for line in dataFile:
            lineNum += 1
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            plateIDMatch = plateIDRE.match(line)
            if plateIDMatch:
                plateID = plateIDMatch.group(1)
                continue

            measDateMatch = measDateRE.match(line)
            if measDateMatch:
                measDate = measDateMatch.group(1)
                continue

            if line.lower().startswith("meas x"):
                if None in (plateID, measDate):
                    raise RuntimeError("Could not parse header")
                return (lineNum, plateID, measDate)

    raise RuntimeError("Could not parse header")

def readFile(filePath):
    """Read a measurement file

    Inputs:
    - filePath: path to data file (see format below).

    Returns:
    - dataArr: a structured array with fields:
        - measPos (x,y)
        - nomPos (x,y)
        - measDia
        - nomDia
        - measRoundness
    - plateID
    - measDate

    The file format is 4 lines of header followed by space-separated floating-point data.
    The measurements are in mm.
    Here is an example header:
    Plug Plate: 2247
    Date: 2005-06-23

    Meas X    Meas Y     Nom X     Nom Y     Err X     Err Y  Meas D   Nom D   Round
    """
    numHeaderLines, plateID, measDate = readHeader(filePath)
    inDtype = [
        ("measPosX", float),
        ("measPosY", float),
        ("nomPosX", float),
        ("nomPosY", float),
        ("posErrX", float),
        ("posErrY", float),
        ("measDia", float),
        ("nomDia", float),
        ("measRoundness", float),
    ]
    dataArr = numpy.loadtxt(filePath, dtype=inDtype, skiprows=numHeaderLines)
    desDType = [
        ("measPos", float, (2,)),
        ("nomPos", float, (2,)),
        ("measDia", float),
        ("nomDia", float),
        ("measRoundness", float),
    ]
    outArr = numpy.zeros(dataArr.shape, dtype=desDType)
    outArr["measPos"][:,0] = dataArr["measPosX"]
    outArr["measPos"][:,1] = dataArr["measPosY"]
    outArr["nomPos"][:,0] = dataArr["nomPosX"]
    outArr["nomPos"][:,1] = dataArr["nomPosY"]
    outArr["measDia"] = dataArr["measDia"]
    outArr["nomDia"] = dataArr["nomDia"]
    outArr["measRoundness"] = dataArr["measRoundness"]
    return outArr, plateID, measDate
