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
#import re
import Tkinter
import numpy
#import scipy.optimize
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import RO.Constants
import RO.StringUtil
import RO.Wdg
from . import fitData
from . import PlateMeas
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
        plateMeas = PlateMeas(filePath)
        # fileDir, fileName = os.path.split(filePath)

        # dataArr, plateID, measDate = readFile(filePath)
        # measDate = str(measDate) # meas date may be None
        # if plateID not in fileName:
        #     raise RuntimeError("File name = %s does not match plate ID = %s" % (fileName, plateID))

        # # fit translation, rotation and scale; raise an exception if fit fails since we can't use the data
        # fitTransRotScale = fitData.ModelFit(
        #     model = fitData.TransRotScaleModel(),
        #     measPos = dataArr["measPos"],
        #     nomPos = dataArr["nomPos"],
        #     doRaise=True,
        # )
        # xyOff, rotAngle, scale = fitTransRotScale.model.getTransRotScale()
        # fitPos = fitTransRotScale.getFitPos()
        # residPosErr = fitTransRotScale.getPosError()
        # residRadErr = fitTransRotScale.getRadError()

        # # Identify manga holes by their larger diameter.
        # mangaNomDia = 3.28
        # largeDia = dataArr["measDia"] > 3.
        # mangaInds = numpy.flatnonzero(largeDia)
        # nonMangaInds = numpy.flatnonzero(numpy.invert(largeDia))
        # ## ad-hoc adjust nominal diameter for manga holes.
        # dataArr["nomDia"][mangaInds] = mangaNomDia
        # diaErr = dataArr["measDia"] - dataArr["nomDia"]


        # newDataToSave = numpy.zeros(dataArr.shape, dtype=self.savedDataArr.dtype)
        # newDataToSave["measPos"] = dataArr["measPos"]
        # newDataToSave["fitPos"] = fitPos
        # self.savedDataArr = numpy.concatenate((self.savedDataArr, newDataToSave))
        # self.numSavedPlates += 1

        # # handle quadrupole (if it cannot be fit then display what we already got)
        # fitQuadrupole = fitData.ModelFit(
        #     model = fitData.QuadrupoleModel([0.01, 0.0]),
        #     measPos = dataArr["measPos"],
        #     nomPos = fitTransRotScale.getFitPos(),
        #     doRaise = False,
        # )
        # quadrupoleMag, quadrupoleAng = fitQuadrupole.model.getMagnitudeAngle()
        # quadrupoleResidPosErr = fitQuadrupole.getPosError()
        # quadrupleResidRadErr = fitQuadrupole.getRadError()

        if len(plateMeas.fileName) > 9:
            dispFileName = plateMeas.fileName + "\n         "
        else:
            dispFileName = plateMeas.fileName

        # residRadErrRMS_nonManga = fitData.arrayRMS(residRadErr[nonMangaInds])
        # diaErrRMS_nonManga = fitData.arrayRMS(diaErr[nonMangaInds])
        # residRadErrRMS_manga = fitData.arrayRMS(residRadErr[mangaInds])
        # diaErrRMS_manga = fitData.arrayRMS(diaErr[mangaInds])
        # maxDiaErr = numpy.max(diaErr)
        # quadrupleResidRadErrRMS = fitData.arrayRMS(quadrupleResidRadErr)
        self.logWdg.addMsg("%-9s %10s %5d  %8.3f  %8.3f   %8.6f  %8.3f %8.4f  %8.4f   %8.4f      %8.4f     %8.3f   %8.2f  %8.2f    %8.4f" % \
            (dispFileName, plateMeas.measDate, len(plateMeas.dataArr), plateMeas.xyOff[0], plateMeas.xyOff[1], plateMeas.scale, plateMeas.rotAngle, \
            plateMeas.residRadErrRMS_nonManga, plateMeas.diaErrRMS_nonManga, plateMeas.residRadErrRMS_manga, plateMeas.diaErrRMS_manga, plateMeas.maxDiaErr, \
            plateMeas.quadrupoleMag * 1.0e6, plateMeas.quadrupoleAng, plateMeas.quadrupleResidRadErrRMS))

        # posErrDType = [
        #     ("nomPos", float, (2,)),
        #     ("residPosErr", float, (2,)),
        #     ("quadrupoleResidPosErr", float, (2,)),
        # ]
        # posErrArr = numpy.zeros(dataArr.shape, dtype=posErrDType)
        # posErrArr["nomPos"] = dataArr["nomPos"]
        # posErrArr["residPosErr"] = residPosErr
        # posErrArr["quadrupoleResidPosErr"] = quadrupoleResidPosErr

        self.graphWdg.addData(plateMeas.fileName, plateMeas.posErrArr)

        # save fit data to a file
        inBriefName = plateMeas.fileName
        if inBriefName.startswith("D"):
            inBriefName = inBriefName[1:]
        outName = "fit_%s_%s.txt" % (inBriefName, plateMeas.measDate)
        outPath = os.path.join(plateMeas.fileDir, outName)
        with file(outPath, "w") as outFile:
            outFile.write("PlugPlate %s\nMeasDate %s\n" % (plateMeas.plateID, plateMeas.measDate))
            outFile.write("# Lengths in mm, angles in deg\n")
            outFile.write("FitOffset %8.3f %8.3f\nFitScale %8.6f\nFitRotAngle %8.3f\n" % \
                (plateMeas.xyOff[0], plateMeas.xyOff[1], plateMeas.scale, plateMeas.rotAngle))
            outFile.write("QPMagnitude %10.8f\nQPAngle %8.2f\n" % (plateMeas.quadrupoleMag, plateMeas.quadrupoleAng))
            outFile.write("ResidRadErrRMS (MANGA) %10.4f (%10.4f)\nDiaErrRMS (MANGA) %8.4f (%8.4f)\nMaxDiaErr %8.3f\nQPResidRadErrRMS %8.4f\n" % \
                (plateMeas.residRadErrRMS_nonManga, plateMeas.residRadErrRMS_manga, plateMeas.diaErrRMS_nonManga, plateMeas.diaErrRMS_manga, plateMeas.maxDiaErr, plateMeas.quadrupleResidRadErrRMS))
            outFile.write("DataTable\n")
            outFile.write("    NomX     NomY    MeasX    MeasY   ResidX   ResidY ResidRad   NomDia   DiaErr QPResidX QPResidY QPResidRad\n")
#                         |        |        |        |        |        |        |        |        |        |        |        |          |
            for ind, dataRow in enumerate(plateMeas.dataArr):
                outFile.write("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %10.3f\n" % (
                    dataRow[ "nomPos"][0], dataRow[ "nomPos"][1],
                    dataRow["measPos"][0], dataRow["measPos"][1],
                    plateMeas.residPosErr[ind][0], plateMeas.residPosErr[ind][1],
                    plateMeas.residRadErr[ind],
                    dataRow["nomDia"],
                    plateMeas.diaErr[ind],
                    plateMeas.quadrupoleResidPosErr[ind][0], plateMeas.quadrupoleResidPosErr[ind][1],
                    plateMeas.quadrupleResidRadErr[ind],
                ))

# # match the Plug Plate: plateID
# plateIDRE = re.compile(r"^Plug Plate: ([0-9a-zA-Z_]+) *(?:#.*)?$", re.IGNORECASE)
# # matching Date: YYYY-MM-DD
# measDateRE = re.compile(r"^Date: ([0-9-]+) *(?:#.*)?$", re.IGNORECASE)

