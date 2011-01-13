#!/usr/bin/env python
from __future__ import with_statement
"""Fit dipole in plug plate measurements

History:
2010-06-17 ROwen    2.0: first version using Python instead of Igor
2011-01-13 ROwen    2.1: Show another digit of radial position error.
                    Display radial position error as Pos Err instead of Rad Err.
                    Fit and show dipole moment as additional data.
"""
import math
import os.path
import re
import sys
import Tkinter
import numpy
import scipy.optimize
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import RO.Constants
import RO.Wdg

__version__ = "2.1"

plateIDRE = re.compile(r"^Plug Plate: ([0-9a-zA-Z_]+) *(?:#.*)?$", re.IGNORECASE)
measDateRE = re.compile(r"^Date: ([0-9-]+) *(?:#.*)?$", re.IGNORECASE)

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
    
        self.plotFig = matplotlib.figure.Figure(figsize=(7, 6), frameon=True)
        self.figCanvas = FigureCanvasTkAgg(self.plotFig, self)
        figCnvWdg = self.figCanvas.get_tk_widget()
        figCnvWdg.grid(row=1, column=0, columnspan=5, sticky="news")
        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(4, weight=1)
        self.plotAxis = None
        self.resetPlotAxis()

    def resetPlotAxis(self):
        if self.plotAxis:
            self.plotFig.delaxes(self.plotAxis)
        self.plotAxis = self.plotFig.add_subplot(
            1, 1, 1,
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
            "nomPos" and "residPosErr"
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
        self.plotAxis.quiverkey(q, 0.9, 1.02, 0.01, "0.01 mm")
        self.figCanvas.draw()
        

class FitPlugPlateMeasWdg(Tkinter.Frame):
    def __init__(self, master, filePathList):
        Tkinter.Frame.__init__(self, master)
        
        self.logWdg = RO.Wdg.LogWdg(
            master = self,
            width = 138,
            height = 50,
        )
        self.logWdg.text["font"] = "Courier 12"
        self.logWdg.grid(row=0, column=0, sticky="nsew")
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.logWdg.addOutput("Plug Plate Fitter %s\n" % (__version__,))
        
        self.graphTL = RO.Wdg.Toplevel(
            master = self,
            title = "Position Errors",
            geometry = "+950+50",
            closeMode = RO.Wdg.tl_CloseDisabled,
        )
        self.graphWdg = GraphWdg(
            master = self.graphTL,
        )
        self.graphWdg.pack(side="left", expand=True, fill="both")

        if RO.OS.PlatformName == "mac":
            self.tk.createcommand('::tk::mac::OpenDocument', self._macOpenDocument)
        
        if filePathList:
            self.processFileList(filePathList)
        else:
            self.logWdg.addOutput("Drop files on the application icon to process them\n")
    
    def processFile(self, filePath):
        fileName = os.path.basename(filePath)
        
        dataArr, plateID, measDate = readFile(filePath)
        if plateID not in fileName:
            raise self.logWdg.addOutput("File name = %s does not match plate ID = %s" % (fileName, plateID),
                severity = RO.Constants.sevError)
        coeffs = fitData(dataArr["measPos"], dataArr["nomPos"])
        # negate the offset and invert the scale to match results from the Igor-based analysis system
        xyOff = - coeffs[0:2]
        scale = 1.0 / math.sqrt(coeffs[2]**2 + coeffs[3]**2)
        rotAngle = math.atan2(coeffs[2], coeffs[3]) * 180.0 / math.pi
        
        fitPos = computeFitPos(coeffs, dataArr["measPos"])
        residPosErr = fitPos - dataArr["nomPos"]
        residRadPosErr = numpy.sqrt(residPosErr[:,0]**2 + residPosErr[:,1]**2)
    
        diaErr = dataArr["measDia"] - dataArr["nomDia"]
        
        # handle dipole
        dipoleCoeffs, nomAng = fitDipole(fitPos, dataArr["nomPos"], doRaise=False)
        dipoleFitPos = computeDipoleFitPos(dipoleCoeffs, fitPos, nomAng)
        dipoleResidPosErr = dipoleFitPos - dataArr["nomPos"]
        dipoleMag = dipoleCoeffs[0]
        dipoleAng = dipoleCoeffs[1] * 180.0 / math.pi
        
        self.logWdg.addOutput("%-8s  %10s %5d  %8.3f  %8.3f   %8.6f  %8.3f %8.4f  %8.4f  %8.4f   %8.3f   %8.3f     %8.4f\n" % \
            (fileName, measDate, len(dataArr), xyOff[0], xyOff[1], scale, rotAngle, \
            rms(residRadPosErr), rms(diaErr), numpy.max(diaErr), \
            dipoleMag, dipoleAng, rms(dipoleResidPosErr)))

        posErrDType = [
            ("nomPos", float, (2,)),
            ("residPosErr", float, (2,)),
        ]
        posErrArr = numpy.zeros(dataArr.shape, dtype=posErrDType)
        posErrArr["nomPos"] = dataArr["nomPos"]
        posErrArr["residPosErr"] = residPosErr

        self.graphWdg.addData(fileName, posErrArr)
    
    def processFileList(self, filePathList):
        self.logWdg.addOutput("""
File       Meas Date  Holes  Offset X  Offset Y   Scale     Rotation  Pos Err   Dia Err   Dia Err  Dipole Mag  Dipole Ang  DPl Pos Err
                                mm        mm                  deg     RMS mm    RMS mm    Max mm      1e-3        deg        RMS mm
""")
        for filePath in filePathList:
            self.processFile(filePath)

        self.logWdg.addOutput("\nReady for more files\n")

    def _macOpenDocument(self, *filePathList):
        """Handle Mac OpenDocument event
        """
        self.processFileList(filePathList)

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

def computeRadSqErr(coeffs, measPos, nomPos):
    """Compute the radial error squared for a particular set of coefficients
    
    Inputs:
    - coeffs: see computeFitPos
    - measPos: array of measured x,y positions
    - nomPos: array of nominal x,y positions
    """
    fitPos = computeFitPos(coeffs, measPos)
    residPosErr = fitPos - nomPos
    return residPosErr[:,0]**2 + residPosErr[:,1]**2

def computeFitPos(coeffs, measPos):
    """Compute the fit position for a given set of coefficients
    
    Inputs:
    - coeffs: coefficients for the model (see below)
    - measPos: array of measured x,y positions

    The model is: fit x  =  c0  +  meas x * (c3  -c2)
                      y     c1          y   (c2   c3)
    """
    rotMat = numpy.array(((coeffs[3], -coeffs[2]), (coeffs[2], coeffs[3])), dtype=float)
    return coeffs[0:2] + numpy.dot(measPos, rotMat)

def computeDipoleRadSqErr(coeffs, measPos, nomPos, nomAng):
    """Compute the radial error squared for a particular set of coefficients
    
    Inputs:
    - coeffs: see computeFitPos
    - measPos: array of measured x,y positions
    - nomPos: array of nominal x,y positions
    - nomAng: array of angles to nominal x,y positions = arctan2(nom y, nom x);
        precomputed to save time
    """
    fitPos = computeDipoleFitPos(coeffs, measPos, nomAng)
    residPosErr = fitPos - nomPos
    return residPosErr[:,0]**2 + residPosErr[:,1]**2

def computeDipoleFitPos(coeffs, measPos, nomAng):
    """Compute the fit position for a given set of coefficients
    
    Inputs:
    - coeffs: coefficients for the model (see below)
    - measPos: array of measured x,y positions
    - nomAng: array of angles to nominal position

    The model is: fit x  =  meas x * (1 + (c0 * 1e-3 * cos(2 * (nom ang - c1))))
                      y  =  meas y * (1 + (c0 * 1e-3 * cos(2 * (nom ang - c1))))
    """
    return measPos * (1.0 + (coeffs[0] * 1.0e-3 * numpy.cos(2.0 * (nomAng - coeffs[1]))))[:,numpy.newaxis]

def fitData(measPos, nomPos, doRaise=False):
    """Fit measured data to nominal data using the model described in computeFitPos
    
    Inputs:
    - measPos: array of measured x,y positions
    - nomPos: array of nominal x,y positions
    
    Fit the coefficients that minimize radial error squared = (fitX - nomX)**2 + (fitY - nomY)**2
    where fitX,Y is as computed by computeFitPos

    Returns coeffs
    """
    initialCoeffs = [0.0, 0.0, 0.0, 1.0]
    coeffs, status = scipy.optimize.leastsq(
        computeRadSqErr,
        initialCoeffs,
        args=(measPos, nomPos),
        maxfev = 5000,
    )
    if status not in range(5):
        if doRaise:
            raise RuntimeError("fit failed")
        else:
            coeffs[:] = numpy.nan
    return coeffs

def fitDipole(measPos, nomPos, doRaise=False):
    """Fit measured data to nominal data using the model described in computeDipoleFitPos
    
    Inputs:
    - measPos: array of measured x,y positions
    - nomPos: array of nominal x,y positions
    
    Fit the coefficients that minimize radial error squared = (fitX - nomX)**2 + (fitY - nomY)**2
    where fitX,Y is as computed by computeFitPos
    
    Returns:
    - coeffs
    - nomAng: array of angle to nominal position (to use in calling computeDipoleFitPos)
    """
    initialCoeffs = [0.0, 0.0]
    nomAng = numpy.arctan2(nomPos[:,1], nomPos[:,0])
    coeffs, status = scipy.optimize.leastsq(
        computeDipoleRadSqErr,
        initialCoeffs,
        args=(measPos, nomPos, nomAng),
        maxfev = 5000,
    )
    if status not in range(5):
        if doRaise:
            raise RuntimeError("fit failed")
        else:
            coeffs[:] = numpy.nan
    return coeffs, nomAng
    

def rms(arr):
    """Return the RMS of an array"""
    return math.sqrt(numpy.sum(arr**2) / len(arr))


if __name__ == "__main__":
    filePathList = sys.argv[1:]
    # strip first argument if it starts with "-", as happens when run as a Mac application
    if filePathList and filePathList[0].startswith("-"):
        filePathList = filePathList[1:]

    root = Tkinter.Tk()
    root.title("FitPlugPlateMeas")
    
    fitPlugPlateWdg = FitPlugPlateMeasWdg(master = root, filePathList=filePathList)
    fitPlugPlateWdg.pack(side="left", expand=True, fill="both")
    root.mainloop()
