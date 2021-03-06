import os
import numpy
import re

from . import fitData

# match the Plug Plate: plateID
plateIDRE = re.compile(r"^Plug Plate: ([0-9a-zA-Z_]+) *(?:#.*)?$", re.IGNORECASE)
# matching Date: YYYY-MM-DD
measDateRE = re.compile(r"^Date: ([0-9-]+) *(?:#.*)?$", re.IGNORECASE)

__all__ = ["PlateMeas"]

class PlateMeas(object):
    def __init__(self, pathToFile=None, fileName=None, fileLines=None):
        """inputs:
        pathToFile: path to a "D" file (output by the cmm measurement routine)
        fileName: name of the file
        fileLines: lines of file to be parsed
        """
        if not pathToFile:
            assert not None in [fileName, fileLines]
        desDType = [
            ("measPos", float, (2,)),
            ("fitPos", float, (2,)),
        ]
        self.savedDataArr = numpy.zeros(0, dtype=desDType)
        if pathToFile:
            fileDir, fileName = os.path.split(pathToFile)
            self.fileDir = fileDir
            self.fileName = fileName
            with open(pathToFile, "r") as f:
                fileLines = f.readlines()
        else:
            self.fileDir = None
            self.fileName = fileName
        self.processFile(fileLines)

    def processFile(self, fileLines):
        self.dataArr, plateID, measDate = readFile(fileLines)
        try:
            self.plateID = int(plateID)
        except ValueError:
            # must be string or something in it
            plateIDint=""
            for char in plateID:
                try:
                    int(char)
                    plateIDint += char
                except:
                    break
            self.plateID = int(plateIDint)
        self.measDate = None if measDate is None else str(measDate) # meas date may be None
        if plateID not in self.fileName:
            raise RuntimeError("File name = %s does not match plate ID = %s" % (self.fileName, self.plateID))

        # fit translation, rotation and scale; raise an exception if fit fails since we can't use the data
        fitTransRotScale = fitData.ModelFit(
            model = fitData.TransRotScaleModel(),
            measPos = self.dataArr["measPos"],
            nomPos = self.dataArr["nomPos"],
            doRaise=True,
        )
        self.xyOff, self.rotAngle, self.scale = fitTransRotScale.model.getTransRotScale()
        fitPos = fitTransRotScale.getFitPos()
        residPosErr = fitTransRotScale.getPosError()
        residRadErr = fitTransRotScale.getRadError()
        self.residPosErr = residPosErr
        self.residRadErr = residRadErr

        # Identify manga holes by their larger diameter.
        mangaNomDia = 3.28
        largeDia = self.dataArr["measDia"] > 3.
        mangaInds = numpy.flatnonzero(largeDia)
        nonMangaInds = numpy.flatnonzero(numpy.invert(largeDia))
        ## ad-hoc adjust nominal diameter for manga holes.
        self.dataArr["nomDia"][mangaInds] = mangaNomDia
        self.measDia = self.dataArr["measDia"]
        self.diaErr = self.dataArr["measDia"] - self.dataArr["nomDia"]


        newDataToSave = numpy.zeros(self.dataArr.shape, dtype=self.savedDataArr.dtype)
        newDataToSave["measPos"] = self.dataArr["measPos"]
        newDataToSave["fitPos"] = fitPos
        self.savedDataArr = numpy.concatenate((self.savedDataArr, newDataToSave))
        # self.numSavedPlates += 1

        # handle quadrupole (if it cannot be fit then display what we already got)
        fitQuadrupole = fitData.ModelFit(
            model = fitData.QuadrupoleModel([0.01, 0.0]),
            measPos = self.dataArr["measPos"],
            nomPos = fitTransRotScale.getFitPos(),
            doRaise = False,
        )
        self.quadrupoleMag, self.quadrupoleAng = fitQuadrupole.model.getMagnitudeAngle()
        quadrupoleResidPosErr = fitQuadrupole.getPosError()
        self.quadrupoleResidPosErr = quadrupoleResidPosErr
        quadrupleResidRadErr = fitQuadrupole.getRadError()
        self.quadrupleResidRadErr = quadrupleResidRadErr

        self.residRadErrRMS_nonManga = fitData.arrayRMS(residRadErr[nonMangaInds])
        self.diaErr_nonManga = self.diaErr[nonMangaInds]
        self.diaErrRMS_nonManga = fitData.arrayRMS(self.diaErr_nonManga)
        self.diaErr_manga = self.diaErr[mangaInds]
        self.residRadErrRMS_manga = fitData.arrayRMS(residRadErr[mangaInds])
        self.diaErrRMS_manga = fitData.arrayRMS(self.diaErr_manga)
        self.maxDiaErr = numpy.max(self.diaErr)
        self.quadrupleResidRadErrRMS = fitData.arrayRMS(quadrupleResidRadErr)

        posErrDType = [
            ("nomPos", float, (2,)),
            ("residPosErr", float, (2,)),
            ("quadrupoleResidPosErr", float, (2,)),
        ]
        self.posErrArr = numpy.zeros(self.dataArr.shape, dtype=posErrDType)
        self.posErrArr["nomPos"] = self.dataArr["nomPos"]
        self.posErrArr["residPosErr"] = residPosErr
        self.posErrArr["quadrupoleResidPosErr"] = quadrupoleResidPosErr

    def export(self):
        xPos = []
        yPos = []
        measx = self.dataArr["measPos"][:,0]
        measy = self.dataArr["measPos"][:,1]
        xErr = []
        yErr = []
        radErr = []
        qpXErr = []
        qpYErr = []
        qpRadErr = []
        for meas in self.posErrArr:
            xPos.append(meas[0][0])
            yPos.append(meas[0][1])
            xErr.append(meas[1][0])
            yErr.append(meas[1][1])
            radErr.append((meas[1][0]**2+meas[1][1]**2)**0.5)
            qpXErr.append(meas[2][0])
            qpYErr.append(meas[2][1])
            qpRadErr.append((meas[2][0]**2+meas[2][1]**2)**0.5)
        return {
            "fileName": self.fileName,
            "plateID": self.plateID,
            "measDate": self.measDate,
            "xyOff": self.xyOff,
            "rotAngle": self.rotAngle,
            "scale": self.scale,
            "quadrupoleMag": self.quadrupoleMag,
            "quadrupoleAngle": self.quadrupoleAng,
            "residRadErrRMS": self.residRadErrRMS_nonManga,
            "diaErrAll": list(self.diaErr),
            "diaErr": self.diaErr_nonManga,
            "diaErrRMS": self.diaErrRMS_nonManga,
            "residRadErrRMS_manga": self.residRadErrRMS_manga,
            "diaErrRMS_manga": self.diaErrRMS_manga,
            "maxDiaErr": self.maxDiaErr,
            "quadrupleResidRadErrRMS": self.quadrupleResidRadErrRMS,
            "measDia": self.measDia,
            "nomDia": list(self.dataArr["nomDia"]),
            "xPos": xPos,
            "yPos": yPos,
            "measx": measx,
            "measy": measy,
            "xErr": xErr,
            "yErr": yErr,
            "radErr": radErr,
            "qpXErr": qpXErr,
            "qpYErr": qpYErr,
            "qpRadErr": qpRadErr,
        }


def readHeader(fileLines):
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
    for line in fileLines:
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
            if plateID is None:
                raise RuntimeError("Could not parse plateID in header")
            return (lineNum, plateID, measDate)

    raise RuntimeError("Could not parse header")

def readFile(fileLines):
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
    numHeaderLines, plateID, measDate = readHeader(fileLines)
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
    dataArr = numpy.loadtxt(fileLines, dtype=inDtype, skiprows=numHeaderLines)
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