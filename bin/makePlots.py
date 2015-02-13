import pickle
import datetime
import numpy
import numpy.random
import matplotlib
import itertools
matplotlib.use('Agg')
# matplotlib.use("PDF")
import seaborn
import matplotlib.pyplot as plt

def sampleFromHist(nSamples, histogram):
        """nSamples: number of samples to return
            histogram: output of matplotlib.hist
        """
        hist, edges, foo = histogram
        # determine bin width
        width = edges[1]-edges[0]
        # create array of bin centers
        centers = edges[:-1] + width/2.
        prob = hist*width
        return numpy.random.choice(centers, size=nSamples, p=prob)

def sampleRandomAngle(nSamples):
    """from 0 to 2*pi
    """
    return numpy.random.random(nSamples)*2*numpy.pi

def sampleFiberConc(nSamples):
    """Fiber core concentricity
    """
    return numpy.random.random(nSamples)*0.003

def sampleFerruleID(nSamples):
    """flat between 0.192 and 0.202
    """
    return (.202 - .192)*numpy.random.random(nSamples) + .192

def sampleFerruleConc(nSamples):
    return numpy.random.random(nSamples)*0.009

def sampleFerruleOD(nSamples):
    return (2.157 - 2.151)*numpy.random.random(nSamples) + 2.151

class PlotManager(object):
    def __init__(self, picklePath, figsize=(8,4), figFormat="tiff", radErrLim=((0, .02))):
        self.figsize = figsize
        self.figFormat = figFormat
        self.data = pickle.load(open(picklePath, "rb"))
        self.radErrLim = radErrLim
        self.cmap = seaborn.cubehelix_palette(start=8, light=1, as_cmap=True)
        self.figDir="figs"

    def plotHist(self, plotName, errType="radErr", rms=False):
        # get only values between dates
        xKey = "date"
        yKey = errType
        if rms:
            xKey += "RMS"
            yKey += "RMS"
        yearSpan1=[2005,2008]
        yearSpan2=[2013,2016]
        if "rad" in errType:
            bins=50
        else:
            bins=40
        plt.figure()
        histOut = []
        for yearSpan in [yearSpan1, yearSpan2]:
            keepInds = numpy.nonzero(numpy.bitwise_and(self.data[xKey]>yearSpan[0], self.data[xKey]<yearSpan[1]))
            errs = self.data[yKey][keepInds]
            histOut.append(plt.hist(errs, bins=bins, alpha=0.5, normed=True))
        plt.legend(["2005-2008", "2013+"])
        plt.title(plotName)
        plt.xlabel("Err (mm)")
        plt.ylabel("Normalized Counts")
        plt.savefig(self.figDir+"/"+plotName+"."+self.figFormat, format=self.figFormat)
        return histOut


    def plotScatter(self, plotName, errType="radErr", rms=False, alpha=0.1):
        plt.figure()
        xKey = "date"
        yKey = errType
        if rms:
            xKey += "RMS"
            yKey += "RMS"
        plt.plot(self.data[xKey], self.data[yKey], ".k", alpha=alpha)
        plt.title(plotName)
        if "rad" in errType:
            plt.ylim(self.radErrLim)
        plt.xlabel("year")
        plt.ylabel("Err (mm)")
        plt.savefig(self.figDir+"/"+plotName+"."+self.figFormat, format=self.figFormat)

    def plot2DHist(self, plotName, errType="radErr", rms=False):
        xKey = "date"
        yKey = errType
        if rms:
            xKey += "RMS"
            yKey += "RMS"
        # create histogram, normalize by # of measurements
        H,yedges,xedges=numpy.histogram2d(self.data[yKey], self.data[xKey],bins=(150,50))
        # throw away all above error lim
        # import pdb; pdb.set_trace()
        colSums = numpy.sum(H, axis=0)
        # replace any zeros with 1 to avoid division error
        colSums[numpy.nonzero(colSums == 0)] = 1
        colSums.shape = (1, len(colSums))
        divisor = numpy.dot(numpy.ones((H.shape[0],1)),colSums)
        # print "divisor shape, ", divisor.shape
        # print "divisor", divisor
        normed = H / divisor
        # print numpy.sum(normed, axis=0)
        # plt.imshow(normed)
        # plt.savefig(self.figDir+"/"+"/home/csayres/imshow.tiff", format="tiff")
        X,Y = numpy.meshgrid(xedges,yedges)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.pcolormesh(X,Y,normed,cmap=self.cmap)
        if "rad" in errType:
            ax.set_ylim(self.radErrLim)
        ax.set_title(plotName)
        ax.set_xlabel("year")
        ax.set_ylabel("Err (mm)")
        # ax.set_aspect("equal")
        plt.savefig(self.figDir+"/"+plotName+"."+self.figFormat, format=self.figFormat)

    def fullErrorHist(self, plotName, diaHists, radHists):
        fiberOD = 0.190 # perfect
        nomHoleDia = 2.167
        nSamples = 100000
        plt.figure()
        for diaHist, radHist in itertools.izip(diaHists, radHists):
            holePosErr = [sampleFromHist(nSamples, radHist), sampleRandomAngle(nSamples)]
            # note ferrule always gets pushed to one edge of the hole
            ferrulePosErr = [
                (sampleFromHist(nSamples, diaHist)+nomHoleDia-sampleFerruleOD(nSamples))/2., # nom hole dia needed because diaHist is error not diameter
                sampleRandomAngle(nSamples)
            ]
            ferruleConcErr = [sampleFerruleConc(nSamples), sampleRandomAngle(nSamples)]
            # ferrulePosErr = [[0]*nSamples, sampleRandomAngle(nSamples)]
            # fiber errors
            fiberPosErr = [
                (sampleFerruleID(nSamples)-fiberOD)/2.,
                sampleRandomAngle(nSamples)
            ]
            fiberConcErr = [sampleFiberConc(nSamples), sampleRandomAngle(nSamples)]

            # convert all to cartesian and sum errors
            xyErr = numpy.zeros((nSamples, 2))
            for err in [holePosErr, ferrulePosErr, ferruleConcErr, fiberPosErr, fiberConcErr]:
                err = numpy.asarray(err).T
                # convert from r, theta to x, y
                xyErr += numpy.asarray([err[:,0]*numpy.cos(err[:,1]), err[:,0]*numpy.sin(err[:,1])]).T
            # determine final radial error
            fullRadialErr = numpy.linalg.norm(xyErr, axis=1)
            plt.hist(fullRadialErr, bins=50, alpha=0.5, normed=True)
        plt.legend(["2005-2008", "2013+"])
        plt.title(plotName)
        plt.xlabel("Err (mm)")
        plt.ylabel("Normalized Counts")
        plt.savefig(self.figDir+"/"+plotName+"."+self.figFormat, format=self.figFormat)

def makeSubsetPickle():
    plateMeasurements = pickle.load(open("plateMeasurements.p", "rb"))
    date = []
    dateRMS = []
    radErr = []
    radErrRMS = []
    diaErr = []
    diaErrRMS = []
    minDate = datetime.datetime(year=2004, month=1, day=1)
    minYear = minDate.year
    for meas in plateMeasurements:
        if meas["measDate"] == "None":
            continue
        _radErr = numpy.asarray(meas["radErr"])
        _diaErr = numpy.asarray(meas["diaErr"])
        if len(_radErr) != len(_diaErr):
            continue
        # only keep measurements for rad errors < 40 microns
        keepInds = numpy.nonzero(numpy.bitwise_and(_radErr<0.04, numpy.abs(_diaErr)<0.04))
        _radErr = _radErr[keepInds]
        _diaErr = _diaErr[keepInds]
        radErr.extend(list(_radErr))
        diaErr.extend(list(_diaErr))
        _date = datetime.datetime.strptime(meas["measDate"], "%Y-%m-%d")
        decimalDate = (_date - minDate).total_seconds() / float(3.15569*10**7) + minYear
        nMeas = len(_radErr)
        dateRMS.append(decimalDate)
        date.extend([decimalDate]*nMeas)
        diaErrRMS.append(meas["diaErrRMS"])
        radErrRMS.append(meas["residRadErrRMS"])
    outDict = {
        "date": numpy.asarray(date),
        "dateRMS": numpy.asarray(dateRMS),
        "radErr": numpy.asarray(radErr),
        "radErrRMS": numpy.asarray(radErrRMS),
        "diaErr": numpy.asarray(diaErr),
        "diaErrRMS": numpy.asarray(diaErrRMS)
    }
    pickle.dump(outDict, open("plateMeasurementsSubset.p", "wb"))


if __name__ == "__main__":
    # makeSubsetPickle()
    pm = PlotManager("plateMeasurementsSubset.p")
    pm.plotScatter("radRMSErr", errType="radErr", rms=True, alpha=0.1)
    pm.plotScatter("radErr", errType="radErr", rms=False, alpha=0.005)
    pm.plotScatter("diaRMSErr", errType="diaErr", rms=True, alpha=0.1)
    pm.plotScatter("diaErr", errType="diaErr", rms=False, alpha=0.005)

    diaHists = pm.plotHist("diaErr_1DHist", errType="diaErr", rms=False)
    # pm.plotHist("diaRMSErr_1DHist", errType="diaErr", rms=True)
    radHists = pm.plotHist("radErr_1DHist", errType="radErr", rms=False)
    pm.plotHist("radRMSErr_1DHist", errType="radErr", rms=True)

    pm.plot2DHist("radRMSErr_2DHist", errType="radErr", rms=True)
    pm.plot2DHist("radErr_2DHist", errType="radErr", rms=False)

    # pm.plot2DHist("diaRMSErr_2DHist", errType="diaErr", rms=True)
    pm.plot2DHist("diaErr_2DHist", errType="diaErr", rms=False)

    pm.fullErrorHist("includingAllErrSources", diaHists, radHists)

    # pm = PlotManager("plateRMSMeas.p")
    # pm.plotScatter("plateRMSErr", alpha=0.1)
    # pm.plot2DHist("plateRMSErr_2DHist")
    # pm.plotHist("plateRMSErr_1DHist")


    # pm = PlotManager("allMeas.p")
    # pm.plotScatter("allMeasErr", alpha=0.005)
    # pm.plot2DHist("allMeasErr_2DHist")
    # pm.plotHist("allMeasErr_1DHist")
