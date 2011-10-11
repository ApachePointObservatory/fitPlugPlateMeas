#!/usr/bin/env python
from __future__ import with_statement
"""Fit plug plate measurements

History:
2011-10-11 ROwen    Moved this code from a library module.
"""
import sys
import Tkinter
import fitPlugPlateMeas

filePathList = sys.argv[1:]
# strip first argument if it starts with "-", as happens when run as a Mac application
if filePathList and filePathList[0].startswith("-"):
    filePathList = filePathList[1:]

root = Tkinter.Tk()
root.title("FitPlugPlateMeas")

fitPlugPlateWdg = fitPlugPlateMeas.FitPlugPlateMeasWdg(master=root, filePathList=filePathList)
fitPlugPlateWdg.pack(side="left", expand=True, fill="both")
root.mainloop()
