from __future__ import with_statement
"""The Python portion of the script that builds FitPlugPlateMeas

Usage:
% python setup.py [--quiet] py2app

History:
2010-06-16 ROwen
2010-06-24 ROwen    Modified for new package layout.
2011-10-11 ROwen    Modified for newer package layout.
                    Deleted unused inclModules and inclPackages.
"""
import os
import shutil
import subprocess
import sys
from setuptools import setup

# add various bits to the path (but preferably use eups instead)
pkgRoot = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
pythonRoot = os.path.join(pkgRoot, "python", "fitPlugPlateMeas")
pathList = [pythonRoot]
if os.path.isfile("extraPaths.txt"):
    with file("extraPaths.txt") as pathFile:
        for pathStr in pathFile:
            pathStr = pathStr.strip()
            if not pathStr or pathStr.startswith("#"):
                continue
            pathList.append(pathStr)
sys.path = pathList + sys.path

import fitPlugPlateMeas

appName = "FitPlugPlateMeas"
mainScript = "runFitPlugPlateMeas.py"
mainProg = os.path.join(pkgRoot, "bin", mainScript)
appPath = os.path.join("dist", appName + ".app")
iconFile = "%s.icns" % appName
versStr = fitPlugPlateMeas.__version__

# see plistlib for more info
plist = dict(
    CFBundleName                = appName,
    CFBundleExecutable          = appName,
    CFBundleShortVersionString  = versStr,
    CFBundleGetInfoString       = "%s %s" % (appName, versStr),
    CFBundleDocumentTypes       = [
        dict(
            CFBundleTypeName = "File",
            CFBundleTypeRole = "Editor",
            LSItemContentTypes = [
                "public.plain-text",
                "public.text",
                "public.data",
            ],
        ),
        dict(
            CFBundleTypeName = "Folder",
            CFBundleTypeRole = "Viewer",
            LSItemContentTypes = [
				"public.folder",
            ],
        ),
    ],
    LSPrefersPPC                = False,
)

setup(
    app = [mainProg],
    setup_requires = ["py2app"],
    options = dict(
        py2app = dict (
            plist = plist,
            iconfile = iconFile,
        )
    ),
)

print "*** Creating disk image ***"
appName = "%s_%s_Mac" % (appName, versStr)
destFile = os.path.join("dist", appName)
args=("hdiutil", "create", "-srcdir", appPath, destFile)
retCode = subprocess.call(args=args)
