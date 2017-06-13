#!/usr/bin/env python
from __future__ import with_statement
""" Grab and process all historical measurement files

History:
2015-01-22 CS    Initial.
"""
import os
import glob
import time
import pickle

from fitPlugPlateMeas import PlateMeas

measDir = "/nfsmount/shopdc0/meas"
assert os.path.exists(measDir)
plateList = []

def recursiveDive(plateList, path):
    print("current dir: ", path)
    dFiles = glob.glob(os.path.join(path, "D905*"))
    for f in dFiles[1:]:
        try:
            plate = PlateMeas(f).export()
            if not plate["measDate"]:
                continue
            plateList.append(plate)
        except Exception:
            pass
            # print("Exceptoin!: %s"%str(e))
    # now recurse into lower directories
    directories = [x[0] for x in os.walk(path)]
    print("directories", directories)
    for name in directories[1:]:
        recursiveDive(plateList, os.path.join(path, name))
t1 = time.time()
recursiveDive(plateList, measDir)
print("took: ", time.time()-t1)
pickle.dump( plateList, open( "plateMeasurements.p", "wb" ) )
