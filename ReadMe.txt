Fitting plug plate measurements from the UW shop's CMM

Select the CMM measurement results you wish to process and drop them onto the application icon. It will then fit each set of data and display a report. In addition it shows a graph of the residuals (which includes a pop-up menu allowing you to select which plate you wish to display).

Once one batch has been processed you may drop additional measurement files to process more. However, the graph data keeps accumulating until you quite the application and eventually the application will run out of memory.



To build the application see BuildForMac/README.html



Data files from the CMM are named D<plate number><plate letter><measurement number>, e.g. D2130_1 or D2130A1. The format is as per the following example:

Plug Plate: 2130A
Date: 2005-03-15

    Meas X    Meas Y     Nom X     Nom Y     Err X     Err Y  Meas D   Nom D   Round
   -80.754  -215.282   -80.743  -215.281    -0.011    -0.001   2.176   2.167   0.004
...

Notes:
- Units are mm
- Meas = measured, Nom = nominal, Err = Meas - Nom
- D is diameter, Round is a measure of roundness
- The plug plate ID doesn't not include the measurement number.
- The files are written using Windows line endings (\r\n), but the data processing code should accept any standard line endings.
- The data processing software should accept blank lines and comment lines beginning with #
