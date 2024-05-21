#!/usr/bin/env python

# This python script checks the output file for an example to
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main BDISTRIB directory.

import sys

sys.path = ["../"] + sys.path
from testsCommon import *

desiredTolerance = 0.001

numFailures = 0

f = readOutputFile()

variableName = "svd_s_transferMatrix"
data = f.variables[variableName][()]

# Compare to analytically expected values:
analyticalResults = [1.0 for i in range(84)]
numFailures += arrayShouldBe(
    variableName, data[0, :], analyticalResults, desiredTolerance
)

f.close()
exit(numFailures > 0)
