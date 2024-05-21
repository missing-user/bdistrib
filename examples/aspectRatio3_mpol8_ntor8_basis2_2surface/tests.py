#!/usr/bin/env python

# This python script checks the output file for an example to 
# see if the results are close to expected values.  This script may be
# run directly, and it is also called when "make test" is run from the
# main BDISTRIB directory.

import sys

sys.path = ["../"] + sys.path
from testsCommon import *

desiredTolerance = 0.03

numFailures = 0

f = readOutputFile()

variableName = 'svd_s_transferMatrix'
data = f.variables[variableName][()]

# These numbers come from the basis_option=3 case
numFailures += arrayShouldBe(variableName, data[0,:], [\
    0.607490335213963, 0.589192794088282, 0.542536334050193, 0.506806588519051, \
    0.44969004300756, 0.421239961356771, 0.366379203520474,                   \
    0.356139843449531, 0.339653064334853, 0.337884534065549,                  \
    0.323331043328159, 0.304949710632313, 0.303763176927845,                  \
    0.302275984828589, 0.281557489896902, 0.26547778259879,                   \
    0.259746212290225, 0.252131370836729, 0.246352424562248,                  \
    0.232035845459978, 0.22711477828712, 0.222139657894696,\
    0.212945659529801, 0.211710030928354, 0.200517306665367,\
    0.194091595983469, 0.193002934838878, 0.19276386093305,\
    0.189876665657361, 0.189295702076502, 0.178499591403917,\
    0.176309725923651, 0.175862089899159, 0.172983808229027,\
    0.163123259361317, 0.156019550965374, 0.151237007076328,\
    0.15089698157323, 0.14640717200485, 0.140841016182582, 0.137944047942987,\
    0.129857448150553, 0.124614049622635, 0.122297370992557,\
    0.118322140642292, 0.11667763071986, 0.115561501592206,\
    0.115541809320068, 0.108735031852018, 0.107736100281378,\
    0.107666091096947, 0.106849116086679, 0.102193094508361,\
    0.101362864181465, 0.096054811407493, 0.095718858742488,\
    0.0913824272101467, 0.0855127374624614, 0.0829313737562351,\
    0.0825729854572384, 0.0821661621001919, 0.0806149750099409,\
    0.0716973300407983, 0.0697452533949301, 0.068481772367013,\
    0.0684766517446092, 0.0684760471521765, 0.0679727583467491,\
    0.0666672446383673, 0.0648334613561275, 0.0648311654853511,\
    0.0632662537191218, 0.0592642730667506, 0.0592529042062994,\
    0.056939679871764, 0.0552586229206479, 0.0548211789890563,\
    0.0523793053135483, 0.0523204426419841, 0.0477730674352245,\
    0.0453999067870541, 0.044896643885103, 0.0448358597333031,\
    0.0445560807217712, 0.0410800303639223, 0.0404724225783292,\
    0.0404723961077637, 0.038705079402562, 0.0387049843592777,\
    0.0372965373402629, 0.0371101845447544, 0.0367442622367256,\
    0.0362419164269342, 0.0359634635744499, 0.0359630495185585,\
    0.0325027843787135, 0.0325006985078774, 0.0304057816109254,\
    0.0298752131280538, 0.0285978064998139, 0.0285868555786331,\
    0.0279297587416107, 0.0246118165567521, 0.0245099383399632,\
    0.0244565464377698, 0.0241591323967747, 0.0238665219965405,\
    0.0238665144835234, 0.0229860734139177, 0.0229860523816096,\
    0.0216093152307436, 0.0216092829113719, 0.02119655850403,\
    0.0204857169602324, 0.0202654237022347, 0.0198579804844702,\
    0.0198579556764902, 0.0178539291211975, 0.0178536784471443,\
    0.0167739484242147, 0.0161047019943563, 0.0156977203582643,\
    0.0156976950086827, 0.0140612638953482, 0.0139343819282577,\
    0.0139343726869448, 0.013492490364908, 0.0134924590779882,\
    0.013462441230402, 0.0134557993059233, 0.0127848307320213,\
    0.0127847516587829, 0.0118833129083892, 0.0118831943671912,\
    0.0112046329698903, 0.011154866731615], desiredTolerance, requireSameLength=False)


#numFailures += shouldBe(variableName, data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(variableName, data[0,0], 0.500, desiredTolerance)
#numFailures += arrayShouldBe(variableName, data[0,:], 0.500, desiredTolerance, requireSameLength=False)

f.close()
print "numFailures: ",numFailures
exit(numFailures > 0)
