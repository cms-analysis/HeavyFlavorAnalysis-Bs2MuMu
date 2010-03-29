#! /bin/tcsh

setenv MODE BsToMuMu      && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BdToMuMu      && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list

setenv MODE BsToKPi       && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BsToKK        && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BsToPiPi      && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BsToMuMuGa    && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list

setenv MODE BuTo3MuNu     && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list

setenv MODE BdToKPi       && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BdToPiPi      && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BdToMuMuPi0   && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list

setenv MODE LambdaBToPPi  && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE LambdaBToPK   && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list

setenv MODE BsToJPsiPhi   && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BdToJPsiKstar && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BdToJPsiKs    && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
setenv MODE BuToJPsiKplus && mkPyFiles -t ./validation-XXXX.py -n 99999 -s $MODE -f ../Winter10/files/$MODE.list
