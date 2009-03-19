#  This example program will convert a native format tipsy file
#  into a standard format file.
#

import sys
from pytipsy import ifTipsy, ofTipsy

inp = ifTipsy()      # The input file
out = ofTipsy()      # The output file

# Make sure we have two parameters, a native and a standard file.
if len(sys.argv) != 3:
    print >>sys.stderr, "Usage: %s <native> <standard>" % sys.argv[0]
    sys.exit(1)

# Open the native file and abort if there is an error.
inp.open(sys.argv[1],"native")
if not inp.is_open():
    print >>sys.stderr, "Unable to open Native binary %s" % argv[1]
    sys.exit(2)

# Open the output file.
out.open(sys.argv[2],"standard")
if not out.is_open():
    print >>sys.stderr, "Unable to create Standard binary %s" % sys.argv[2]
    sys.exit(2)

# Read the from the input and write it to the output.
for i in inp: out.write(i)

# Close the files.
out.close()
inp.close()

