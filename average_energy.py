import sys
import numpy as np

try:
    infile  = open(sys.argv[1],'r')
except:
    print("Provide a file name!")
    sys.exit()


lines = infile.readlines()	
infile.close()

data_ET = []
data_E = []

for line in lines:
    if line.startswith('///') or len(line) < 4 :
        continue
    dummy = line.split()
    data_ET.append(float(dummy[0]))
    data_E.append(float(dummy[1]))

print("Sum ET average: "  + str(np.average(data_ET)))
print("Sum ET std: "      + str(np.std(data_ET)))
print("Sum ET variance: " + str(np.var(data_ET)))
print("Sum E average: "   + str(np.average(data_E)))
print("Sum E std: "       + str(np.std(data_E)))
print("Sum E variance: "  + str(np.var(data_E)))