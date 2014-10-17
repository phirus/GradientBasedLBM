import sys
from paraview.simple import *

print 'Number of files to be converted:', len(sys.argv) - 1

for x in range(1, len(sys.argv)):    
    inputFile = str(sys.argv[x])
    outputFile = inputFile[:-1] + 'u'
    print x,': Converting ', inputFile, '  ->  ', outputFile
    r = LegacyVTKReader( FileNames= str(sys.argv[x]) )
    w = XMLUnstructuredGridWriter()
    w.FileName = outputFile
    w.UpdatePipeline()