import numpy as np
import mpmath as mp
import os as os
import columns as col

energy = [ 1, 2, 5, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 100, 200, 500, 1000, 2000 ]
Z      = [ 28, 29, 30, 31, 32, 33, 34 ]

def read_probs( path ):
	"""Reads the probability data from the file at path."""
	
	readfile = open(path, 'r')
	
	#read and discard header:
	readfile.readline()
	
	
	result = np.empty([18,71], dtype = mp.mpf)
	
	i = 0
	
	for line in readfile:
		
		data = line.split()
		
		for j in range(len(data)):
			
			result[j][i] = mp.mpf(data[j])
			
		i += 1
			
	readfile.close()
	
	return result

mp.mp.dps = 7

folder = ""

for E in energy:
		
	z28 = read_probs( ".output/" + folder + "/z28/" + folder + "-E{0}-z28.txt".format(E) )
	z29 = read_probs( ".output/" + folder + "/z29/" + folder + "-E{0}-z29.txt".format(E) )
	z30 = read_probs( ".output/" + folder + "/z30/" + folder + "-E{0}-z30.txt".format(E) )
	z31 = read_probs( ".output/" + folder + "/z31/" + folder + "-E{0}-z31.txt".format(E) )
	z32 = read_probs( ".output/" + folder + "/z32/" + folder + "-E{0}-z32.txt".format(E) )
	z33 = read_probs( ".output/" + folder + "/z33/" + folder + "-E{0}-z33.txt".format(E) )
	z34 = read_probs( ".output/" + folder + "/z34/" + folder + "-E{0}-z34.txt".format(E) )
		
	writefile = open( ".output/MCHF/temp", 'w' )
	
	writefile.write("b p_T p_P Ic_TT Ic_PP Ic_TP ie_pTT ie_pTI ie_pII ie_pTP ie_pIP ie_pPP p_TT p_TI p_II p_TP p_PI p_PP\n")
		
	for i in range(71):
			
		line = np.empty(18, dtype = mp.mpf)
			
		for j in range(18):
				
			line[j] = sum( [ z28[j][i], z29[j][i], z30[j][i], z31[j][i], z32[j][i], z33[j][i], z34[j][i] ] )/7
				
		writefile.write( "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17}\n".format(*line) )
			
	writefile.close()
		
	readfile = open( ".output/"+ folder + "/temp",'r')
	writefile = open( ".output/" + folder + "/" + folder + "-E{0}.txt".format(E), 'w' )
		
	col.form_columns(readfile, writefile)
				
	readfile.close()
	writefile.close()
		
	os.remove(".output/" + folder + "/temp")
			
