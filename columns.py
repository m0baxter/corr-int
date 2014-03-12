
def colum_max( arr, i ):
	"""Finds the length of the longest element in the given column of arr."""
	
	longest = 0
	
	for j in range( len(arr) ):
		
		temp = len(arr[j][i])

		if temp > longest:
			longest = temp
			
	return longest

	
def column_max_array( arr ):
	"""Returns an array of the lenghts of the longest elements in each column or arr."""
	
	long_arr = []
	
	for i in range( len(arr[0]) ):
		
		long_arr.append( colum_max(arr,i) )
	
	return long_arr


def form_columns( readfile, writefile ):
	"""Formats the file into columns and stores it in writefile."""
	
	#read the file and make it into an array:
	full = readfile.read()
	arr=[row.split() for row in full.split("\n")]
	del arr[-1]
	
	#Make a list of the longest word in each column:
	longest = column_max_array( arr )
	
	cols = len(longest)
	rows = len(arr)
	
	for i in range( rows ):
		
		line = ""
		
		for j in range( cols - 1 ):
				
			line += arr[i][j] + (3 + longest[j] - len(arr[i][j]) ) * " "
			
		line += arr[i][cols - 1] + "\n"
		
		writefile.write(line)
		
	return

if __name__ == "__main__":
	
	readpath  = ".output/MCHF-un/z{0}/He2pHe_MCHF_resp_E{1}_z{0}.txt"
	writepath = ".output/MCHF/z{0}/He2pHe_MCHF_resp_E{1}_z{0}.txt"
	
	readfile  = open( readpath, "r")
	writefile = open( writepath, "w")
				
	form_columns(readfile, writefile)
				
	readfile.close()
	writefile.close()

