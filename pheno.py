import re

def main():

	f = open("Input_data/pheno_in.txt", "r") # open imput file 
	holder = [] # holds raw data
	
	for line in f: 
		holder.append(line)

	phenoData = [] # hold final data
	outF = open("Output/Data/pheno_out.txt", "w") # opens output file

	for entry in range(len(holder)): # cycles through each data entry
		items = holder[entry].strip().split("\t") # splits data by tab
		
		if re.search("_CC\t*$", holder[entry]) or re.search("^e", holder[entry]): # if header 
			lane = holder[entry] # sets lane to header name

		else:
			if re.search("^l", holder[entry]) == None: # if not col names 
				head = "" # for earlier hiseq run
				
				if re.search("^062118_CC", lane):
					head = "L6."
				elif re.search("^061518_CC", lane):
					head = "L3."
				elif re.search("^062518_CC", lane):
					head = "L2."

				items[5] = head + items[5] # forming header for counts compatbility
				print(items[5])
				items += [lane] # add new col entry for that row the desc lane
				
			else:
				items += ["lane\n"] # col name lane for new col
			outF.write(",".join(items)) # write row to outfile and join by comma

	f.close()
	outF.close()

main()