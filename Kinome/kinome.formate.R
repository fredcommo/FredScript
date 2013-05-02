kinome.formate<-function(X,File.infos=infos)

# sort values and build a data.frame in a same format as on the slide

# input : col datas 954 values

# output: 384 rows with well, kinase info, duplicates, mean, sd


{
	# Call of rename1.R and rename2.R

	source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.rename1.R")
	source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.rename2.R")
	
	X.rename<-kinome.rename1(X)
	output<-kinome.rename2(X.rename,File.infos)
	
	output
}