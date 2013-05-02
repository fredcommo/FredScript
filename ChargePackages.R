
setwd("path_old_library")
pack.list <- list.files()
for(i in 1:length(pack.list)){
	install.packages(pack.list[i])
	cat(pack.list[i], "installed\n")
	}

