# Access to staging version
synapseRepoServiceEndpoint(endpoint="https://repo-staging.sagebase.org/repo/v1")
synapseAuthServiceEndpoint(endpoint="https://auth-staging.sagebase.org/auth/v1")
synapseLogin('frederic.commo@sagebase.org', 'Se@ttle7')

# Create a project
metProject <- Project(list(name = 'metPoject', description = 'Correlations between MET exprs. and ETs family genes'))
# Then create an entity
metProject <- creatEntity(metProject)
# Get the synapse Id
mainId = propertyValue(metProject, "id")  
# Create folders attached to the main entity 
# Doesn't work !!!!
metPlots<-Folder(list(name = "Correlation plots", parentId=myId, description = "Correlation tables as heatmaps"))
metTables<-Folder(list(name = "Correlation tables", parentId=myId, description = "Correlation tables"))
# Create new entities (use the web client)
metPlots <- createEntity(metPlots)
metTables <- createEntity(metTables)
# Get the entity Ids
plotsId = 'syn1586720'
tablesId = 'syn1586721'
# Push all the plots in the directory as data entities
setwd("/gluster/home/jcommo/MetProject/Plots")
listFiles = list.files()
for (l in 15:length(listFiles)){
  tmpName = listFiles[l]
  myData <- Data(list(name = tmpName, parentId = plotsId))
  myData <- addFile(myData, tmpName)
  myData <- storeEntity(myData)
}

# Push all the tables in the directory as data entities
setwd("/gluster/home/jcommo/MetProject/Tables")
listFiles = list.files()
for (l in 1:length(listFiles)){
  tmpName = listFiles[l]
  myData <- Data(list(name = tmpName, parentId = tablesId))
  myData <- addFile(myData, tmpName)
  myData <- storeEntity(myData)
}
