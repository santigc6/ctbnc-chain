library("jsonlite")

for(i in list(3,4,5,6,10,15,20)){
	dataset_path =  paste("data/networks_and_trajectories_binary_data_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))
    
	dataset_path =  paste("data/networks_and_trajectories_ternary_data_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))

    
    
    dataset_path =  paste("data/networks_and_trajectories_binary_data_01_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))
    
	dataset_path =  paste("data/networks_and_trajectories_ternary_data_01_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))

    
    dataset_path =  paste("data/networks_and_trajectories_binary_data_02_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))
    
	dataset_path =  paste("data/networks_and_trajectories_ternary_data_02_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))
    
    dataset_path =  paste("data/networks_and_trajectories_binary_data_04_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))
    
	dataset_path =  paste("data/networks_and_trajectories_ternary_data_04_",i,sep="")
    load(paste(dataset_path,".RData", sep=""))
    write(toJSON(networks), paste(dataset_path, ".json", sep=""))


}
