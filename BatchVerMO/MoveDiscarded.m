function MoveDiscarded(path2cnmfe,datafolder)
currentFolder = pwd;
cd(datafolder)
mkdir ActuallyUsedInCNMFE
load(path2cnmfe)
Filename=[];
destination_Filename=[];
for i=1:length(neuron_batchMO)
    filename=neuron_batchMO(i).FileOrigin.name;
    destination_filename=['ActuallyUsedInCNMFE' '/' neuron_batchMO(i).FileOrigin.name];
    movefile(filename,destination_filename)
end
cd(currentFolder)