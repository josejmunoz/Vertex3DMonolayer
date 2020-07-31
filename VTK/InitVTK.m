function InitVtk()
% Creates Folder "VTKResults" and deletes old *.vtk files
% INPUT:
% REMARK:
str0='VTKResults'; % Nodal
newSubFolder = strcat(pwd,Esc,str0);
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);
delete *.*
cd '..'