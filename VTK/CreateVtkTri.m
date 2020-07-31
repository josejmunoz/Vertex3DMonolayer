function CreateVtkTri(step,Y,Cell,Ablated,Y0,X,X0,NameFile)
% Cv = bar connectivty of vertices
% s  = stresses of vertex bars
ncell=length(Cell);
% str0='VTKResults';
str0=NameFile;
str1='T';
str2='.vtk';
str3 = num2str(step);
str13=strcat(str0,str1,str3);
newSubFolder = strcat(pwd,Esc,str0);
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);
nameout=strcat(str13,str2);
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');  %%% check the version no
fprintf(file,'%s\n','Vertex_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nvert=size(Y,1);
nodes=nvert+2*ncell;
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
for i=1:size(Y,1)
    fprintf(file,' %f %f %f\n',Y(i,1),Y(i,2),Y(i,3));
end
for i=1:ncell
    fprintf(file,' %f %f %f\n',X(Cell{i}.xTopp,:));
    fprintf(file,' %f %f %f\n',X(Cell{i}.xBott,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bar element
nele=Cell{1}.nTri;
%triCell=24;
% nele=triCell*ncell; 
celsize=4*nele;
nelenw=nele;

fprintf(file,'%s %d %d\n','CELLS',nelenw,celsize);
for iCell=1:ncell
    %top polygon
    nVertex=length(Cell{iCell}.vtop);
    lnodd=[Cell{iCell}.vtop(:,1)' Cell{iCell}.vtop(1,1)]-1;
    for iV=1:nVertex
        fprintf(file,'%d %d %d %d\n',3,lnodd(iV:iV+1),nvert+2*(iCell-1));
    end
    %bottom polygon
    nVertex=length(Cell{iCell}.vbot);
    lnodd=[Cell{iCell}.vbot(:,1)' Cell{iCell}.vbot(1,1)]-1;
    for iV=1:nVertex
        fprintf(file,'%d %d %d %d\n',3,lnodd(iV:iV+1),nvert+2*(iCell-1)+1);
    end
    %print end of line
    % lateral triangles
    %-----------------Malik end comment -------------------------------
    %---------------------------------Malik Added (begin)--------------

    for i=1:size(Cell{iCell}.Tri,1)
        fprintf(file,'%d %d %d %d\n',3,(Cell{iCell}.Tri(i,1)-1)...
                                      ,(Cell{iCell}.Tri(i,2)-1)...
                                      ,(Cell{iCell}.Tri(i,3)-1));
    end 
    %---------------------------------Malik Added (end) ---------------

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'%s %d\n','CELL_TYPES',nelenw);
for i=1:ncell
    nTriT=length(Cell{i}.vtop);                             %%% Malik Added
    nTriB=length(Cell{i}.vbot);                             %%% Malik Added
    nTriL=size(Cell{i}.Tri,1);                              %%% Malik Added
    for j=1:(nTriT+nTriB+nTriL) % apical + basal + lateral
        fprintf(file,'%d\n',5);
    end
end
%%%%%%Results%%%%%%%%%%
%------------------------ Malik begin comment (temporary)------------------
% fprintf(file,'%s %d\n','POINT_DATA',nodes);
% fprintf(file,'%s \n','VECTORS  Displacements float');
% for k=1:nvert
%     fprintf(file,'%f   %f   %f \n',Y(k,:)-Y0(k,:));
% end
% for k=1:ncell
%     iTop=Cell{k}.xTopp;
%     iBot=Cell{k}.xBott;
%     fprintf(file,'%f   %f   %f \n',X(iTop,:)-X0(iTop,:));
%     fprintf(file,'%f   %f   %f \n',X(iBot,:)-X0(iBot,:));
% end
%------------------------ Malik end comment (temporary)--------------------


% ADD RELATIVE VOLUME CHANGE
fprintf(file,'%s %d \n','CELL_DATA',nelenw);
fprintf(file,'%s \n','SCALARS RelVolChange double');
fprintf(file,'%s \n','LOOKUP_TABLE default');
%
if isfield(Cell{1},'Vol')
    for i=1:ncell
        ntri=ones(size(Cell{i}.Tri,1)+length(Cell{i}.vtop)+length(Cell{i}.vbot),1);      %%% Malik Added
        fprintf(file,'%f\n', (Cell{i}.Vol-Cell{i}.Vol0)/Cell{i}.Vol0*ntri);
    end
else
    for i=1:ncell
        ntri=ones(size(Cell{i}.Tri,1)+length(Cell{i}.vtop)+length(Cell{i}.vbot),1);     %%% Malik Added
        fprintf(file,'%f\n', 0*ntri);
    end
end
fclose(file);
cd '..'