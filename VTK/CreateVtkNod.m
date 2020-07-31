function CreateVtkNod(step,X0,X,lnod,Ablated,L,L0,s,Set,NameFile)
% Prints output for owunded and unwounded cells
% INPUT:
% step = step number
% X0   = initial nodal coordinates
% X    = current nodal coordinates
% lnod = nodal network connectivity
% s    = stress values
% L    = resting lengths
% REMARK:
% Same number of points in wounded and unwounded files
% Different number of elements in wounded and unwounded files
% str0='VTKResults';
str0=NameFile;
str2='.vtk';
str3 = num2str(step);
newSubFolder = strcat(pwd,Esc,str0);
if ~exist(newSubFolder, 'dir')
    mkdir(newSubFolder);
end
cd(newSubFolder);
% Write non-ablated rod elements
nameout=strcat(str0,'N',str3,str2);
file=fopen(nameout,'w');
fprintf(file,'%s\n','# vtk DataFile Version 3.98');
fprintf(file,'%s\n','Delaunay_vtk');
fprintf(file,'%s\n','ASCII');
fprintf(file,'%s\n','DATASET UNSTRUCTURED_GRID');
nele=size(lnod,1);
nodes=size(X0,1);
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
nelew=0;
if Ablated.Exist
    nameoutw=strcat(str0,'NW',num2str(step),str2);
    filew=fopen(nameoutw,'w');
    fprintf(filew,'%s\n','# vtk DataFile Version 3.98');
    fprintf(filew,'%s\n','Delaunay_vtk');
    fprintf(filew,'%s\n','ASCII');
    fprintf(filew,'%s\n','DATASET UNSTRUCTURED_GRID');
    fprintf(filew,'%s %d %s\n','POINTS',nodes,'float');
    % Count number of ablated elements
    for e=1:nele
        if Ablated.Nodal(e)>0
            nelew=nelew+1;
        end
    end
end
nelenw=nele-nelew;
%%%%%%%%%%%%
dim=size(X,2);
for i=1:nodes
    if dim==2
        fprintf(file,' %f %f %f\n',X(i,1),X(i,2),0);
    else
        fprintf(file,' %f %f %f\n',X(i,1),X(i,2),X(i,3));
    end
end
if Ablated.Exist
    for i=1:nodes
        if dim==2
            fprintf(filew,' %f %f %f\n',X(i,1),X(i,2),0);
        else
            fprintf(filew,' %f %f %f\n',X(i,1),X(i,2),X(i,3));
        end
    end
end
if Ablated.Exist
    fprintf(file,'%s %d %d\n','CELLS',nelenw,nelenw*(size(lnod(1,:),2)+1));
    fprintf(filew,'%s %d %d\n','CELLS',nelew,nelew*(size(lnod(1,:),2)+1));
    lnodd=lnod-1;
    for j=1:size(lnodd,1)
        if Ablated.Nodal(j)>0
            fprintf(filew,'%d %d %d\n',size(lnodd(j,:),2),lnodd(j,1),lnodd(j,2));
        else
            fprintf(file,'%d %d %d\n',size(lnodd(j,:),2),lnodd(j,1),lnodd(j,2));
        end
    end
    fprintf(file,'%s %d\n','CELL_TYPES',nelenw);
    for j=1:nelenw
        fprintf(file,'%d\n',3);
    end
    fprintf(filew,'%s %d\n','CELL_TYPES',nelew);
    for j=1:nelew
        fprintf(filew,'%d\n',3);
    end
else
    fprintf(file,'%s %d %d\n','CELLS',nelenw,nelenw*(size(lnod(1,:),2)+1));
    lnodd=lnod-1;
    for j=1:size(lnodd,1)
        fprintf(file,'%d %d %d\n',size(lnodd(j,:),2),lnodd(j,1),lnodd(j,2));
    end
    fprintf(file,'%s %d\n','CELL_TYPES',nelenw);
    for j=1:size(lnod,1)
        fprintf(file,'%d\n',3);
    end
end
%%%%%%%Results%%%%%%%%%%
fprintf(file,'%s %d\n','POINT_DATA',nodes);
fprintf(file,'%s \n','VECTORS  Displacements float');
for k=1:nodes
    fprintf(file,'%f   %f   %f \n',X(k,1)-X0(k,1),X(k,2)-X0(k,2),X(k,3)-X0(k,3));
end
%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'%s %d \n','CELL_DATA',nelenw);
fprintf(file,'%s \n','SCALARS stressTot double ');
fprintf(file,'%s \n','LOOKUP_TABLE default');
if Ablated.Exist
    % Also point data and add cell header
    fprintf(filew,'%s %d\n','POINT_DATA',nodes);
    fprintf(filew,'%s \n','VECTORS  Displacements float');
    for k=1:nodes
        fprintf(filew,'%f   %f   %f \n',X(k,1)-X0(k,1),X(k,2)-X0(k,2),X(k,3)-X0(k,3));
    end
    fprintf(filew,'%s %d \n','CELL_DATA',nelew);
    fprintf(filew,'%s \n','SCALARS stressTot double ');
    fprintf(filew,'%s \n','LOOKUP_TABLE default');
    for i=1:size(lnod,1)
        if Ablated.Nodal(i)
            fprintf(filew,'%f\n',s(i));
        else
            fprintf(file,'%f\n',s(i));
        end
    end
    fprintf(file,'%s \n','SCALARS Contractility double 1');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    fprintf(filew,'%s \n','SCALARS Contractility double 1');
    fprintf(filew,'%s \n','LOOKUP_TABLE default');
    for i=1:size(lnod,1)
        wound=Ablated.Nodal(i);
        if Ablated.Nodal(i)
            fprintf(filew,'%f\n',CheckAblated(wound,Set));
        else
            fprintf(file,'%f\n',CheckAblated(wound,Set));
        end
    end
    fclose(filew);
else
    for i=1:size(lnod,1)
        fprintf(file,'%f\n',s(i));
    end
    fprintf(file,'%s \n','SCALARS RestLength double ');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    for i=1:size(lnod,1)
        fprintf(file,'%f\n',(L(i)-L0(i))/L0(i));
    end
end
fclose(file);
cd '..'