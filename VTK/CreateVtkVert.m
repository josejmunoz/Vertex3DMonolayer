function CreateVtkVert(step,Y,Cell,Cv,L,L0,s,Ablated,Y0,Set,NameFile)
% Plots Vertex database
% Cv = bar connectivty of vertices
% s  = stresses of vertex bars
ncell=length(Cell);
% str0='VTKResults';
str0=NameFile;
str1='V';
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
nodes=size(Y,1);
fprintf(file,'%s %d %s\n','POINTS',nodes,'float');
for i=1:size(Y,1)
    fprintf(file,' %f %f %f\n',Y(i,1),Y(i,2),Y(i,3));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bar element
% barCell=24; % bars/cell
% nele=barCell*ncell; % Sum of bar elements
nele=Cell{1}.nElem;
if Ablated.Exist
    filew=fopen(strcat(str0,str1,'W',str3,str2),'w');
    fprintf(filew,'%s\n','# vtk DataFile Version 3.98');  %%% check the version no
    fprintf(filew,'%s\n','Vertex_vtk');
    fprintf(filew,'%s\n','ASCII');
    fprintf(filew,'%s\n','DATASET UNSTRUCTURED_GRID');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nodes=size(Y,1);
    fprintf(filew,'%s %d %s\n','POINTS',nodes,'float');
    for i=1:size(Y,1)
        fprintf(filew,' %f %f %f\n',Y(i,1),Y(i,2),Y(i,3));
    end
    nelew=sum(Ablated.Cv);
    nelenw=nele-nelew;
    fprintf(file,'%s %d %d\n','CELLS',nelenw,3*nelenw);
    fprintf(filew,'%s %d %d\n','CELLS',nelew,3*nelew);
    lnodd=Cv-1;
    for j=1:size(lnodd,1)
        if Ablated.Cv(j)>0
            fprintf(filew,'%d %d %d\n',size(lnodd(j,:),2),lnodd(j,1),lnodd(j,2));
        else
            fprintf(file,'%d %d %d\n',size(lnodd(j,:),2),lnodd(j,1),lnodd(j,2));
        end
    end
else
    nelenw=nele;
    nelew=0;
    fprintf(file,'%s %d %d\n','CELLS',nelenw,3*nelenw);
    lnodd=Cv-1;
    for j=1:size(lnodd,1)
        fprintf(file,'%d %d %d\n',size(lnodd(j,:),2),lnodd(j,1),lnodd(j,2));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(file,'%s %d\n','CELL_TYPES',nelenw);
if Ablated.Exist
    fprintf(filew,'%s %d\n','CELL_TYPES',nelew);
    for j=1:nelenw
        fprintf(file,'%d\n',3);
    end
    for j=1:nelew
        fprintf(filew,'%d\n',3);
    end
else
    for j=1:size(Cv,1)
        fprintf(file,'%d\n',3);
    end
end
%%%%%%Results%%%%%%%%%%
fprintf(file,'%s %d\n','POINT_DATA',nodes);
fprintf(file,'%s \n','VECTORS  Displacements float');
for k=1:nodes
    fprintf(file,'%f   %f   %f \n',Y(k,1:3)-Y0(k,1:3));
end

fprintf(file,'%s %d \n','CELL_DATA',nelenw);
fprintf(file,'%s \n','SCALARS stressTot float');
fprintf(file,'%s \n','LOOKUP_TABLE default');
if Ablated.Exist
    fprintf(filew,'%s %d \n','CELL_DATA',nelew);
    fprintf(filew,'%s \n','SCALARS stressTot float');
    fprintf(filew,'%s \n','LOOKUP_TABLE default');
    for i=1:size(Cv,1)
        if Ablated.Cv(i)>0
            fprintf(filew,'%f\n',s(i));
        else
            fprintf(file,'%f\n',s(i));
        end
    end
    % Contractility TODO
    fprintf(file,'%s \n','SCALARS Contractility double 1');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    fprintf(filew,'%s \n','SCALARS Contractility double 1');
    fprintf(filew,'%s \n','LOOKUP_TABLE default');
    %---------------------------------Malik (begin commnet)----------------
%     for i=1:ncell
%         ntop=length(Cell{i}.vtop); 
%         for e=1:ntop % top
%             wound=Ablated.Cell{i}.Top(e);
%             epsCW=CheckAblated(wound,Set);
%             if Ablated.Cell{i}.Ablated
%                 fprintf(filew,'%f\n',epsCW);
%             else
%                 fprintf(file,'%f\n',epsCW);
%             end
%         end
%         nbot=length(Cell{i}.vbot); 
%         for e=1:ntop % Bottom
%             wound=Ablated.Cell{i}.Bot(e);
%             epsCW=CheckAblated(wound,Set);
%             if Ablated.Cell{i}.Ablated
%                 fprintf(filew,'%f\n',epsCW);
%             else
%                 fprintf(file,'%f\n',epsCW);
%             end
%         end
%         for e=1:ntop % LAteral
%             woundV=Ablated.Cell{i}.LatV(e);
%             epsCW=CheckAblated(woundV,Set);
%             woundD=Ablated.Cell{i}.LatD(e);
%             epsCWd=CheckAblated(woundD,Set);
%             if Ablated.Cell{i}.Ablated
%                 fprintf(filew,'%f\n',epsCW);
%                 fprintf(filew,'%f\n',epsCWd);
%             else
%                 fprintf(file,'%f\n',epsCW);
%                 fprintf(file,'%f\n',epsCWd);
%             end
%         end
%     end
%     fclose(filew);
    %---------------------------------Malik (end commnet)------------------
    %---------------------------------Malik (Added begin)------------------

    for i=1:ncell
        nElem=length(Cell{i}.CCv); 
        for e=1:nElem % top
            wound=Ablated.Cell{i}.eType(e);
            epsCW=CheckAblated(wound,Set);
            if Ablated.Cell{i}.Ablated
                fprintf(filew,'%f\n',epsCW);
            else
                fprintf(file,'%f\n',epsCW);
            end
        end
    end 
    fclose(filew);
    %---------------------------------Malik (Added end) -------------------
    
    
    
    
else
    for i=1:size(Cv,1)
        fprintf(file,'%f\n',s(i));
    end
    fprintf(file,'%s \n','SCALARS RestLength double ');
    fprintf(file,'%s \n','LOOKUP_TABLE default');
    for i=1:size(Cv,1)
        fprintf(file,'%f\n',(L(i)-L0(i))/L0(i));
    end
    
    
    
    
end
fclose(file);
cd '..'