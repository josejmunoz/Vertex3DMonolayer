function Cell=Volume(Cell,X,Y)
ncell=size(Cell,1);
%%% top area
for i=1:ncell
    vtop=0;
    areatop=0;
    Yc=Y(Cell{i}.vtop,:);
    Yc(end+1,:)=Yc(1,:);
    x= X(Cell{i}.xTopp,:);
    for t=1: size(Cell{i}.vtop,1)
        Yi=Yc(t,:);
        Yj=Yc(t+1,:);
        a=[Yi;Yj;x];
        ab=[Yi(1:2),1;Yj(1:2),1;x(1:2),1];
        vtop=vtop+det(a)/6;
        areatop=polyarea(a(:,1),a(:,2));
    end
    Cell{i}.Vol=vtop;
    Cell{i}.areatop=areatop;
end

% Bottom
for i=1:ncell
    vbot=0;
    areabot=0;
    Yc=Y(Cell{i}.vbot,:);
    Yc(end+1,:)=Yc(1,:);
    x= X(Cell{i}.xBott,:);
    for t=size(Cell{i}.vbot,1):-1:1
        Yi=Yc(t,:);
        Yj=Yc(t+1,:);
        a=[Yj;Yi;x];
        b=[Yi(1:2),1;Yj(1:2),1;x(1:2),1];
        vbot=vbot+det(a)/6;
        areabot=polyarea(a(:,1),a(:,2));
    end
    Cell{i}.Vol=Cell{i}.Vol+vbot;
    Cell{i}.areabot=areabot;
end


%----------------------Malik ADDED ----------------------------------------
for i=1:ncell
    vlat=0;
    arealat=0;
    Tris=Cell{i}.Tri;
    for t=1:size(Tris,1)
        YTri=Y(Tris(t,:),:);
        T1=det(YTri)/6;
        Area=cross(YTri(1,:)-YTri(3,:),YTri(2,:)-YTri(3,:));
        vlat=vlat+T1;
        arealat=arealat+norm(Area)/2;
    end 
    Cell{i}.Vol=Cell{i}.Vol+vlat;
    Cell{i}.arealat=arealat;
end 
%----------------------Malik ADDED ----------------------------------------



end % end 

 