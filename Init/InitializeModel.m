function [X,Xn,X0,C,Cv,T,Cell,Y,Yn,Y0,L,Ln,L0,N,Set,xExternal,Ablated,Ymid]= InitializeModel(Set)
% function that initializes tissue structure   

%% Input
%     - Set.CellCentres

% INPUT:
% nx,ny,nz=number of divisions on x,y,z
% Prisms  % =0; (default) Meshed with tets (Diagonalised nodal netwowrk),
%           =1 Meshed with triangular vertical prisms

% OUTPUT:
% X=nodal coordinates
% C=connectivty of bars
% Cv=connectivty of vertex bars
% Y  = vertex coordinates
% T   = connectivty of triangles on top and bottom surface
% Cell{i} Database of each cell i:
%         vbot= list of vertices at bottom surface 
%         vtop= list of vertices at top surface 
%         xBott=Node at center of bottom surface 
%         xTopp=Node at center of top surface 
%         xBot=Nodal elements at bottom surface
%         xTop=Nodal elements at top surface 
% cq1=Positions of vertices at bottom
% cq2=Positions of vertices at top
% c1  =Connection top nodes-vertex for each cell
% c12 =Connection bottom nodes-vertex for each cell

% Define coordinates
if ~isempty(Set.CellCentres)
    load(Set.CellCentres);
    X=X*Set.umPerPixel;
    [X,C,Cell,Cv,N,T,Y,xExternal,Ymid]= MeshgenData(Set,X);
    if (isfield(Set,'nx') || isfield(Set,'nx')) && isfield(Set,'CellCentres')
        if Set.nx>0 || Set.ny>0
            warning('Set.nx>0 or Set.ny>0, but file in Set.CellCentres is used. Set.nx and Set.ny ignored. ');
            Set.nx=0;
            Set.ny=0;
        end
    end
else
    [X]= MeshgenGrid(Set);
     X=X(1:size(X,1)/2,:);
     Xb=X;
     Xt=X; Xt(:,3)=ones(size(X,1),1)*Set.h;
%      Xt=X+rand(size(X))/50; Xt(:,3)=ones(size(X,1),1)*Set.h;
     XX.Xb=Xb;
     XX.Xt=Xt;
     [X,C,Cell,Cv,N,T,Y,xExternal,Ymid]= MeshgenData(Set,XX);

end

Cell=Volume(Cell,X,Y);
for i=1:length(Cell)
        Cell{i}.Vol0=Cell{i}.Vol;
end 
% General dimensions
Set.nodes=size(X,1);
Set.dim=size(X,2);
Set.nvert=size(Y,1);
Set.nMidY=length(Ymid(Ymid>0));

% Define rest length
L.V.L=zeros(size(Cv,1),1);
L.D.L=zeros(size(C,1),1);
for e=1:size(C,1)
    L.D.L(e)=norm(X(C(e,1),:)-X(C(e,2),:));
end
for e=1:size(Cv,1)
    L.V.L(e)=norm(Y(Cv(e,1),:)-Y(Cv(e,2),:));
end

Ln=L; 
L0=L; 
X0=X;
Y0=Y;% update vertices
Xn=X;
Yn=Y;

[Ablated,Set]=InitializeAblation(Cell,C,Cv,Set,Y);

end

%%    duplicate tissue  
%     xx1=X(:,1);
%     xx2=xx1+max(xx1);
%     yy1=X(:,2);
%     yy2=yy1+max(yy1);
%     X=[X;
%        [xx2 yy1];
%        [xx1 yy2];
%        [xx2 yy2]];
% [X,C,Cell,Cv,L,N,T,Y,xExternal]= MeshgenData(Set,X)
%     [X,C,Cell,Cv,L,N,T,Y,xExternal]=Remodel(Set,X); % Remodel if Set.Remodel>0

