function [X,C,Cell,Cv,N,T,Y,xExternal,Ymid]= MeshgenGrid(Set)
%% Input 
%    - Set.Prisms- % =0; (default) Meshed with tets (Diagonalised nodal netwowrk),
                   % =1 Meshed with triangular vertical prisms
%    - Set.nx      -Set.ny   -Set.nz  -Set.h
%    - Set.XRand : % Randomise initial positions. Recomended value ~0.2
%    - Set.midY

% Meshgen builds division of paralepiped into prisms with triangle base/top
% (not tetrahedra or Delaunay)

nx=Set.nx;
ny=Set.ny;
nz=Set.nz;
Prisms=Set.Prisms;
XRand=Set.XRand;

t=0;
dim=3;
if ~exist('nz','var')
    dim=2;
    nz=0;
elseif nz==0
    dim=2;
end


X=zeros((nx+1)*(ny+1)*(nz+1),dim);
nodes=size(X,1);
for z=0:nz
    for i=0:nx
        if dim==2
            t=t+1;clear
            X(t,:)=[z,i];
        else
            for j=0:ny
                t=t+1;
                X(t,:)=[i,j,z*Set.h];
                if XRand>0 && j>0 && i>0 && j<ny && i<nx
                    X(t,1)=i+XRand*(2*rand(1)-1); % Randomise positions
                    X(t,2)=j+XRand*(2*rand(1)-1); % Randomise positions
                end
            end
        end
    end
end

if nargout ==1
    return
end 
%X(:,1:2)=X(:,1:2)+0.15*rand(size(X,1),2);

%%%%%%%%%%%%%%%%%%%%%%%%% 3D mesh generator
idx = reshape(1:(ny+1)*(nx+1)*(nz+1),[ny+1,nx+1,nz+1]);
P1 = idx(1:end-1,1:end-1,1:end);
P1=P1(:);
P2 = idx(1:end-1,2:end,1:end);
P2=P2(:);
P3 = idx(2:end,1:end-1,1:end);
P3=P3(:);
P4 = idx(2:end,2:end,1:end);
P4=P4(:);
% %
T=[P1 P2 P3 % all left bottom triangles of each square
    P2 P4 P3]; % all right top triangles of each square
% Complete bar connectivity C
idx = reshape(1:(ny+1)*(nx+1)*(nz+1),[ny+1,nx+1,nz+1]);
P1x = idx(1:end,1:end-1,1:end);
P1x=P1x(:);
P2x = idx(1:end,2:end,1:end);
P2x=P2x(:);
P1y = idx(1:end-1,1:end,1:end);
P1y=P1y(:);
P3y = idx(2:end,1:end,1:end);
P3y=P3y(:);
C=[P1x P2x % On plane vertical
    P1y P3y % On plane horizotal
    P2 P3]; % On plane diagonal
if ~Prisms && nz==1
    % Add top -Bottom bars
    P1 = idx(1:end,1:end,nz);
    P1=P1(:);
    P2 = idx(1:end,1:end,nz+1);
    P2=P2(:);
    Cl=[P1 P2]; % Top-bottom vertical
    cdd=[P1;P2];
    cdx=zeros((ny+1)*nx,2); % Top-bottom diagonal along x
    for i=1:(ny+1)*nx
        cdx(i,:)=[cdd(i,:),cdd(i+((nx+1)*(ny+1)+(ny+1)),:)];
    end
    cdy=zeros((nx+1)*ny,2); % Top-bottom diagonal along y
    ii=0;
    for i=1:(((nx+1)*(ny+1))-1)
        if i==(ny+1)+ii
            ii=ii+(nx+1);
            i=i+1;
        end
        cdy(i,:)=[cdd(i,:),cdd(cdd(i,:)+(nx+1)*(ny+1)+1,:)];
        b = cdy(any(cdy,2),:);
    end
    C=[C
        Cl % bottom-top bars
        b
        cdx];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VERTICES
%%%%%% centers of the cells
% Define shape functions for computing vertex positions
nvert=size(T,1);
N=ones(nvert,size(T,2))/3;
Y=zeros(nvert,dim);
for i=1:nvert
    Y(i,:)=1/3*sum(X(T(i,:),:),1);
end
%%%%%%exclude boundaries
VnodAll=1:nodes/2; % Internal nodes at the bottom
for i=1:nodes/2
    if (X(i,1)<eps || X(i,2)<eps  ||...
            X(i,1)==nx || X(i,2)==ny)
        VnodAll(i)=0;
    end
end
VnodAll(VnodAll==0)=[];
ncell=length(VnodAll) ;
Cell=cell(ncell,1);
%tableA =[]%zeros(28,3)
for i=1:ncell
    Vnod=VnodAll(i);
    xTriBot1=[Vnod+ny+1  , Vnod+1     , Vnod;
        Vnod+1     , Vnod-(ny)  , Vnod;
        Vnod-(ny)  , Vnod-(ny+1), Vnod;
        Vnod-(ny+1), Vnod-1     , Vnod;
        Vnod-1     , Vnod+ny    , Vnod;
        Vnod+(ny)  , Vnod+ny+1  , Vnod];

    Cell{i}.xBot=[Vnod+ny+1  , Vnod;
        Vnod+1     ,  Vnod;
        Vnod-(ny)  ,  Vnod;
        Vnod-(ny+1),  Vnod;
        Vnod-1     ,  Vnod;
        Vnod+(ny)  ,  Vnod];

    %%%% Find the external vertices around the free internal node Top
    xTriTop1=[(Vnod+ny+1)+(nx+1)*(ny+1)  , (Vnod+1)+(nx+1)*(ny+1)     ,Vnod+((nx+1)*(ny+1));
        (Vnod+1)+(nx+1)*(ny+1)     , (Vnod-(ny))+(nx+1)*(ny+1)  , Vnod+((nx+1)*(ny+1));
        (Vnod-(ny))+(nx+1)*(ny+1)  , (Vnod-(ny+1))+(nx+1)*(ny+1), Vnod+((nx+1)*(ny+1));
        (Vnod-(ny+1))+(nx+1)*(ny+1), (Vnod-1)+(nx+1)*(ny+1)     , Vnod+((nx+1)*(ny+1));
        (Vnod-1)+(nx+1)*(ny+1)     , (Vnod+ny)+(nx+1)*(ny+1)    , Vnod+((nx+1)*(ny+1));
        (Vnod+(ny))+(nx+1)*(ny+1)  , (Vnod+ny+1)+(nx+1)*(ny+1)  , Vnod+((nx+1)*(ny+1))];


    Cell{i}.xTop=[(Vnod+ny+1)+(nx+1)*(ny+1)  , Vnod+((nx+1)*(ny+1));
        (Vnod+1)+(nx+1)*(ny+1)     ,  Vnod+((nx+1)*(ny+1));
        (Vnod-(ny))+(nx+1)*(ny+1)  ,  Vnod+((nx+1)*(ny+1));
        (Vnod-(ny+1))+(nx+1)*(ny+1),  Vnod+((nx+1)*(ny+1));
        (Vnod-1)+(nx+1)*(ny+1)     ,  Vnod+((nx+1)*(ny+1));
        (Vnod+(ny))+(nx+1)*(ny+1)  ,  Vnod+((nx+1)*(ny+1))];

    % for cell{i}.vbot
    css=[];
    ti=sort(T,2);
    tt=sort(xTriBot1,2);
    k=1;
    I=ones(size(ti,1),1);
    for j=1:size(tt,1)
        for ii=1:size(ti,1)
            if tt(j,:)==ti(ii,:)
                css(k,:)=ii;
                I(ii)=0;
                k=k+1 ;
            else
            end
            Cell{i,:}.vbot=css;
        end
    end
    csc=[];
    ti=sort(T,2);
    tt=sort(xTriTop1,2);
    k=1;
    I=ones(size(ti,1),1);
    for j=1:size(tt,1)
        for ii=1:size(ti,1)
            if tt(j,:)==ti(ii,:)
                csc(k,:)=ii;
                I(ii)=0;
                k=k+1 ;
            else
            end
            Cell{i,:}.vtop=csc;
        end
    end
end
% Cv. List of vertex bars
%------------------------------------ Build vertex connectivity 
nv=0; % Total numberof vertex bars
barCell=28; % bar/Cell
Cv=zeros(ncell*barCell,3);
kk=0;                                % NUMBER OF TRIANGLES 
nMidY=0;                             % mid-plane vertices counter
if Set.ModelTop==2 || Set.ModelTop==3 
    IndicatorFaces=zeros(4*ncell,3);     %  to know which face is already triangularized         
    FaceCunt=0;                          % Lateral face counter    
end 
for i=1:ncell
    cellCv=zeros(80);   % cellular bar elements   %%%%% ADDED MALIK
    nk=0;                                         %%%%% ADDED MALIK
    % Top vertex bars
    Yt=Cell{i}.vtop;
    Yt(end+1)=Cell{i}.vtop(1);
    for e=1:length(Yt)-1
        nv=nv+1;
        Cv(nv,:)=[Yt(e) Yt(e+1) 2];
        cellCv(nk+1)=nv;                          
        nk=nk+1;                                  
        kk=kk+1;                                  
    end
    % Bottom vertex bars
    Yt=Cell{i}.vbot;
    Yt(end+1)=Cell{i}.vbot(1);
    for e=1:length(Yt)-1
        nv=nv+1;
        Cv(nv,:)=[Yt(e) Yt(e+1) 1];
        cellCv(nk+1)=nv;                          
        nk=nk+1;                                  
        kk=kk+1;                                  
    end
    % Lateral Vertex bars
    Yb=Cell{i}.vbot;
    Yb(end+1)=Cell{i}.vbot(1);
    Yt=Cell{i}.vtop;
    Yt(end+1)=Cell{i}.vtop(1);
    %  triangulate lateral surfaces
    Cell{i}.Tri=[];
    if Set.ModelTop==1
%------------------------Without mid-plane vertices------------------------
       for e=1:length(Yt)-1
            % vertical
            Cv(nv+1,:)=[Yt(e) Yb(e) 3];
            cellCv(nk+1)=nv+1;                       
            nk=nk+1;                                 
            nv=nv+1;                                 
            % diagonal
            Cv(nv+1,:)=[Yb(e) Yt(e+1) 3];
            cellCv(nk+1)=nv+1;                       
            nk=nk+1;                                 
            nv=nv+1;  
            % triangles 
            FaceTri=[Yt(e)   Yb(e)  Yt(e+1);
                     Yb(e)   Yb(e+1)  Yt(e+1)];
            Cell{i}.Tri=[Cell{i}.Tri ;FaceTri];
            kk=kk+size(FaceTri,1);
       end
%--------------------------------------------------------------------------
    else
%------------------------With mid-plane vertices---------------------------
        for e=1:length(Yt)-1
            iii=Yb(e);  jjj=Yb(e+1);
            III=Yt(e);  JJJ=Yt(e+1);
            % Check if the face already has vertex. 
            [VFace,FaceSt]=CheckFace(iii,jjj,IndicatorFaces);
            if FaceSt==1 % The face already triangularized
                FaceTri=[iii   jjj  VFace;
                         III   iii  VFace;
                         JJJ   III  VFace;
                         jjj   JJJ  VFace];
                Cv(nv+1,:)=[iii VFace 3];
                Cv(nv+2,:)=[jjj VFace 3];
                Cv(nv+3,:)=[III VFace 3];
                Cv(nv+4,:)=[JJJ VFace 3];
                Cv(nv+5,:)=[iii III 3];
                cellCv(nk+1:nk+5)=nv+1:nv+5;
                nv=nv+5;
                nk=nk+5;                            
            elseif FaceSt==0 % the face was not triangularized
                % Add Vertex in the center     
                AddY=sum(Y([iii jjj III JJJ],:))/4;  % interpolate 
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCunt=FaceCunt+1;
                IndicatorFaces(FaceCunt,:)=[iii jjj IAddY];
                % bulid bar elements 
                Cv(nv+1,:)=[iii IAddY 3];
                Cv(nv+2,:)=[jjj IAddY 3];
                Cv(nv+3,:)=[III IAddY 3];
                Cv(nv+4,:)=[JJJ IAddY 3];
                Cv(nv+5,:)=[iii III 3];
                cellCv(nk+1)=nv;
                cellCv(nk+1:nk+5)=nv+1:nv+5;
                nv=nv+5;
                nk=nk+5;
                FaceTri=[iii   jjj  IAddY;
                         III   iii  IAddY;
                         JJJ   III  IAddY;
                         jjj   JJJ  IAddY];
            end
            Cell{i}.Tri=[Cell{i}.Tri ;FaceTri];
            kk=kk+size(FaceTri,1);
        end
%--------------------------------------------------------------------------
    end 
    cellCv(nk+1:end)=[];
    Cell{i}.CCv=cellCv;
    Cell{i}.eType=Cv(cellCv,3);
end



%% Cell center nodes and total number of elements
Cv(nv+1:end,:)=[];
Cv=Cv(:,[1 2]);
for i=1:ncell
    Cell{i}.xBott=VnodAll(i);
    Cell{i}.xTopp=VnodAll(i)+((nx+1)*(ny+1));
    Cell{i}.nTri=kk;
    Cell{i}.nElem=size(Cv,1);
end
xExternal=[];
Set.nMidY=nMidY;
 T=[T;zeros(Set.nMidY,3)];
 N=[N;zeros(Set.nMidY,3)];
Ymid=zeros(size(Y,1)-nMidY,1);  
YYmid=length(Ymid)+1:size(Y,1);
Ymid=[Ymid; YYmid'];                    
end 


function [VFace,FaceSt]=CheckFace(i,j,IndicatorFaces)
         TEST=ismember(IndicatorFaces(:,1:2),[i j]); 
         RowTEST=sum(TEST,2);
         [r,~]=find(RowTEST==2); 
         if isempty(r)
             % faces is not trianglized
             VFace=[];
             FaceSt=0;
         else  % faces is already trianglized
             VFace=IndicatorFaces(r,3);
             FaceSt=1;
         end 
end

%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  2D MESH GENERATOR
% if nz==0
%     nelem=3*nx*ny+nx+ny;
%     C = zeros(nelem,2);
%     %% for vertical
%     k=0;
%     for i=1:nx*(ny+1)%6
%         for j=i+(ny+1)
%             k=k+1;
%             C(k,:) =[i,j];
%         end
%     end
%     %%%for horizontal
%     for j = 1: nx+1
%         for i = 1:ny
%             n1 = i + ( j - 1 ) * ( nx + 1 );
%             n2 = i + 1 + ( j - 1 ) * ( nx + 1 );
%             k=k+1;
%             C(k,:) = [n1,n2];
%         end
%     end
%     %% for diagonal
%     for j = 1: nx
%         for i = 1:ny
%             n1 = i + ( j - 1 ) * ( nx + 1 );
%             n3 = i + 1 + j * ( nx + 1 );
%             k=k+1;
%             C(k,:) =[n1,n3];
%         end
%     end
%     T=[]; % TO DO
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END MESH GENERATOT FOR 2D

