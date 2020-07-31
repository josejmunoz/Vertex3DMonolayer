function [Cv,Cell,Y,Set,Ymid]=MatrixCvSameTB(Cell,Y,Set,Ablated)
% Build vertex connectivity with the same top-bottom triangulation 
%% Output
% - Cv(i,:) = vertices connecting vertex bar element i
% - Cell  .CCv= Vector of the global IDs of the bar-elements in the matrix Cv
%         .eType= Vector of size(CCv) definse the type of the bar elements
%                = 1 bottom element , =2 top elements , =3 Lateral 
%         .Tri= matrix with triangles that defines the lateral surface
%         .nTri: Total number of surface triangles 
%         .nElem: Total number of bar-elements (size(Cv,1))
% - Ymid: vector for size(Y) 
%      - Ymid(i)=0 --> Top/bottom vertex 
%      - Ymid(i)=1 --> mid-plane vertex or vertex on the lateral side  
% - Set.nYmid: total number of mid-plane vertices
% - Set.nvert: Totla number of vertices. 
%% Initialize 
ncell=length(Cell);
Cv=zeros(ncell*24,3);
k=0;
kk=0;
IndicatorFaces=zeros(4*ncell,3);     %  to know which face is already triangularized   
FaceCunt=0;
nMidY=0;
if ~isempty(Ablated) && ~Ablated.Exist
    Ablated=[];
end 
if Set.ModelTop==2 % 
%% Addtional vertices (dofs) in mid-plane or lateral surfaces 
%                 -------
%                 |\   /|
%                 | \ / |
%                 |  *  |
%                 | / \ |
%                 |/   \|
%                 -------
    for i=1:ncell
        cellCv=zeros(80);   % cellular bar elements 
        nk=0;
        
        vtop=[Cell{i}.vtop' Cell{i}.vtop(1)];
        vbot=[Cell{i}.vbot' Cell{i}.vbot(1)];
        % Top
        for v=1:length(vtop)-1 % loop around all vertices
            Cv(k+1,:)=[vtop(v) vtop(v+1) 2]; % horiz top
            cellCv(nk+1)=k+1;                        
            nk=nk+1;                                  
            k=k+1;
            kk=kk+1;
        end
        % Bottom 
        for v=1:length(vbot)-1 % loop around all vertices
            Cv(k+1,:)=[vbot(v) vbot(v+1) 1]; % horiz bottom
            cellCv(nk+1)=k+1;                        
            nk=nk+1;                                  
            k=k+1;
            kk=kk+1;
        end
        % Lateral
        Cell{i}.Tri=[];
        
        for v=1:length(vbot)-1 % loop around all vertices
            if ~isempty(Ablated) && ismember(vbot(v),Ablated.VRingBot(:,1)) && ismember(vbot(v+1),Ablated.VRingBot(:,1)) && Set.YmidWound==0 
                Tri=[vbot(v)   vbot(v+1)  vtop(v+1);
                    vtop(v+1)   vtop(v)  vbot(v)];
                Cv(k+1,:)=[vbot(v) vtop(v) 3];
                Cv(k+2,:)=[vbot(v) vtop(v+1) 3];
                cellCv(nk+1:nk+2)=k+1:k+2;
                k=k+2;
                nk=nk+2; 
            else 
                % Check if the face already has vertex. 
                [VFace,FaceSt]=CheckFace(vbot(v),vbot(v+1),IndicatorFaces);
                if FaceSt==1 % The face already triangularized
                    Tri=[vbot(v)   vbot(v+1)  VFace;
                             vtop(v)   vbot(v)  VFace;
                             vtop(v+1)   vtop(v)  VFace;
                             vbot(v+1)   vtop(v+1)  VFace];
                    Cv(k+1,:)=[vbot(v) VFace 3];
                    Cv(k+2,:)=[vbot(v+1) VFace 3];
                    Cv(k+3,:)=[vtop(v) VFace 3];
                    Cv(k+4,:)=[vtop(v+1) VFace 3];
                    Cv(k+5,:)=[vbot(v) vtop(v) 3];
                    cellCv(nk+1:nk+5)=k+1:k+5;
                    k=k+5;
                    nk=nk+5;                         

                elseif FaceSt==0 % the face was not triangularized
                    % Add Vertex in the center     
                    AddY=sum(Y([vbot(v) vbot(v+1) vtop(v) vtop(v+1)],:))/4;  % interpolate 
                    IAddY=1+size(Y,1);
                    Y=[Y;AddY];
                    nMidY=nMidY+1;
                    FaceCunt=FaceCunt+1;
                    IndicatorFaces(FaceCunt,:)=[vbot(v) vbot(v+1) IAddY];
                    % bulid bar elements 
                    Cv(k+1,:)=[vbot(v) IAddY 3];
                    Cv(k+2,:)=[vbot(v+1) IAddY 3];
                    Cv(k+3,:)=[vtop(v) IAddY 3];
                    Cv(k+4,:)=[vtop(v+1) IAddY 3];
                    Cv(k+5,:)=[vbot(v) vtop(v) 3];
                    cellCv(nk+1)=k;
                    cellCv(nk+1:nk+5)=k+1:k+5;
                    k=k+5;
                    nk=nk+5;
                    Tri=[vbot(v)   vbot(v+1)  IAddY;
                             vtop(v)   vbot(v)  IAddY;
                             vtop(v+1)   vtop(v)  IAddY;
                             vbot(v+1)   vtop(v+1)  IAddY];
                end
            end 
            Cell{i}.Tri=[Cell{i}.Tri ;Tri];
            kk=kk+size(Tri,1);
        end 
        cellCv(nk+1:end)=[];
        Cell{i}.CCv=cellCv;
        Cell{i}.eType=Cv(cellCv,3);
    end
elseif   Set.ModelTop==1 
%%   triangulate without mid-plane vertices     
%                 -------
%                 |\    |
%                 | \   |
%                 |  \  |
%                 |   \ |
%                 |    \|
%                 -------
    for i=1:ncell
        cellCv=zeros(80);   % cellular bar elements 
        nk=0;
        
        vtop=[Cell{i}.vtop' Cell{i}.vtop(1)];
        vbot=[Cell{i}.vbot' Cell{i}.vbot(1)];
        % Top
        for v=1:length(vtop)-1 % loop around all vertices
            Cv(k+1,:)=[vtop(v) vtop(v+1) 2]; % horiz top
            cellCv(nk+1)=k+1;                        
            nk=nk+1;                                  
            k=k+1;
            kk=kk+1;
        end
        % Bottom 
        for v=1:length(vbot)-1 % loop around all vertices
            Cv(k+1,:)=[vbot(v) vbot(v+1) 1]; % horiz bottom
            cellCv(nk+1)=k+1;                        
            nk=nk+1;                                  
            k=k+1;
            kk=kk+1;
        end
        % Lateral
        Cell{i}.Tri=[];
        for v=1:length(vbot)-1 % loop around all vertices
            Cv(k+1,:)=[vtop(v) vbot(v) 3]; % vertical
            Cv(k+2,:)=[vbot(v) vtop(v+1) 3]; % diagonal
            cellCv([nk+1 nk+2])=[k+1 k+2];                        
            nk=nk+2;
            k=k+2;
            Tri=[vtop(v) vbot(v)   vtop(v+1);
                 vbot(v) vbot(v+1) vtop(v+1)];
           Cell{i}.Tri=[Cell{i}.Tri; Tri];
           kk=kk+size(Tri,1);
        end
        cellCv(nk+1:end)=[];
        Cell{i}.CCv=cellCv;
        Cell{i}.eType=Cv(cellCv,3);
    end
end 


for n=1:ncell                             
Cell{n}.nTri=kk;
Cell{n}.nElem=k;
end 

Set.nMidY=nMidY;           
Set.nvert=size(Y,1);
Ymid=zeros(size(Y,1)-nMidY,1);  
YYmid=length(Ymid)+1:size(Y,1);
Ymid=[Ymid; YYmid'];                     

Cv(k+1:end,:)=[];
% eType=Cv(:,3);
Cv=Cv(:,[1 2]);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
