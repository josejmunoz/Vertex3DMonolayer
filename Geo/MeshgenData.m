function [X,C,Cell,Cv,N,T,Y,xExternal,Ymid]= MeshgenData(Set,X)
%% Input 
%    - Set.DiffTB
%    - Set.RemodelTolF
%    - Set.RemodelDelta

% [X,C,Cell,Cv,L,N,T,Y,xExternal]=Remodel(Set,X);

if isstruct(X) % Different top-bottom
    Xb=(1:size(X.Xb,1))';
    Xt=(size(X.Xb,1)+1:size(X.Xb,1)+size(X.Xt,1))';
    X=[X.Xb;
       X.Xt];
    nnodes=length(Xt);
else          % same top-bottom
    nnodes=size(X,1);
    X=X(:,1:2);
    Xb=(1:nnodes)';
    Xt=(nnodes+1:2*nnodes)';
    X=[X zeros(nnodes,1)
       X Set.h*ones(nnodes,1)];
end
Set.nodes=2*nnodes;
Set.dim=3;
CellOld=[];
% [Set,Tb,Tt,X,X0,Xb,Xt]=XbXt(Cell,Set,T,X,X0);
Xbt=zeros(size(Xt));
Xbt(Xb)=Xt; % Transformation from top bottom to top

if  Set.ModelTop ==3 
    %% Different Top/Bottom
    [Cb,Edges]=delaunayF([],X(Xb,1),X(Xb,2),Set.RemodelTolF,Set.RemodelDelta); % Delaunay with filtering on boundary
    [Ct,~]    =delaunayF([],X(Xt,1),X(Xt,2),Set.RemodelTolF,Set.RemodelDelta); % Delaunay with filtering on boundary
    Cb=[Xb(Cb(:,1)) Xb(Cb(:,2)) Xb(Cb(:,3))]; % Global numbers
    Ct=[Xt(Ct(:,1)) Xt(Ct(:,2)) Xt(Ct(:,3))]; % Global numbers
    % Search triangles per node
    [nTriNodeB,TriNodeB,nTriNodeT,TriNodeT]=ConnectivityTri(Cb,Set.nodes,Ct); 
    % Find internal/external nodes
    [xInternalB,xExternal]=ExternalInternal(Cb,Edges);
    % xExternalT=xExternalB+Set.nodes/2;                               
    xInternalT=xInternalB+Set.nodes/2;
    % Vertex positions and updated x
    [nTrianglesB,Y,nTrianglesT]=GetY(Cb,Ct,X);
    Set.nvert=nTrianglesB+nTrianglesT;
    [Cell]=xBottvbot(Cb,CellOld,xInternalB,nTriNodeB,TriNodeB,Y,Ct,xInternalT,nTriNodeT,TriNodeT);
     T=[Cb; Ct];
     N=ones(nTrianglesB+nTrianglesT,3)/3;
     [Cv,Cell,Y,Set,Ymid]=MatrixCvDiffTB(Cell,Cb,Ct,Y,Set,[]);
     C=MatrixC(Cb,nTriNodeB,TriNodeB,Xb,Xbt,Xt,Ct,nTriNodeT,TriNodeT);
     T=[T;zeros(Set.nMidY,3)];
     N=[N;zeros(Set.nMidY,3)];   
else 
    %% Same Top/Bottom
    [Cb,Edges]=delaunayF([],X(Xb,1),X(Xb,2),Set.RemodelTolF,Set.RemodelDelta); % Delaunay with filtering on boundary
    Ct=[Xt(Cb(:,1)) Xt(Cb(:,2)) Xt(Cb(:,3))];
    % Search triangles per node
    [nTriNode,TriNode]=ConnectivityTri(Cb,Set.nodes);
    % Find internal/external nodes
    [xInternal,xExternal]=ExternalInternal(Cb,Edges);
    % Vertex positions and updated x
    [nTriangles,Y]=GetY(Cb,Ct,X);
    Set.nvert=2*nTriangles;  
    [Cell]=xBottvbot(Cb,CellOld,xInternal,nTriNode,TriNode,Y);
    T=[Cb; Ct];
    N=ones(2*nTriangles,3)/3;
    C=MatrixC(Cb,nTriNode,TriNode,Xb,Xbt,Xt);
    [Cv,Cell,Y,Set,Ymid]=MatrixCvSameTB(Cell,Y,Set,[]);
         T=[T;zeros(Set.nMidY,3)];
    N=[N;zeros(Set.nMidY,3)];
end 
if min(size(Cb))==0
    error('Could not generate cells with given data. Check parameters and constants.')
end
% C=MatrixC(Cb,Ct,nTriNodeB,nTriNodeT,TriNodeB,TriNodeT,Xb,Xbt,Xt,Set.DiffTB);
end





