function [g,K,Cell]=gKVol(Cell,X,Y,Yr,Set,T,N,Ymid)
% Computes JAcobian of volume conservation
Cell=Volume(Cell,X,Y);
% K(i,j)= derivative of g(i) wrt to x(j)
ncell=size(Cell,1);
% Residual
dim=size(X,2);
dofx=Set.nodes*Set.dim;
if Set.yRelaxation || Set.midY                              %%% Added Malik 
    dimg=(Set.nvert+Set.nodes)*dim;
%     Rv=find(Yr>0);
%     Tr=zeros(length(Rv),3);
%     T=[T;Tr];
else
    dimg=Set.nodes*dim;
end
% T=[T;zeros(Set.nMidY,3)];                                   %%% Added Malik 
g=zeros(dimg,1); % Local cell residual
if Set.Sparse
    sk=0;
    si=zeros(floor((ncell*12*6*3)^2/100),1); % Each vertex affecting 6 nodes, but reduced by an empirical factor of 100
    sj=si;
    sv=si;
    K=sparse(zeros(dimg)); % Also used in sparse
else
    K=zeros(dimg); % Also used in sparse
end
%     % Analytical residual g and Jacobian K
for i=1:ncell
    lambdaV=Set.lambdaV;
%     if min(abs(Set.Ablation(:)-i))==0
%         lambdaV=Set.lambdaV*Set.ReductionlambdaV;
%     end
    fact=lambdaV*(Cell{i}.Vol-Cell{i}.Vol0)/Cell{i}.Vol0^2;
    ge=zeros(dimg,1); % Local cell residual
    % Top triangles
    Yt=Y(Cell{i}.vtop,:);
    Yt(end+1,:)=Yt(1,:);
    Yb=Y(Cell{i}.vbot,:);
    Yb(end+1,:)=Yb(1,:);
    xb= X(Cell{i}.xBott,:);
    xt= X(Cell{i}.xTopp,:);
    ft=Cell{i}.vtop; % Top Vertex numbers
    ft(end+1,:)=ft(1,:);
    fb=Cell{i}.vbot; % bot Vertex numbers
    fb(end+1,:)=fb(1,:);
    %top triangles
    for j=1:size(Cell{i}.vtop,1)
        Y1=Yt(j,:);
        Y2=Yt(j+1,:);
        Y3=xt;
        nY=[ft(j) % vertex number Y1
            ft(j+1)
            0];
        nx=[T(ft(j),:) % Nodes associated to vertex 1
            T(ft(j+1),:)
            Cell{i}.xTopp 0 0];
        [gs,Ks]=gKDet(Y1,Y2,Y3);
        gt=gTriangleVol(nY,gs,N,Yr,Ymid); % Der of factor
        ge=AssemblegTriangleVol(ge,gt,nx,nY,Yr,Set.nodes,Ymid);
        if nargout==3
            Kt=KTriangleVol(nY,Ks/6,N,Yr,Ymid); % der(Vol)=der(det)/6
            if Set.Sparse
                [si,sj,sv,sk]=AssembleKTriangleVolSparse(Kt*fact,nx,nY,Yr,Set.nodes,si,sj,sv,sk,Ymid);
            else
                K=AssembleKTriangleVol(K,Kt*fact,nx,nY,Yr,Set.nodes,Ymid);
            end
        end
    end
    %bottom triangles (opposite direction than top, since external normal)
    for j=1:size(Cell{i}.vbot,1)
        Y1=Yb(j+1,:);
        Y2=Yb(j,:);
        Y3=xb;
        nY=[fb(j+1) % vertex number Y1
            fb(j)
            0];
        nx=[T(fb(j+1),:) % Nodes associated to vertex 1
            T(fb(j),:)
            Cell{i}.xBott 0 0];
        [gs,Ks]=gKDet(Y1,Y2,Y3);
        gt=gTriangleVol(nY,gs,N,Yr,Ymid); % Der of factor
        ge=AssemblegTriangleVol(ge,gt,nx,nY,Yr,Set.nodes,Ymid);
        if nargout==3
            Kt=KTriangleVol(nY,Ks/6,N,Yr,Ymid); % Der of gi
            if Set.Sparse
                [si,sj,sv,sk]=AssembleKTriangleVolSparse(Kt*fact,nx,nY,Yr,Set.nodes,si,sj,sv,sk,Ymid);
            else
                K=AssembleKTriangleVol(K,Kt*fact,nx,nY,Yr,Set.nodes,Ymid);
            end
        end
    end
    
%Lateral traingles
%----------------------- Malik begin comment --------------------------------
%     for j=1:size(Cell{i}.vtop,1)
%         % Triangle T1
%         Y1=Yt(j,:);
%         Y2=Yb(j,:);
%         Y3=Yt(j+1,:);
%         nY=[ft(j) % vertex number Y1
%             fb(j)
%             ft(j+1)];
%         [gs,Ks]=gKDet(Y1,Y2,Y3);
%         gt=gTriangleVol(nY,gs,N,Yr); % Der of factor
%         ge=AssemblegTriangleVol(ge,gt,T(nY,:),nY,Yr,Set.nodes);
%         if nargout==2
%             Kt=KTriangleVol(nY,Ks/6,N,Yr); % Der of gi
%             if Set.Sparse
%                 [si,sj,sv,sk]=AssembleKTriangleVolSparse(Kt*fact,T(nY,:),nY,Yr,Set.nodes,si,sj,sv,sk);
%             else
%                 K=AssembleKTriangleVol(K,Kt*fact,T(nY,:),nY,Yr,Set.nodes);
%             end
%         end
%         % Triangle T2
%         Y1=Yb(j,:);
%         Y2=Yb(j+1,:);
%         Y3=Yt(j+1,:);
%         nY=[fb(j) % vertex number Y1
%             fb(j+1)
%             ft(j+1)];
%         [gs,Ks]=gKDet(Y1,Y2,Y3);
%         gt=gTriangleVol(nY,gs,N,Yr); % Der of factor
%         ge=AssemblegTriangleVol(ge,gt,T(nY,:),nY,Yr,Set.nodes);
%         if nargout==2
%             Kt=KTriangleVol(nY,Ks/6,N,Yr); % Der of gi
%             if Set.Sparse
%                 [si,sj,sv,sk]=AssembleKTriangleVolSparse(Kt*fact,T(nY,:),nY,Yr,Set.nodes,si,sj,sv,sk);
%             else
%                 K=AssembleKTriangleVol(K,Kt*fact,T(nY,:),nY,Yr,Set.nodes);
%             end
%         end
%     end
    %-----------------------Malik end comment -----------------------------
        Tris=Cell{i}.Tri;
    for t=1:size(Tris,1)
        nY=Tris(t,:);
        Y1=Y(nY(1),:);
        Y2=Y(nY(2),:);
        Y3=Y(nY(3),:);
        [gs,Ks]=gKDet(Y1,Y2,Y3);
        gt=gTriangleVol(nY,gs,N,Yr,Ymid); % Der of factor
        ge=AssemblegTriangleVol(ge,gt,T(nY,:),nY,Yr,Set.nodes,Ymid);
        if nargout==3
            Kt=KTriangleVol(nY,Ks/6,N,Yr,Ymid); % Der of gi
            if Set.Sparse
                [si,sj,sv,sk]=AssembleKTriangleVolSparse(Kt*fact,T(nY,:),nY,Yr,Set.nodes,si,sj,sv,sk,Ymid);
            else
                K=AssembleKTriangleVol(K,Kt*fact,T(nY,:),nY,Yr,Set.nodes,Ymid);
            end
        end
    end 
 
    g=g+ge*fact/6; % Volume contribution of each triangle is det(Y1,Y2,Y3)/6
    if nargout==3
        if Set.Sparse
            K=K+lambdaV*sparse((ge)*(ge'))/6/6/Cell{i}.Vol0^2;
        else
            K=K+lambdaV*(ge)*(ge')/6/6/Cell{i}.Vol0^2;
        end
    end
end
if Set.Sparse
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg)+K;
end
end
%%
function Ymat=Cross(y)
Ymat=[0 -y(3) y(2)
    y(3) 0 -y(1)
    -y(2) y(1) 0];
end
%%
function [gs,Ks]=gKDet(Y1,Y2,Y3)
% Returns residual and  Jacobian of det(Y)=y1'*cross(y2,y3)
% gs=[der_y1 det(Y) der_y2 det(Y) der_y3 det(Y)]
% Ks=[der_y1y1 det(Y) der_y1y2 det(Y) der_y1y3 det(Y)
%     der_y2y1 det(Y) der_y2y2 det(Y) der_y2y3 det(Y)
%     der_y3y1 det(Y) der_y3y2 det(Y) der_y3y3 det(Y)]
dim=length(Y1);
gs=[cross(Y2,Y3) % der_Y1 (det(Y1,Y2,Y3)) 
    cross(Y3,Y1)
    cross(Y1,Y2)];
Ks=[ zeros(dim) -Cross(Y3)   Cross(Y2) % g associated to der wrt vertex 1
    Cross(Y3)   zeros(dim) -Cross(Y1)
    -Cross(Y2)   Cross(Y1)  zeros(dim)];
end
%%
function   ge=AssemblegTriangleVol(ge,gt,nx,nY,Yr,nnodes,Ymid)
% Assembles volume residual of a triangle of vertices (27 components)
dim=3;
nodes=size(nx,2); % 3 nodes/triangle
% correct nx according to relaxed vertices
for I=1:3 % loop on 3 vertices of triangle
    if nY(I)>0
        if Yr(nY(I))>0 || Ymid(nY(I))>0                   %%% Added Malik   % relaxed vertex  or mid-plane 
            nx(I,:)=0;
            nx(I,1)=nY(I)+nnodes;
        end
    end
    for n=1:nodes % Loop on nodes
        if nx(I,n)>0 || n==1
            idofg=(nx(I,n)-1)*dim+1:nx(I,n)*dim; % global dof
            idofl=(I-1)*dim*nodes+(n-1)*dim+1:(I-1)*dim*nodes+n*dim;
            ge(idofg)=ge(idofg)+gt(idofl);
        end
    end
end
end
%%
function Ke= AssembleKTriangleVol(Ke,Kt,n,nY,Yr,nnodes,Ymid)
% Assembles volume Jacobian of a triangle of vertices (27x27 components)
dim=3;
nodes=size(n,2); % 3 nodes/triangle
for I=1:3 % loop on 3 vertices of triangle
    if nY(I)>0
        if Yr(nY(I))>0 || Ymid(nY(I)) >1                   %%% Added Malik   % relaxed vertex or mid-plane
            n(I,:)=0;
            n(I,1)=nY(I)+nnodes;
        end
    end
end
for I=1:3 % loop on 3 vertices of triangle
    for i=1:nodes % Loop on nodes
        for J=1:3 % loop on 3 vertices of triangle
            for j=1:nodes % Loop on nodes
                if (n(I,i)>0 || i==1) && (n(J,j)>0 || j==1) && (I~=J)
                    idofg=((n(I,i)-1)*dim+1):n(I,i)*dim; % global dof
                    idofl=((I-1)*dim*nodes+(i-1)*dim)+1:(I-1)*dim*nodes+i*dim;
                    jdofg=((n(J,j)-1)*dim+1):n(J,j)*dim; % global dof
                    jdofl=((J-1)*dim*nodes+(j-1)*dim)+1:(J-1)*dim*nodes+j*dim;
                    Ke(idofg,jdofg)=Ke(idofg,jdofg)+Kt(idofl,jdofl);
                end
            end
        end
    end
end
end
%%%
function [si,sj,sv,sk]= AssembleKTriangleVolSparse(Kt,n,nY,Yr,nnodes,si,sj,sv,sk,Ymid)
% Assembles volume Jacobian of a triangle of vertices (27 components)
dim=3;
nodes=size(n,2); % 3 nodes/triangle
for I=1:3 % loop on 3 vertices of triangle
    if nY(I)>0
        if Yr(nY(I))>0 || Ymid(nY(I)) >0               %%% Added Malik  % relaxed vertex or mid-plane
            n(I,:)=0;
            n(I,1)=nY(I)+nnodes;
        end
    end
end
for I=1:3 % loop on 3 vertices of triangle
    for i=1:nodes % Loop on nodes
        for J=1:3 % loop on 3 vertices of triangle
            for j=1:nodes % Loop on nodes
                if (n(I,i)>0 || i==1) && (n(J,j)>0 || j==1) && I~=J
                    idofg=((n(I,i)-1)*dim+1):n(I,i)*dim; % global dof
                    idofl=((I-1)*dim*nodes+(i-1)*dim)+1:(I-1)*dim*nodes+i*dim;
                    jdofg=((n(J,j)-1)*dim+1):n(J,j)*dim; % global dof
                    jdofl=((J-1)*dim*nodes+(j-1)*dim)+1:(J-1)*dim*nodes+j*dim;
                    for d=1:dim
                        si(sk+1:sk+dim)=idofg;
                        sj(sk+1:sk+dim)=jdofg(d);
                        sv(sk+1:sk+dim)=Kt(idofl,jdofl(d));
                        sk=sk+dim;
                    end
                end
            end
        end
    end
end
end
%
function Ke=KTriangleVol(nY,Ks,N,Yr,Ymid)
% Returns residual associated to a triangle, with contributions
% of each one of the 3 nodes/vertex, total 9 nodal contributions=27
% components.
% INPUT:
% N(I,:) = nodal shape function of vertex I
% nY = vertex number . nY=0 : "vertex" is a node
% Ks(I,J,,:,:) = matrix contribution associaed to vertices I and J
% Yr(j)   : >0, vertex j is relaxed
% OUTPUT:
% Ke  = volume Jacobian associated to face (81*dim*dim components, non-zero 54*dim*dim components)
nnodes=size(N,2); % Nodes per triangle
dim=3;
Ke=zeros(length(nY)*dim*nnodes);
for I=1:3 % loop on 3 vertices of triangle
    for i=1:nnodes % Loop on 3 nodes
        if nY(I)==0
            pi=1;
        elseif Yr(nY(I))>0 || Ymid(nY(I)) >0               %%% Added Malik
            pi=1;
            nY(I)=0;
        else
            pi=N(nY(I),i);
        end
        for J=1:3 % loop on 3 vertices of triangle
            for j=1:nnodes % Loop on 3 nodes
                if nY(J)==0
                    pj=1;
                elseif Yr(nY(J))>0 || Ymid(nY(J)) >0      %%% Added Malik
                    pj=1;
                    nY(J)=0;
                else
                    pj=N(nY(J),j);
                end
                dy=Ks((I-1)*dim+1:I*dim,(J-1)*dim+1:J*dim)*pi*pj; % 3x3 matrix
                if (nY(I)>0 || i==1) && (nY(J)>0 || j==1) && I~=J
                    idof=((I-1)*dim*nnodes+(i-1)*dim)+1:(I-1)*dim*nnodes+i*dim;
                    jdof=((J-1)*dim*nnodes+(j-1)*dim)+1:(J-1)*dim*nnodes+j*dim;
                    Ke(idof,jdof)=dy;
                end
            end
        end
    end
end
end
%%
function ge=gTriangleVol(nY,gs,N,Yr,Ymid)
% Returns residual associated to a triangle, with contributions
% of each one of the 3 nodes/vertex, total 9 nodal contributions=27
% components.
% INPUT:
% N(I,:) = nodal shape function of vertex I
% nY = vertex number . nY=0 : "vertex" is a node
% gs(I,:) = vector contribution associaed to vertex I
% Yr(j)   : >0, vertex j is relaxed
% OUTPUT:
% ge  = volume residual associated to face (9*dim components)
nnodes=size(N,2); % Nodes per triangle
dim=3;
ge=zeros(length(nY)*dim*nnodes,1);
for I=1:3 % loop on 3 vertices of triangle
    for n=1:nnodes % Loop on 3 nodes
        if nY(I)==0
            p=1;
        elseif Yr(nY(I))>0 || Ymid(nY(I)) >0               %%% Added Malik
            p=1;
            nY(I)=0;
        else
            p=N(nY(I),n);
        end
        dy=gs(I,:)*p; % I-th vertex, n-th node
        if nY(I)>0 || n==1 % if "vertex"=nodes or relaxed, just one loop
            idof=((I-1)*dim*nnodes+(n-1)*dim)+1:((I-1)*dim*nnodes+n*dim);
            ge(idof)=dy;
        end
    end
end
end


