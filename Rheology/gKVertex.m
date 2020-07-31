function [Ener,gy,Ky,Stress]=gKVertex(Ablated,Cell,Cv,Mat,L,L0,N,Y,Ymid,T,Set)
% Computes vertex residual/Jacobain
% T = triangle connectivity. Triangle number=Vertex number
% OUTPUT:x
% gy = vertex residual
% Ky = vertex Jacobain
dim=size(Y,2);
if Set.yRelaxation
    dimg=(Set.nodes+Set.nvert)*dim;
else
    dimg=Set.nodes*dim;
end
gy=zeros(dimg,1);
ncell=length(Cell);
bartot=Cell{1}.nElem;
Stress=zeros(bartot,1);
Ener=zeros(bartot,3); % energy due to [Elastic, tissue comtract, Wound contraactility]
if Set.Sparse
    sk=0;
    ndof=dim^2*(Set.nvert*10+(2+36)*bartot); % Assuming relaxed and unrelaxed. 10=1 (nodes) + 3x3 (nodal)
    si=zeros(ndof,1);
    sj=zeros(ndof,1);
    sv=zeros(ndof,1);
else
    Ky=zeros(dimg);
end
nv=0;
for c=1:ncell
    k=Mat.k;
    k0=Mat.k0;
    if Ablated.Cell{c}.Ablated
        k=k*Set.ReductionK;
        k0=k0*Set.ReductionK;
    end
    %-------------------------Malik Added begin ---------------------------
    cellCv=Cell{c}.CCv;
    eType=Cell{c}.eType;
    for e=1:length(cellCv)
        nv=cellCv(e);
        IJ=Cv(nv,:);
        etype=eType(e);
        IJr=Ablated.Yr(IJ);
        IJrZ=Ablated.YrZ(IJ);
        ye=[Y(IJ(1),:) ; Y(IJ(2),:)];
        wound=Ablated.Cell{c}.eType(e);
         if Ymid(IJ(1))>0
             IJr(1)=IJ(1);
         end 
         if Ymid(IJ(2))>0
             IJr(2)=IJ(2);
         end
         if etype ==1
            EpsC=Mat.EpsCB;
         elseif etype ==2
             EpsC=Mat.EpsCT;
         elseif etype ==3
             EpsC=Mat.EpsCL;
         end 
        if Set.Sparse    %%%%% 
            [Ener(nv,:),gy,Stress(nv),si,sj,sv,sk]=gKAssembleVertexSparse(gy,IJ,IJr,IJrZ,k0,k,Mat,EpsC,L(nv),L0(nv),N,T,ye,Set,wound,si,sj,sv,sk);
        else
            [Ener(nv,:),gy,Ky,Stress(nv)]=gKAssembleVertex(gy,Ky,IJ,IJr,k0,k,Mat,EpsC,L(nv),L0(nv),N,T,vdim,ye,Set,wound);
        end  
    end
    %-------------------------Malik Added end -----------------------------


    
    
end
if Set.Sparse
    Ky=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg);
end
end
%%
function  [Ener,gy,Ky,Stress]=gKAssembleVertex(gy,Ky,IJ,IJr,IJrZ,k0,k,Mat,epsC,L,L0,N,T,ye,Set,wound)
% epsC : contractility
dim=size(ye,2);
vdim=size(T,2)*dim; % Nodal dimesions of each vertex
[ge,Ke,Stress,Ener]=gkElemVor(k0,k,Mat,epsC,L,L0,IJ,IJr,IJrZ,L0,N,vdim,ye,Set,wound);
aux1=zeros(1,vdim);
aux2=aux1;
for d=1:dim
    aux1(d:dim:vdim)=dim*T(IJ(1),:)-dim+d;
    aux2(d:dim:vdim)=dim*T(IJ(2),:)-dim+d;
end
if IJr(1)>0
    aux1=[((IJr(1)-1)*dim+1:IJr(1)*dim)+Set.dim*Set.nodes zeros(1,2*dim)];
elseif IJrZ(1)>0
    aux1(dim:dim:vdim)=[IJrZ(1)*dim+Set.dim*Set.nodes 0 0];
end
if IJr(2)>0
    aux2=[((IJr(2)-1)*dim+1:IJr(2)*dim)+Set.dim*Set.nodes zeros(1,2*dim)];
elseif IJrZ(2)>0
    aux2(dim:dim:vdim)=[IJrZ(2)*dim+Set.dim*Set.nodes 0 0];
end
aux=[aux1 aux2];
for i=1:2*vdim
    if aux(i)>0
        gy(aux(i))=gy(aux(i))+ge(i);
        for j=1:2*vdim
            if aux(j)>0
                Ky(aux(i),aux(j))=Ky(aux(i),aux(j))+Ke(i,j);
            end
        end
    end
end
end
%%
function  [Ener,gy,Stress,si,sj,sv,sk]=gKAssembleVertexSparse(gy,IJ,IJr,IJrZ,k0,k,Mat,epsC,L,L0,N,T,ye,Set,wound,si,sj,sv,sk)
% epsC : contractility
dim=size(ye,2);
vdim=size(T,2)*dim; % Nodal dimesions of each vertex
[ge,Ke,Stress,Ener]=gkElemVor(k0,k,Mat,epsC,L,L0,IJ,IJr,IJrZ,N,vdim,ye,Set,wound);
aux1=zeros(1,vdim);
aux2=aux1;

if IJr(1)>0
    aux1=[((IJr(1)-1)*dim+1:IJr(1)*dim)+Set.dim*Set.nodes zeros(1,2*dim)];
elseif IJrZ(1)>0
    aux1(dim:dim:vdim)=[IJrZ(1)*dim+Set.dim*Set.nodes 0 0];
else 
    for d=1:dim
        aux1(d:dim:vdim)=dim*T(IJ(1),:)-dim+d;
    end
end

if IJr(2)>0
    aux2=[((IJr(2)-1)*dim+1:IJr(2)*dim)+Set.dim*Set.nodes zeros(1,2*dim)];
elseif IJrZ(2)>0
    aux2(dim:dim:vdim)=[IJrZ(2)*dim+Set.dim*Set.nodes 0 0];
else 
    for d=1:dim
        aux2(d:dim:vdim)=dim*T(IJ(2),:)-dim+d;
    end    
end

aux=[aux1 aux2];
for i=1:2*vdim
    if aux(i)>0
        gy(aux(i))=gy(aux(i))+ge(i);
        for j=1:2*vdim
            if aux(j)>0
                sk=sk+1;
                si(sk)=aux(i);
                sj(sk)=aux(j);
                sv(sk)=Ke(i,j);
            end
        end
    end
end
end
%%
% Elemental function
function [ge,Ke,Stress,Ener]=gkElemVor(k0,k,Mat,epsC,L,L0,IJ,IJr,IJrZ,N,vdim,ye,Set,wound)
% Ec0   : Material contractility
% IJ     : vertex numbers
dim=size(ye,2);
if IJr(1)>0
    ee1=kron([1 0 0],eye(dim));
elseif IJrZ(1)>0
    ee1=kron(N(IJ(1),:),eye(dim));
    ee1(dim,:)=[0 0 1 zeros(1,2*dim)];
else
    ee1=kron(N(IJ(1),:),eye(dim));
end
if IJr(2)>0
    ee2=kron([1 0 0],eye(dim));
elseif IJrZ(2)>0
    ee2=kron(N(IJ(2),:),eye(dim));
    ee2(dim,:)=[ 0 0 1 zeros(1,2*dim)];
else
    ee2=kron(N(IJ(2),:),eye(dim));
end
ee=[ee1 zeros(dim,vdim)
    zeros(dim,vdim) ee2];
[ge1,Ke1,Stress1,ener1]=gKElastic(k0,Mat,epsC,L0,ye,Set,wound,1); % Spring branch
[ge2,Ke2,Stress2,ener2]=gKElastic(k ,Mat,epsC,L ,ye,Set,wound,2); % Active branch
Ke=ee'*(Ke1+Ke2)*ee;
ge=ee'*(ge1+ge2);
Stress=Stress1+Stress2;
Ener=ener1+ener2;
end