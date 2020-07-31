function [Ener,g,K,s]=gKNodal(Ablated,X,C,Mat,L,L0,Set)
% OUTPUT:
% g  = global residual
% K  = global JAcobian
% s  = internal elemental forces
dim=size(X,2);
nodes=Set.nodes;
nvert=Set.nvert;
nelem=size(C,1);
s=zeros(nelem,1);
Ener=zeros(nelem,3); % energy due to [Elastic, tissue comtract, Wound contraactility]
if Set.yRelaxation || Set.midY                              %%% Added Malik
    dimg=dim*(nodes+nvert);
else
    dimg=dim*nodes;
end
g=zeros(dimg,1);
if Set.Sparse
    sk=0;
    si=zeros(dim*nodes*dim*nodes,1);
    sj=zeros(dim*nodes*dim*nodes,1);
    sv=zeros(dim*nodes*dim*nodes,1);
else
    K=zeros(dimg);
end
for e=1:nelem
    wound=Ablated.Nodal(e);
    k=Mat.k;
    k0=Mat.k0;
    if wound>0
      k=k*Set.ReductionK;
      k0=k0*Set.ReductionK;
    end
    [ge1,Ke1,s1,e1]=gKElastic(k0,Mat,Mat.EpsC,L0(e),X(C(e,:),:),Set,wound,1); % Spring branch
    [ge2,Ke2,s2,e2]=gKElastic(k ,Mat,Mat.EpsC,L(e),X(C(e,:),:),Set,wound,2); % Active branch
    s(e)=s1+s2;
    Ener(e,:)=e1+e2;
    if Set.Sparse
        [g,si,sj,sv,sk]=AssembleSparse(g,ge1+ge2,Ke1+Ke2,C(e,:),si,sj,sv,sk);
    else
        [g,K]=Assemble(g,ge1+ge2,K,Ke1+Ke2,C(e,:));
    end
end
if Set.Sparse
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimg,dimg);
end
end
%% Assembling
function [g,K]=Assemble(g,ge,K,Ke,Ce)
nodes=length(Ce); % nodes=2
dim=length(ge)/nodes;
for i=1:nodes
    idof=(i-1)*dim+1:i*dim;
    inode=Ce(i);
    igdof=(inode-1)*dim+1:dim*inode;
    g(igdof)=...
        g(igdof)+ge(idof);
    for j=1:nodes
        jdof=(j-1)*dim+1:j*dim;
        jnode=Ce(j);
        jgdof=(jnode-1)*dim+1:dim*jnode;
        K(igdof,jgdof)=...
            K(igdof,jgdof) + Ke(idof,jdof);
    end
end
end
%% Assembling Sparse
function [g,si,sj,sv,sk]=AssembleSparse(g,ge,Ke,Ce,si,sj,sv,sk)
nodes=length(Ce); % nodes=2
dim=length(ge)/nodes;
for i=1:nodes
    idof=(i-1)*dim+1:i*dim;
    inode=Ce(i);
    igdof=(inode-1)*dim+1:dim*inode;
    g(igdof)=g(igdof)+ge(idof);
    for j=1:nodes
        jnode=Ce(j);
        jdof=(j-1)*dim+1:j*dim;
        jgdof=(jnode-1)*dim+1:dim*jnode;
        for d=1:dim
            si(sk+1:sk+dim)=igdof(d);
            sj(sk+1:sk+dim)=jgdof;
            sv(sk+1:sk+dim)=Ke(idof(d),jdof);
            sk=sk+dim;
        end
    end
end
end
