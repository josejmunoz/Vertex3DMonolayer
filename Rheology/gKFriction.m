function [gf,Kf]=gKFriction(Cell,N,T,Y,Yn,Set)
% Set.nu=0.1;
% OUTPUT:
% g  = global residual
% K  = global JAcobian
% s  = internal elemental forces
dim=size(Y,2);
nodes=Set.nodes;
nvert=Set.nvert;
if Set.yRelaxation || Set.midY                              %%% Added Malik
    dimg=dim*(nodes+nvert);
else
    dimg=dim*nodes;
end
gf=zeros(dimg,1);
Kf=zeros(dimg);


for n=1:length(Cell)
    vb=Cell{n}.vbot;
    Yb=Y(vb,:);
    Ybn=Yn(vb,:);
    y=reshape(Yb',[],1); % row of displacements
    yn=reshape(Ybn',[],1);
    
    Tv=T(vb,:);
    Tv=reshape(Tv',[],1);
    dofs=dim.*(kron(Tv,[1;1;1])-1);
    dofs=kron(ones(length(Tv),1),[1;2;3])+dofs;
    
    
    gc=(Set.eta/Set.dt).*(y-yn);
    Kc=(Set.eta/Set.dt).*ones(size(y));
    
    ee=[];
    for i=1:length(vb)
         eev=kron(N(vb(i),:),eye(dim));
         ee=blkdiag(ee,eev);
    end 

    gc=ee'*gc;
    Kc=ee'*diag(Kc)*ee;
    
    % Assembly
    for i=1:length(dofs)
        gf(dofs(i))=gf(dofs(i))+gc(i);
        for j=1:length(dofs)
            Kf(dofs(i),dofs(j))=Kf(dofs(i),dofs(j))+Kc(i,j);
        end 
    end 
    
    
end 