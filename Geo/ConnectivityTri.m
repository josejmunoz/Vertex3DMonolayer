function [nTriNodeB,TriNodeB,nTriNodeT,TriNodeT]=ConnectivityTri(Cb,nnodes,Ct)                 
% or                [nTriNodeB,TriNodeB]=ConnectivityTri(Cb,nnodes)
% INPUT:
% Ct(i,:) : List of nodes forming triangle i
% OUTPUT:
% TriNode(i,:)  = list of triangles that nodes  i is connected to.
% nTriNode(i)   = number of triangles node i is connected t

ntriB=size(Cb,1); % = number of triangles Bottom
nnodes=nnodes/2;  % number of nodes in each layer
nTriNodeB=zeros(nnodes,1);
TriNodeB=zeros(nnodes,10);

for i=1:ntriB
    nTriNodeB(Cb(i,:))=nTriNodeB(Cb(i,:))+1;
    for s=1:3
        TriNodeB(Cb(i,s),nTriNodeB(Cb(i,s)))=i;
    end
end

if nargout > 2
    ntriT=size(Ct,1); % = number of triangles Top
    nTriNodeT=zeros(nnodes,1);
    TriNodeT=zeros(nnodes,10);
    for i=1:ntriT
        nTriNodeT(Ct(i,:)-nnodes)=nTriNodeT(Ct(i,:)-nnodes)+1;
        for s=1:3
            TriNodeT(Ct(i,s)-nnodes,nTriNodeT(Ct(i,s)-nnodes))=i+size(Cb,1);
        end
    end
end 

end