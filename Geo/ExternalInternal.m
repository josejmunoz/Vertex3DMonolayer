function [xInternal,xExternal]=ExternalInternal(T,Edges)
% Find list of external and internal nodes according to 
% triangulation in T
% INPUT:
% T(i,:)       : Nodes forming triangle i
% Edges(i,s)   : Triangle adjacent to side s of triangle i. Side s is
% opposite to node s. if Esges(i,s)==0, side is external
% OUTPUT:
% xInternal=list of internal nodes
% xExternal=list of external nodes
nnodes=max(max(T));
xInternal=1:nnodes;
Sides=[2 3
       3 1
       1 2];
for i=1:size(T,1)
    for s=1:3
        if Edges(i,s)==0
            xInternal(T(i,Sides(s,:)))=0;
        end
    end
end
xExternal=1:nnodes;
xExternal=xExternal(xInternal==0);
xInternal(xInternal==0)=[];
end