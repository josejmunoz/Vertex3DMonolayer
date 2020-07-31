function [gp]=gPropulsion(Cell,Set)




dim=3;
nodes=Set.nodes;
nvert=Set.nvert;
if Set.yRelaxation || Set.midY                              %%% Added Malik
    dimg=dim*(nodes+nvert);
else
    dimg=dim*nodes;
end
gp=zeros(dimg,1);



for n=1:length(Cell)
    if ismember(n,Set.PropulsiveCells.Region1)
        dof=dim.*(kron(Cell{n}.xBott ,[1;1;1])-1)+[1;2;3]; 
        gp(dof)=gp(dof)+Set.mu1';
    elseif ismember(n,Set.PropulsiveCells.Region2)
        dof=dim.*(kron(Cell{n}.xBott ,[1;1;1])-1)+[1;2;3]; 
        gp(dof)=gp(dof)+Set.mu2';
    end 
end 