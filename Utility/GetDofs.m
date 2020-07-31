function [dof]=GetDofs(Ablated,Set,dofP,Ymid)
% Defines degress of freedom for nodes and vertices, takinginto accoutn
% constraints in dofP.

dof=1:Set.nodes*Set.dim;
dof(dofP)=[];

dofYm=3.*(kron(Ymid(Ymid>0),[1;1;1])-1);
dofYm=kron(ones(length(Ymid(Ymid>0)),1),[1;2;3])+dofYm;
dofYm=dofYm+Set.nodes*Set.dim;

dof=[dof unique([Ablated.dofY' Ablated.dofYZ'])+Set.nodes*Set.dim dofYm'];
dof=unique(dof); % to remove repetition of wounded dofs and mid-plane dofs

end 