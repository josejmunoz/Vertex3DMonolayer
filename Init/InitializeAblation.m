function [Ablated,Set]=InitializeAblation(Cell,C,Cv,Set,Y)
% INPUT
% Cell = Cell database
% C    = connectivity
% X0   = Nodal positions
% Set  = Settings database
% OUTPUT
% Ablated:
%   Ablated.Nodal(e)=1, ablated nodal element e, =0: Non ablated nodal element
%   Ablated.Cell(e)=1, ablated cell e, =0: Non ablated cell
%   Ablated.Cv(e)=0, vertex bar element e is not ablated
%   Ablated.dofY   , list of dof corresponding to relaxed vertex nodes
%   (Startingg from 1)
% Cell{c}.vtopR(j)=j: top vertex j is relaxed
%                 =0: top vertex j is not relaxed (it is attached to nodes)
% Cell{c}.vbotR(j)=j: bottom vertex j is relaxed
%                 =0: bottom vertex j is not relaxed (it is attached to nodes)
% Set.Ablation = cell number being ablated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute total number of vertex elements:

nBar=size(Cv,1);
nelem=size(C,1);
Ablated.Yr=zeros(size(Y,1),1); % List of relaxed vertices
Ablated.YrZ=zeros(size(Y,1),1); % List of relaxed vertices only on Z
Set.Ablation=[];
Ablated.Nodal=zeros(nelem,1);
Ablated.Cv=zeros(nBar,1);
Ablated.Exist=false;
ncell=length(Cell);
Ablated.Cell=cell(ncell,1);
ncell=length(Cell);
for c=1:ncell
     Ablated.Cell{c}.eType=zeros(1,length(Cell{c}.eType));  %% Malik Added
     Ablated.Cell{c}.Ablated=0; % Ablated cell
     Ablated.Cell{c}.wounded=0; % Ablated cell
     Cell{c}.vtopR=zeros(length(Cell{c}.vtop),1);
     Cell{c}.vbotR=zeros(length(Cell{c}.vbot),1);
end
Ablated.dofY=[];
Ablated.dofYZ=[];
Ablated.VRingBot=[];
Ablated.VRingTop=[];
Ablated.NRingBot=[];
Ablated.NRingTop=[];
Ablated.AbNodesBot=[];
Ablated.AbNodesTop=[];
% Initialize Wound Size
Ablated.Height=zeros(Set.Nincr,1);
Ablated.Volume=zeros(Set.Nincr,1);
Ablated.AreaTop=zeros(Set.Nincr,1);
Ablated.AreaBottom=zeros(Set.Nincr,1);
Ablated.AreaLateral=zeros(Set.Nincr,1);
Ablated.RemodelN=zeros(Set.Nincr,2);
Ablated.NRingB=zeros(Set.Nincr,1);

end


