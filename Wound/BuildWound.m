function [Ablated,Cell,X,T,C,Cv,Y,N,L,Ln,L0,Stress,X0,Y0,Xn,Yn,xExternal,Set,x,dof,Ymid]=BuildWound...
                 (Ablated,Cell,X,T,C,Cv,Y,N,L,Ln,L0,Stress,X0,Y0,Xn,Yn,xExternal,Set,Ymid)
%% Updates nodal and vertex database (X;y) and wound database Ablated, Cell 
% when wound is created.
             
%% Define Ablated cells              
[Ablated,Set]=Ablation(Ablated,Cell,C,Cv,X0,Set,Y);
% function to defind wound ring and remove wounded entities
Set.EpsCBW_N=0;
Set.EpsCTW_N=0;
Set.EpsCLW_N=0;
Set.EpsCTW=0;  
Set.EpsCBW=0; 
Set.EpsCLW=0;

%% find top\bottom Ablated nodes and nodes connected to them
[AbNBot,AbNTop,cNBot,cNTop]=FindWoundNodes(Cell,Set.Ablation);

%% build nodal wound ring
[NRingBot,NRingTop]=NodalWoundRing(cNBot,cNTop,AbNBot,AbNTop,T);

%% build vertex wound ring 
[VRingBot,VRingTop]=VertexWoundRing(NRingBot,NRingTop,AbNBot,AbNTop,T);

%% Find Top and bottom ablated (removed) vertices 
[AbVBot,AbVTop,AbVmid]=FindAblatedVertices(Ablated,Cv,AbNBot,AbNTop,T,Set.ModelTop);

%% Remove ablated entities and Update Model parameter 

% -Remove dead Cells
Ablated.Cell(Set.Ablation)=[];
for i=1:length(Cell)
    for j=1:length(Set.Ablation)
        Cell{i}.nTri=Cell{i}.nTri-size(Cell{Set.Ablation(j)}.Tri,1)...
                                  -length(Cell{Set.Ablation(j)}.vbot)...
                                  -length(Cell{Set.Ablation(j)}.vtop);
    end 
end 
Cell(Set.Ablation)=[];
 % - Nodal
 BotDeadTri=any(ismember(T,AbNBot),2);
 TopDeadTri=any(ismember(T,AbNTop),2);

InTri=1:size(N,1);
InTri=InTri(BotDeadTri | TopDeadTri);
N(InTri,:)=zeros(length(InTri),3);  
T(InTri,:)=zeros(length(InTri),3);  
X([AbNBot ; AbNTop],:)=[];
X0([AbNBot ; AbNTop],:)=[];
Xn([AbNBot ; AbNTop],:)=[];

BotDeadBarN=any(ismember(C,AbNBot),2);
TopDeadBarN=any(ismember(C,AbNTop),2);
L.D.L(BotDeadBarN | TopDeadBarN,:)=[]; 
Ln.D.L(BotDeadBarN | TopDeadBarN,:)=[]; 
L0.D.L(BotDeadBarN | TopDeadBarN,:)=[]; 
Stress.D.s(BotDeadBarN | TopDeadBarN)=[];
Ablated.Nodal(BotDeadBarN | TopDeadBarN,:)=[]; %??
C(BotDeadBarN | TopDeadBarN,:)=[];

VRingBot23=VRingBot(:,[2 3]);
VRingTop23=VRingTop(:,[2 3]);

COld=C;
TOld=T;
CellOld=Cell;
xExternalOld=xExternal;
VRingBot23Old=VRingBot23;
VRingTop23Old=VRingTop23;
NRingBotOld=NRingBot;
NRingTopOld=NRingTop;
AbN=[AbNBot;AbNTop];

for i=1:length(AbN)
    T(TOld>AbN(i))=T(TOld>AbN(i))-1;
    C(COld>AbN(i))=C(COld>AbN(i))-1;
    xExternal(xExternalOld>AbN(i))=xExternal(xExternalOld>AbN(i))-1;
    VRingBot23(VRingBot23Old>AbN(i))=VRingBot23(VRingBot23Old>AbN(i))-1;
    VRingTop23(VRingTop23Old>AbN(i))=VRingTop23(VRingTop23Old>AbN(i))-1;
    NRingBot(NRingBotOld>AbN(i))=NRingBot(NRingBotOld>AbN(i))-1;
    NRingTop(NRingTopOld>AbN(i))=NRingTop(NRingTopOld>AbN(i))-1;
    for ncell=1:length(Cell)
        Cell{ncell}.xBott(CellOld{ncell}.xBott>AbN(i))=Cell{ncell}.xBott(CellOld{ncell}.xBott>AbN(i))-1;
        Cell{ncell}.xTopp(CellOld{ncell}.xTopp>AbN(i))=Cell{ncell}.xTopp(CellOld{ncell}.xTopp>AbN(i))-1;
        Cell{ncell}.xBot(CellOld{ncell}.xBot>AbN(i))=Cell{ncell}.xBot(CellOld{ncell}.xBot>AbN(i))-1;
        Cell{ncell}.xTop(CellOld{ncell}.xTop>AbN(i))=Cell{ncell}.xTop(CellOld{ncell}.xTop>AbN(i))-1;
    end 
end

% - vertex
DeadBarV=Ablated.Cv~=0;
IDBars=1:size(Cv,1);
IDBars=IDBars(DeadBarV);
L.V.L(DeadBarV,:)=[]; 
L0.V.L(DeadBarV,:)=[]; 
Ln.V.L(DeadBarV,:)=[]; 
Stress.V.s(DeadBarV)=[];
Cv(DeadBarV,:)=[];
Ablated.Cv(DeadBarV,:)=[];

AbV=[AbVBot'; AbVTop'; AbVmid];

Y(AbV,:)=[];
Y0(AbV,:)=[];
Yn(AbV,:)=[];
Ymid(AbV)=[];
N(AbV,:)=[];  
T(AbV,:)=[];

Ablated.Yr(AbV)=[];   %??
Ablated.YrZ(AbV)=[];  %??

VRingBot1=VRingBot(:,1);
VRingTop1=VRingTop(:,1);

CvOld=Cv;
VRingBot1Old=VRingBot1;
VRingTop1Old=VRingTop1;
AblatedYrOld=Ablated.Yr;
YmidOld=Ymid;
for i=1:length(AbV)
    Cv(CvOld>AbV(i))=Cv(CvOld>AbV(i))-1;
    Ablated.Yr(AblatedYrOld>AbV(i))=Ablated.Yr(AblatedYrOld>AbV(i))-1;
    Ymid(YmidOld>AbV(i))=Ymid(YmidOld>AbV(i))-1;
    VRingBot1(VRingBot1Old>AbV(i))=VRingBot1(VRingBot1Old>AbV(i))-1;
    VRingTop1(VRingTop1Old>AbV(i))=VRingTop1(VRingTop1Old>AbV(i))-1;
    for ncell=1:length(Cell)
        Cell{ncell}.vbot(CellOld{ncell}.vbot>AbV(i))=Cell{ncell}.vbot(CellOld{ncell}.vbot>AbV(i))-1;
        Cell{ncell}.vtop(CellOld{ncell}.vtop>AbV(i))=Cell{ncell}.vtop(CellOld{ncell}.vtop>AbV(i))-1;
        Cell{ncell}.Tri(CellOld{ncell}.Tri>AbV(i))=Cell{ncell}.Tri(CellOld{ncell}.Tri>AbV(i))-1;
    end
end 
for ncell=1:length(Cell)
    for j=1:length(IDBars)
        Cell{ncell}.CCv(CellOld{ncell}.CCv>IDBars(j))=Cell{ncell}.CCv(CellOld{ncell}.CCv>IDBars(j))-1;
    end
end 
VRingBot=[VRingBot1 VRingBot23];
VRingTop=[VRingTop1 VRingTop23];

%% - Cell
Cell{1}.nElem=size(Cv,1);

if Set.FixedBasal==0
    dofY=3.*(kron([VRingBot(:,1);VRingTop(:,1)],[1;1;1])-1);
    dofY=kron(ones(length(VRingBot(:,1))+length(VRingTop(:,1)),1),[1;2;3])+dofY;
    Ablated.dofY=unique(dofY);
else
    dofYBot=3.*(kron(VRingBot(:,1),[1;1])-1);
    dofYBot=kron(ones(length(VRingBot(:,1)),1),[1;2])+dofYBot;

    dofYTop=3.*(kron(VRingTop(:,1),[1;1;1])-1);
    dofYTop=kron(ones(length(VRingTop(:,1)),1),[1;2;3])+dofYTop;
    dofY=[dofYBot; dofYTop];
    Ablated.dofY=unique(dofY);
end 
Set.nvert=size(Y,1);
Set.nodes=size(X,1);

x=[reshape(X',Set.nodes*Set.dim,1);reshape(Y',Set.nvert*Set.dim,1)]; % row of displacements

[~,dofP,~]=BC(Set,X,xExternal);

[dof]=GetDofs(Ablated,Set,dofP,Ymid);

Ablated.VRingBot=VRingBot;
Ablated.VRingTop=VRingTop;
Ablated.NRingBot=NRingBot(1:end-1);
Ablated.NRingTop=NRingTop(1:end-1);
Ablated.AbNodesBot=AbNBot;
Ablated.AbNodesTop=AbNTop;

%% Build wound lateral triangles and redefine wounded bars

YYmid=Ymid(Ymid>0);
wYmid=zeros(size(YYmid));
% find mid-plane wounded vertices
for c=1:ncell 
    wTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1);YYmid]),2) |...
            all(ismember(Cell{c}.Tri,[Ablated.VRingTop(:,1);YYmid]),2);
    wTri=Cell{c}.Tri(wTriRow,:); 
    wTri=reshape(wTri,[],1);
    wYmidRow=ismember(YYmid,wTri);
    wYmid(wYmidRow)=YYmid(wYmidRow);
end

WTri=zeros(3*length(Ablated.VRingBot(:,1)),3);
Tnt=0;

for c=1:ncell     
    unWCCv=Cv(Cell{c}.CCv,:);
    Ablated.Cell{c}.eType=zeros(size(Cell{c}.CCv));
    Ablated.Cell{c}.wounded=0;
    Ablated.Cell{c}.Ablated=0; % Ablated cell
    for jj=1:size(unWCCv,1) % loop over the elements of the unwounded cell (c) 
        if all(ismember(unWCCv(jj,:),[Ablated.VRingBot(:,1); Ablated.VRingTop(:,1); wYmid]));     
            Ablated.Cell{c}.eType(jj)=Cell{c}.eType(jj)+6;
            Ablated.Cell{c}.wounded=1;
        end
    end
    
    if Ablated.Cell{c}.wounded==1
        wwTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1); Ablated.VRingTop(:,1); wYmid]),2);
        wwTri=Cell{c}.Tri(wwTriRow,:); wwTri=flip(wwTri,2);
        nt=size(wwTri,1);
        WTri(Tnt+1:Tnt+nt,:)=wwTri;
        Tnt=Tnt+nt;
    end
end

WTri(Tnt+1:end,:)=[];
Ablated.WTri=WTri;
end 