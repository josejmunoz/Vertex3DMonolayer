function [Cell,Cv,Y,Yn,Y0,Ymid,Ablated,T]=ReorderVertices(Cell,Cv,Y,Yn,Y0,Ymid,Ablated,T,AbVBot,AbVTop,AbVmid)
%% Reorder ablated and wounded vertices\Triangles Y=[Y;Yb]

CellOld=Cell;
CvOld=Cv;
VRingBot1=VRingBot(:,1);
VRingTop1=VRingTop(:,1);
VRingBotOld1=VRingBot(:,1);
VRingTopOld1=VRingTop(:,1);
AV=Ablated.Yr(Ablated.Yr>0);
AY=Y(AV,:); Y(AV,:)=[]; 
nAV=(1:length(AV))+size(Y,1);
AbVBotOld=AbVBot; AbVTopOld=AbVTop; AbVmidOld=AbVmid;

for i=1:length(AV)
    Cv(CvOld>AV(i))=Cv(CvOld>AV(i))-1;
    Cv(CvOld==AV(i))=nAV(i);
    VRingBot1(VRingBotOld1>AV(i))=VRingBot1(VRingBotOld1>AV(i))-1;
    VRingBot1(VRingBotOld1==AV(i))=nAV(i);
    VRingTop1(VRingTopOld1>AV(i))=VRingTop1(VRingTopOld1>AV(i))-1;
    VRingTop1(VRingTopOld1==AV(i))=nAV(i);
    
    AbVBot(AbVBotOld>AV(i))=AbVBot(AbVBotOld>AV(i))-1;
    AbVBot(AbVBotOld==AV(i))=nAV(i);
    AbVTop(AbVTopOld>AV(i))=AbVTop(AbVTopOld>AV(i))-1;
    AbVTop(AbVTopOld==AV(i))=nAV(i);
    AbVmid(AbVmidOld>AV(i))=AbVmid(AbVmidOld>AV(i))-1;
    AbVmid(AbVmidOld==AV(i))=nAV(i);
    
    for ncell=1:length(Cell)
        Cell{ncell}.vbot(CellOld{ncell}.vbot>AV(i))= Cell{ncell}.vbot(CellOld{ncell}.vbot>AV(i))-1;
        Cell{ncell}.vbot(CellOld{ncell}.vbot==AV(i))= nAV(i);
        Cell{ncell}.vtop(CellOld{ncell}.vtop>AV(i))= Cell{ncell}.vtop(CellOld{ncell}.vtop>AV(i))-1;
        Cell{ncell}.vtop(CellOld{ncell}.vtop==AV(i))= nAV(i);
        Cell{ncell}.Tri(CellOld{ncell}.Tri>AV(i))= Cell{ncell}.Tri(CellOld{ncell}.Tri>AV(i))-1;
        Cell{ncell}.Tri(CellOld{ncell}.Tri==AV(i))= nAV(i);
    end 
end 

AT=T(AV,:); T(AV,:)=[]; T=[T;AT];
AY0=Y0(AV,:); Y0(AV,:)=[]; Y0=[Y0;AY0];
AYn=Yn(AV,:); Yn(AV,:)=[]; Yn=[Yn;AYn];
AYmid=Ymid(AV); Ymid(AV)=[]; Ymid=[Ymid; AYmid]; 
Ablated.Yr(AV)=[]; Ablated.Yr=[Ablated.Yr; nAV'];
% N doesb't have to be changed the all the same 
Y=[Y;AY];

VRingBot=[VRingBot1 VRingBot(:,[2 3])];
VRingTop=[VRingTop1 VRingTop(:,[2 3])];