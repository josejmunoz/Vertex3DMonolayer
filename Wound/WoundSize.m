function Ablated=WoundSize(Ablated,Cell,C,Cv,X0,Set,Y,X)
%% Updates wound quantification for current increment:
% Ablated.Height(Set.iIncr)
% Ablated.Volume(Set.iIncr)
% Ablated.AreaTop(Set.iIncr)
% Ablated.AreaBottom(Set.iIncr)
% Ablated.AreaLateral(Set.iIncr)
%
Vw=0; Atw=0; Abw=0; Alw=0;
if Ablated.Exist==0
    [~,FakeSet]=Ablation(Ablated,Cell,C,Cv,X0,Set,Y);
    for i=1:length(FakeSet.Ablation)
        Vw=Vw+Cell{FakeSet.Ablation(i)}.Vol;
        
        vbot=Cell{FakeSet.Ablation(i)}.vbot;
        vbot=[vbot;vbot(1)];
        for j=1:size(vbot,1)-1
            Yi=Y(vbot(j),:);
            Yj=Y(vbot(j+1),:);
            a=[Yi;Yj;X(Cell{FakeSet.Ablation(i)}.xBott,:)];
            Vw=Vw+det(a)/6;
            Abw=Abw+polyarea(a(:,1),a(:,2));  % Top Area 
        end
        
        vtop=Cell{FakeSet.Ablation(i)}.vtop;
        vtop=[vtop;vtop(1)];
        for j=1:size(Cell{FakeSet.Ablation(i)}.vtop,1)-1
            Yi=Y(vtop(j),:);
            Yj=Y(vtop(j+1),:);
            a=[Yi;Yj;X(Cell{FakeSet.Ablation(i)}.xTopp,:)];
            Vw=Vw+det(a)/6;
            Atw=Atw+polyarea(a(:,1),a(:,2));  % Top Area 
        end
    end 
    if Set.iIncr>1
        Ablated.Height(Set.iIncr)=Ablated.Height(Set.iIncr-1);
        if max(Ablated.Volume(1:Set.iIncr-1))>0
            Ablated.Volume(Set.iIncr)=0;
            Ablated.AreaTop(Set.iIncr)=0;
            Ablated.AreaBottom(Set.iIncr)=0;
            Ablated.AreaLateral(Set.iIncr)=0;
        else
            Ablated.Volume(Set.iIncr)=Vw;
            Ablated.AreaTop(Set.iIncr)=Atw;
            Ablated.AreaBottom(Set.iIncr)=Abw;
            Ablated.AreaLateral(Set.iIncr)=Alw;
        end
    else
        Ablated.Height(Set.iIncr)=0;
        Ablated.Volume(Set.iIncr)=Vw;
        Ablated.AreaTop(Set.iIncr)=Atw;
        Ablated.AreaBottom(Set.iIncr)=Abw;
        Ablated.AreaLateral(Set.iIncr)=Alw;
    end
    return
end 

Ybc=sum(Y(Ablated.VRingBot(:,1),:),1)/size(Ablated.VRingBot,1);
Ytc=sum(Y(Ablated.VRingTop(:,1),:),1)/size(Ablated.VRingBot,1);

% Wound Hight
Hw=Ytc(3)-Ybc(3);

% Bottom
YBot=[Y(Ablated.VRingBot(:,1),:);Y(Ablated.VRingBot(1,1),:)];
% Check orientation
v1=YBot(2,:)-YBot(1,:);
v2=YBot(3,:)-YBot(2,:);
if cross(v1,v2)*[0 0 1]'<0
    YBot=flip(YBot,1);
end

for i=1:size(YBot,1)-1
    Yi=YBot(i,:);
    Yj=YBot(i+1,:);
    a=[Yj;Yi;Ybc];
    Vw=Vw+det(a)/6;
    Abw=Abw+polyarea(a(:,1),a(:,2));  % Bottom Area 
end

% Top
YTop=[Y(Ablated.VRingTop(:,1),:);Y(Ablated.VRingTop(1,1),:)];
% Check orientation
v1=YTop(2,:)-YTop(1,:);
v2=YTop(3,:)-YTop(2,:);
if cross(v1,v2)*[0 0 1]'<0
    YTop=flip(YTop,1);
end
for i=1:size(YTop,1)-1
    Yi=YTop(i,:);
    Yj=YTop(i+1,:);
    a=[Yi;Yj;Ytc];
    Vw=Vw+det(a)/6;
    Atw=Atw+polyarea(a(:,1),a(:,2));  % Top Area 
end


% Lateral
for t=1:size(Ablated.WTri,1)
    YTri=Y(Ablated.WTri(t,:),:);
    Vw=Vw+det(YTri)/6;
    Alw=Alw+polyarea(YTri(:,1),YTri(:,2));  % Top Area 
end 

Ablated.Height(Set.iIncr)=Hw;
Ablated.Volume(Set.iIncr)=Vw;
Ablated.AreaTop(Set.iIncr)=Atw;
Ablated.AreaBottom(Set.iIncr)=Abw;
Ablated.AreaLateral(Set.iIncr)=Alw;
Ablated.NRingB(Set.iIncr)=length(Ablated.NRingBot);
% Store cell volume per time step
for i=1:length(Cell)
    if ~isfield(Ablated.Cell{i},'Vol')
        Ablated.Cell{i}.Vol=zeros(Set.Nincr,1);
    end
    Ablated.Cell{i}.Vol(Set.iIncr)=Cell{i}.Vol;
end 
end 