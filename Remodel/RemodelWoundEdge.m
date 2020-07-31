function [Y,Ablated,Set,T,N,VNew2OldRing,Cb,Ct,Edges,YOld]=RemodelWoundEdge(Ablated,Y,Set,T,N,YOld,CellOld,X,Cb,Ct,Edges)
% Checks if wound needs remodelling, and updates Ablated and triangulation
% database accordingly
if ~Ablated.Exist
    VNew2OldRing=[];
    return
end
%% Remove zigzag verices        
VRingBot23=Ablated.VRingBot(:,[2 3]);
VRingTop23=Ablated.VRingTop(:,[2 3]);
Ablated.VRingBot(VRingBot23(:,2)==0,:)=[];
Ablated.VRingTop(VRingTop23(:,2)==0,:)=[];
%% Reorder
% VRingBot1old=Ablated.VRingBot(:,1);
% VRingTop1old=Ablated.VRingTop(:,1);
Yr=YOld(ismember(Ablated.Yr,[Ablated.VRingBot(:,1); Ablated.VRingTop(:,1)]),:);
IdB=1:length(Ablated.VRingBot(:,1));
IdT=1:length(Ablated.VRingTop(:,1));
[~,IdsB]=sort(Ablated.VRingBot(:,1));
[~,IdsT]=sort(Ablated.VRingTop(:,1));
VRingBot1new=zeros(length(Ablated.VRingBot(:,1)),1);
VRingTop1new=zeros(length(Ablated.VRingTop(:,1)),1);
VRingBot1new(IdsB)=IdB;
VRingTop1new(IdsT)=IdT+length(IdB);
if length(Ablated.NRingBot)==4 && length(Ablated.NRingTop)==4 
    ANB = polyarea(X(Ablated.NRingBot,1),X(Ablated.NRingBot,2));
    AVB = polyarea(YOld(Ablated.VRingBot(:,1),1),YOld(Ablated.VRingBot(:,1),2));
    ANT    = polyarea(X(Ablated.NRingTop,1),X(Ablated.NRingTop,2));
    AVT = polyarea(YOld(Ablated.VRingTop(:,1),1),YOld(Ablated.VRingTop(:,1),2));
    toll=0.05;  %%%%%%%%%%%%%%%%%%%%%%%
    if (AVB/ANB<toll || AVT/ANT<toll) 
        d1=norm(X(Ablated.NRingBot(1),:)-X(Ablated.NRingBot(3),:));
        d2=norm(X(Ablated.NRingBot(2),:)-X(Ablated.NRingBot(4),:));
        if d1<=d2
            tri1b=[Ablated.NRingBot(1) Ablated.NRingBot(2) Ablated.NRingBot(3)];
            tri2b=[Ablated.NRingBot(1) Ablated.NRingBot(3) Ablated.NRingBot(4)];
            tri1t=[Ablated.NRingTop(1) Ablated.NRingTop(2) Ablated.NRingTop(3)];
            tri2t=[Ablated.NRingTop(1) Ablated.NRingTop(3) Ablated.NRingTop(4)];
        else
            tri1b=[Ablated.NRingBot(1) Ablated.NRingBot(2) Ablated.NRingBot(4)];
            tri2b=[Ablated.NRingBot(2) Ablated.NRingBot(3) Ablated.NRingBot(4)];
            tri1t=[Ablated.NRingTop(1) Ablated.NRingTop(2) Ablated.NRingTop(4)];
            tri2t=[Ablated.NRingTop(2) Ablated.NRingTop(3) Ablated.NRingTop(4)];
        end
        Ablated.VRingBot=[];
        Ablated.VRingTop=[];
        Cb=[Cb; tri1b; tri2b];
        Ct=[Ct; tri1t; tri2t];
        [nTrianglesB,Y,nTrianglesT]=GetY(Cb,Ct,X);
        N=ones(nTrianglesB+nTrianglesT,3)/3;
        T=[Cb; Ct];
        Set.nvert=size(Y,1);
        Edges=ConnectivityEdges(Cb);
        Ablated.NRingBot=[];
        Ablated.NRingTop=[];
        VNew2OldRing=[];
        Ablated.Exist=false;
    return
    end   
end
%% Compute remodeling Indicator=[perB perL angle];
Indicator=zeros(length(CellOld),2); % Indicator=[perB perL angle];
for n=1:length(CellOld)
    if Ablated.Cell{n}.wounded==0
        Indicator(n,:)=[nan nan];
        continue
    end 
    % Cell perimeter
    vb=[CellOld{n}.vbot; CellOld{n}.vbot(1)];
    vt=[CellOld{n}.vtop; CellOld{n}.vtop(1)];
    perB=0;
    perT=0;
    for i=1:length(vb)-1
        perB=perB+norm(YOld(vb(i),:)-YOld(vb(i+1),:));
    end 
    for i=1:length(vt)-1
        perT=perT+norm(YOld(vt(i),:)-YOld(vt(i+1),:));
    end 
    % vertices on the wouned edge 
    rowWvvb=ismember(Ablated.VRingBot(:,1),vb);   % row logical 
    rowWvvt=ismember(Ablated.VRingTop(:,1),vt);
    Wvvb=Ablated.VRingBot(rowWvvb,1);             % vertix Id 
    Wvvt=Ablated.VRingTop(rowWvvt,1);
    if ~isempty(Wvvb)
        % Length of wounded element
        Lb=norm(YOld(Wvvb(1),:)-YOld(Wvvb(2),:));
        % bot-nodes the wound edge 
        ele=Ablated.VRingBot(rowWvvb,[2 3]);
        n1b=intersect(ele(1,:),ele(2,:));      % common node to be removed from wound ring 
        n2b=ele(1,~ismember(ele(1,:),n1b));    
        n3b=ele(2,~ismember(ele(2,:),n1b));
        % compute angle 
        [AngleB]=FindAngle(n1b,n2b,n3b,X); 
        if Set.ModelTop ==3
            Tri=[n1b n2b n3b]+max(max(Cb));
            check=CheckCompatibility(Tri,Ct,Ablated.VRingTop(:,[2 3]));
        else 
            check=true;
        end 

        if AngleB>130 || ~check
            Indicator(n,1)=nan;
        else 
            Indicator(n,1)=Lb/perB;
        end
    else
        Indicator(n,1)=nan;
    end 
    if ~isempty(Wvvt)
        % Length of wounded element
        Lt=norm(YOld(Wvvt(1),:)-YOld(Wvvt(2),:));    
        % top-nodes the wound edge 
        ele=Ablated.VRingTop(rowWvvt,[2 3]);
        n1t=intersect(ele(1,:),ele(2,:));      % common node to be removed from wound ring 
        n2t=ele(1,~ismember(ele(1,:),n1t));    
        n3t=ele(2,~ismember(ele(2,:),n1t));
        % compute angle 
        [AngleT]=FindAngle(n1t,n2t,n3t,X);
        if Set.ModelTop ==3
             Tri=[n1t n2t n3t]-max(max(Cb));
            [check]=CheckCompatibility(Tri,Cb,Ablated.VRingBot(:,[2 3]));
        else 
            check=true;
        end 

        if AngleT>130 || ~check
            Indicator(n,2)=nan;
        else 
            Indicator(n,2)=Lt/perT;
        end 
    else 
        Indicator(n,2)=nan;
    end 
end 
if Set.ModelTop == 3 
    [ValueB,cellIdB]=min(Indicator(:,1),[],'omitnan');
    [ValueT,cellIdT]=min(Indicator(:,2),[],'omitnan');
else 
    minis=min(Indicator,[],2,'omitnan');
    [ValueB,cellIdB]=min(minis);
    ValueT=ValueB;
    cellIdT=cellIdB;
end 
%%    %% ================================ Bottom
if ValueB<Set.WRemodelThreshold  && length(Ablated.NRingBot)>4
    vb=[CellOld{cellIdB}.vbot; CellOld{cellIdB}.vbot(1)];
    % vertices on the wouned edge 
    rowWvvb=ismember(Ablated.VRingBot(:,1),vb);   % row logical 
    Wvvb=Ablated.VRingBot(rowWvvb,1);             % vertix Id 
    if Wvvb(1)>Wvvb(2)
        Wvvb=flip(Wvvb);
    end 
    %-------- bulid tringle 
    ele=Ablated.VRingBot(rowWvvb,[2 3]);  % the corresponding nodal ele
    n1=intersect(ele(1,:),ele(2,:));      % common node to be removed from wound ring 
    n2=ele(1,~ismember(ele(1,:),n1));    
    n3=ele(2,~ismember(ele(2,:),n1));
    [n1,n2,n3]=OrientNodes(n1,n2,n3,X);
    TriB=[n1 n2 n3];                      % new tringle 
    Cb=[Cb;TriB];
    v1b=Wvvb(1);         % first vertex to be kept  
    v2b=Wvvb(2);         % secound to be removed 
       % --------find the new position 
    dX=X(n2,:)-X(n3,:);                   % reflection vector       
%     Ynew=(1/3).*sum(X(TriB,:),1);         % vertex to be reflected 
    Ynew=(3/7).*sum(X([n2 n3],:),1)+(1/7).*X(n1,:);
    Ynew=Ynew-X(n3,:);
    Yrefl=2*(dot(Ynew,dX)/dot(dX,dX))*dX-Ynew;   % refelect 
    Yrefl=Yrefl+X(n3,:);
%     Yrefl=ProjectY(Yrefl,Ablated.VRingBot(:,1),YOld);
    YOld(v1b,:)=Yrefl;                     % update the position
    %  ---------undate position 
    Yr(VRingBot1new(Ablated.VRingBot(:,1)==v1b),:)=Yrefl;
    % ----------update nodes
    Ablated.VRingBot(Ablated.VRingBot(:,1)==v1b,[2 3])=[n2 n3]; 
    % -------- remove node 
    Ablated.NRingBot(Ablated.NRingBot==n1)=[];
    % ------- vertex from Y 
    Yr(VRingBot1new(Ablated.VRingBot(:,1)==v2b),:)=[];
    % -------- remove vertex from new data 
    VRingBot1new(Ablated.VRingBot(:,1)==v2b,:)=[]; 
    % -------- remove vertex from Old data 
    Ablated.VRingBot(Ablated.VRingBot(:,1)==v2b,:)=[];  
    
    VRingBot1new(Ablated.VRingBot(:,1)>v2b)=VRingBot1new(Ablated.VRingBot(:,1)>v2b)-1;
    VRingTop1new(Ablated.VRingTop(:,1)>v2b)=VRingTop1new(Ablated.VRingTop(:,1)>v2b)-1; 
end 
%%   %% ================================ Top
if ValueT<Set.WRemodelThreshold  && length(Ablated.NRingTop)>4
    vt=[CellOld{cellIdT}.vtop; CellOld{cellIdT}.vtop(1)];
    % vertices on the wouned edge 
    rowWvvt=ismember(Ablated.VRingTop(:,1),vt);
    Wvvt=Ablated.VRingTop(rowWvvt,1);
    if Wvvt(1)>Wvvt(2)
        Wvvt=flip(Wvvt);
    end 
    %% ================================ Top
        %-------- bulid tringle 
    ele=Ablated.VRingTop(rowWvvt,[2 3]);  % the corresponding nodal ele
    n1=intersect(ele(1,:),ele(2,:));      % common node to be removed from wound ring 
    n2=ele(1,~ismember(ele(1,:),n1));    
    n3=ele(2,~ismember(ele(2,:),n1));
    [n1,n2,n3]=OrientNodes(n1,n2,n3,X);
    TriT=[n1 n2 n3];                      % new tringle 
    Ct=[Ct;TriT];
         % --------find the new position 
    v1t=Wvvt(1);         % first vertex to be kept  
    v2t=Wvvt(2);         % secound to be removed 
    dX=X(n2,:)-X(n3,:);                   % reflection vector       
%     Ynew=(1/3).*sum(X(TriT,:),1);         % vertex to be reflected
    Ynew=(3/7).*sum(X([n2 n3],:),1)+(1/7).*X(n1,:);
    Ynew=Ynew-X(n3,:);
    Yrefl=2*(dot(Ynew,dX)/dot(dX,dX))*dX-Ynew;   % refelect 
    Yrefl=Yrefl+X(n3,:);
%     Yrefl=ProjectY(Yrefl,Ablated.VRingTop(:,1),YOld);
    YOld(v1t,:)=Yrefl;                     % update the position
    %  ---------undate position 
    Yr(VRingTop1new(Ablated.VRingTop(:,1)==v1t),:)=Yrefl;
    
    % ----------update nodes
    Ablated.VRingTop(Ablated.VRingTop(:,1)==v1t,[2 3])=[n2 n3]; 
    
    % -------- remove node 
    Ablated.NRingTop(Ablated.NRingTop==n1)=[];
    
    % ------- vertex from Y 
    Yr(VRingTop1new(Ablated.VRingTop(:,1)==v2t),:)=[];
        % -------- remove vertex from new data 
    VRingTop1new(Ablated.VRingTop(:,1)==v2t,:)=[]; 
        % -------- remove vertex from Old data 
    Ablated.VRingTop(Ablated.VRingTop(:,1)==v2t,:)=[];
    
    VRingBot1new(Ablated.VRingBot(:,1)>v2t)=VRingBot1new(Ablated.VRingBot(:,1)>v2t)-1;
    VRingTop1new(Ablated.VRingTop(:,1)>v2t)=VRingTop1new(Ablated.VRingTop(:,1)>v2t)-1;
end 
    %% ======== update variables 
[nTrianglesB,Y,nTrianglesT]=GetY(Cb,Ct,X);
Edges=ConnectivityEdges(Cb);
VRingBot1new=VRingBot1new+nTrianglesB+nTrianglesT;
VRingTop1new=VRingTop1new+nTrianglesB+nTrianglesT;
Y=[Y;Yr];
T=[Cb; Ct];
N=ones(nTrianglesB+nTrianglesT,3)/3;
N=[N; zeros(size(Yr,1),3)];
T=[T; zeros(size(Yr,1),3)];
Set.nvert=size(Y,1);
VNew2OldRingBot=[VRingBot1new Ablated.VRingBot(:,1)];
VNew2OldRingTop=[VRingTop1new Ablated.VRingTop(:,1)];
VNew2OldRing=[VNew2OldRingBot;VNew2OldRingTop];
Ablated.VRingBot=[VRingBot1new Ablated.VRingBot(:,[2 3])];
Ablated.VRingTop=[VRingTop1new Ablated.VRingTop(:,[2 3])];
end 
%%
function [Angle]=FindAngle(n,i,j,X)
V1=X(j,:)-X(n,:);
V2=X(i,:)-X(n,:);
Angle=acosd(dot(V1,V2)/(norm(V1)*norm(V2)));
end
%%
function check=CheckCompatibility(Tri,C,VRing)
%%
 logic=all(ismember(C,Tri),2);
 if any(logic)
     check=true;
 else 
     logic=all(ismember(VRing,Tri),2);
     if sum(logic)==2
              check=true;
     else 
              check=false;
     end 
 end 
end 
%%
function Ynew=ProjectY(Yi,Yedge,Y)
% Yi    :  vertex to be projected
% Yedge :  ordered list of vertices that form the edge front
% Y(v,:): coordinates of vertex v
%
% Find two vertices closer to vertex
ne=length(Yedge);
Dist=zeros(ne,1);
for v=1:ne
    Dist(v)=norm(Yi-Y(Yedge(v),:));
end
[~,Order]=sort(Dist);
Ys=Yedge(Order(1:2)); % Two closer vertices
% Force that two consecutive vertices are taken
if Order(1)==ne
    Y2(1)=1;
else
    Y2(1)=Order(1)+1;
end
if Order(1)==1
    Y2(2)=ne;
else
    Y2(2)=Order(1)-1;
end
if Ys(2)~=Yedge(Y2(1)) && Ys(2)~=Yedge(Y2(2)) % Not consecutive (not in one segment)
    if norm(Y(Yedge(Y2(1)),:)-Yi)<norm(Y(Yedge(Y2(2)),:)-Yi)
        Ys(2)=Yedge(Y2(1));
    else
        Ys(2)=Yedge(Y2(2));
    end
end
% Project on segment
Y12=Y(Ys(2),:)-Y(Ys(1),:);
Yi1=Yi(1:2)-Y(Ys(1),1:2);
s=Yi1*Y12(1:2)'/norm(Y12);
Ynew=Y(Ys(1),:)+s*Y12/norm(Y12);
% clf;plot(Y(Yedge,1),Y(Yedge,2));axis equal;hold on;plot(Yi(1,1),Yi(1,2),'o');plot(Ynew(1,1),Ynew(1,2),'x');hold off
end 
