function [T,iter] = DelaunayR( X,T,delta,xExternal,RemodelExternal,RemodelEntangled)
%% DELAUNAYR : Regulated Dalunay
% Computes Delaunay triangulation from previous coordinates X and
% connectivity T.
% Minimal usage:
%  T=DelaunayR(X)
% INPUT
% X=nodal coordinates. X(in,:)=[x y z] coord of node in.
% T=initial connectivity from which aspect ratios improved.
%   If not present, uses Matlab delaunay connectivity.
% function.
% delta= Relative improvement in aspect ratio when building Delaunay.
%       delta=0, Standard Delaunay.
%       delta>0, Aspects ratios >= (delta+1)*AspectRatioDelaunay are allowed
% OUTPUT
% T = resulting connectivity
% iter = number of iterations (loops on all mesh for improving flips in connecticity)
% LOCAL:
% Sides(e,s) Element to which element i is connected in side s. Side s is opposite to node s.
% SidesElem(i,:) : Global side number per element. Numebering: first Internal, then External
% SidesIntExt(i,:) : List of Internal Sides [elem1 elem2 side1 side2], external sides [elem1 0     side1 0    ]
%
% Generate list of sides
if nargin<5
    RemodelExternal=false;
end
if nargin<4
    xExternal=[];
end
Standard=false;
if ~exist('T','var')
    T=delaunay(X);
    Standard=true;
end
if ~exist('delta','var')
    delta=0;
end
dim=size(X,2);
if dim==2
    sx=[ 1 2 3 1 2]; % extended indexing of sides
else
    sx=[ 1 2 3 4 1 2 3]; % extended indexing of sides
end
SidesL=Lengths(T,X);
% Build list of connected sides
[Sides,SidesElem,SidesIntExt]=BuildSides(T,sx);
% Correct Overlappings: Ensure trianguilation has no non-convex qudrilaterals
if ~Standard && RemodelEntangled
    [Sides,SidesL,SidesElem,SidesIntExt,T]=RemeshNonConvexQ(Sides,SidesL,SidesElem,SidesIntExt,sx,T,X,xExternal,RemodelExternal);
end
% Correct optimality. Keep swapping quadrilaterals while improvement in optimatlity
% of aspect ratios,
T=SetOptimal(delta,Sides,SidesL,SidesElem,SidesIntExt,sx,T,X,RemodelExternal);
end
%%
function [Sides,SidesElem,SidesIntExt]=BuildSides(T,sx)
% Bulds set of sides per element
% INPUT
% T(i,:) : nodes connected to element T
% OUTPUT:
% Sides(e,s) Element to which element i is connected in side s. Side s is opposite to node s.
% SidesElem(i,:) : Global side number per element. Numebering: first Internal, then External
% SidesIntExt(i,:) : List of Internal Sides [elem1 elem2 side1 side2], external sides [elem1 0     side1 0    ]
nelem=size(T,1);
dim=size(T,2)-1;
ns=size(T,2); % Number of sides (faces). 2D=3, 3D=4
Sides=zeros(nelem,ns);
SidesElem=zeros(nelem,ns);
SidesInt=zeros(nelem*ns,4);
SidesExt=zeros(nelem*ns,4);
Is=0; % internal
Es=0; % external
for e=1:nelem
    for s=1:ns
        nfound=Sides(e,s)==0;
        f=1;
        while nfound && f<=nelem
            if f==e
                f=f+1;
                continue;
            end
            for t=1:ns % Loop on adjacent element
                if max(abs(sort(T(e,sx(s+1:s+dim)))-sort(T(f,sx(t+1:t+dim)))))==0
                    Sides(e,s)=f;
                    Sides(f,t)=e;
                    Is=Is+1;
                    SidesInt(Is,:)=[e f s t];
                    SidesElem(e,s)=Is;
                    SidesElem(f,t)=Is;
                    nfound=false;
                    break;
                end
            end
            f=f+1;
        end
        if nfound % Mark as external
            Es=Es+1;
            SidesExt(Es,1:3)=[e 0 s];
            SidesElem(e,s)=-Es;
        end
    end
end
SidesIntExt=[SidesInt(1:Is,:)
    SidesExt(1:Es,:)];
SidesElem(SidesElem<0)=-SidesElem(SidesElem<0)+Is; % Add external at the end
end
%%
function SidesL=Lengths(T,X)
% Builds lengths of each side
% INPUT
% T(i,:)  nodes forming triangle i
% X(i,:)  coordintes of node i
% OUTPUT
% SidesL(i,:) lengtsh of sides of triangle i
nelem=size(T,1);
dim=size(X,2);
ne=3*(dim-1); % Number of segments in 3D=6
SidesL=zeros(nelem,ne); % List of length per segment
if dim==2
    S=[ 2 3
        3 1
        1 2];
else
    S=[1 2
        2 3
        3 1
        1 4
        2 4
        3 4];
end
% Build list of sides / elem
for e=1:nelem
    for s=1:ne
        SidesL(e,s)=norm(X(T(e,S(s,1)),:)-X(T(e,S(s,2)),:));
    end
end
end
%%
function w=AspectRatioL(l)
%%
% l= vector of lengths
% w=aspect ratio as circumradius/inradius
dim=length(l)-1;
if dim==2
    w=0.5*(l(1)+l(2)-l(3))*(l(1)+l(3)-l(2))*(l(3)+l(2)-l(1))/l(1)/l(2)/l(3);
end
w=1/w;
end
%%
function w=AspectRatioP(Xe)
%%
% l= vector of lengths
% w=aspect ratio as Perimeter/Area (in 3D, Area/Volume)
dim=size(Xe,2);
if dim==2
    y12=norm(Xe(2,:)-Xe(1,:));
    y13=norm(Xe(3,:)-Xe(1,:));
    y23=norm(Xe(3,:)-Xe(2,:));
else
    y12=Xe(3);
    y13=Xe(2);
    y23=Xe(1);
end
P=y12+y13+y23;
%     Area=norm(cross(y12,y13))/2;
Area=sqrt(P*(P/2-y12)*(P/2-y13)*(P/2-y23)/2);
w=P^2/Area;
end
%%
function f=GetElemNum(e,SidesElem)
% GetSide number of adjacent element f.
ns=size(SidesElem,1);
f=zeros(ns,1);
for i=1:ns
    if e==SidesElem(i,1)
        f(i)=SidesElem(i,2);
    elseif e==SidesElem(i,2)
        f(i)=SidesElem(i,1);
    else
        error('Adjacent element %i not found in list %i %i',e,SidesElem(i,:));
    end
end
end
%%
function s=GetSideNum(e,f,SidesIntExt)
% GetSide number of adjacent element f.
if e==SidesIntExt(1) && f==SidesIntExt(2)
    s=SidesIntExt(4);
elseif e==SidesIntExt(2) && f==SidesIntExt(1)
    s=SidesIntExt(3);
else
    error('Adjacent elements %i and %i  not found in %i %i %i %i',e,f,SidesIntExt);
end
end
%%
function  [Sides,SidesL,SidesElem,SidesIntExt,T]=RemeshNonConvexQ(Sides,SidesL,SidesElem,SidesIntExt,sx,T,X,xExternal,RemodelExternal)
% Checks that no overlapping of triangles exist in triangulation T. Detects
% non-convex quadrilaterals formed by 2 triangles.
% INPUT:
% T(i,:) nodes forming triangle i
% X(i,:) coordinates of node i
% SidesElem(i,:) : Global side number of element i. 1st Side is opposite to 1st node in element connectivity.
% SidesIntExt(i,:) : List of Internal Sides [elem1 elem2 side1 side2], external sides [elem1 0     side1 0    ]
% OUTPUT:
% T(i,:) nodes forming triangles i
if nargin<9
    RemodelExternal=false;
end
Node1=[1 2 4
    1 2 3
    1 3 4];
Node2=[3 4 2];
Xq=zeros(4,size(X,2));
Tq=1:4;
nchanges=1;
while nchanges>0
    nchanges=0;
    nsides=size(SidesIntExt,1);
    for e=1:nsides % loop on internal sides
        e1=SidesIntExt(e,1);
        e2=SidesIntExt(e,2);
        s1=SidesIntExt(e,3);
        s2=SidesIntExt(e,4);
        if min(e1,e2)==0 || s1<0 % Skip if external edge or adjacent to procesed polytop
            continue;%break;
        elseif (min(Sides(e1,:))==0 || min(Sides(e2,:))==0) && ~RemodelExternal && ...
               (sum(ismember(xExternal,T(e1,:)))>0 || sum(ismember(xExternal,T(e2,:)))>0)
            continue;%break;
        end
        Xq(Node1(s1,:),:)=X(T(e1,:),:);
        Xq(Node2(s1),:)=X(T(e2,s2),:);
        [C,N,Knots]=Convex((1:4),Xq); % N=nodes adjacent to non-convex node in local numbering.
        Tq([Node1(s1,:) Node2(s1)])=[T(e1,:) T(e2,s2)];
        N=Tq(N); % Nodes in global numbering
        % Swap if non-convex, no knots, and 2 adjacent nodes to non-convex node form common edge
        if ~C && ~Knots && sum(ismember(T(e1,:),N([1 3])))==2 && sum(ismember(T(e2,:),N([1 3])))==2
            nchanges=nchanges+1;
            Tnc=find(sum(ismember(T,N(2)),2)>0); % Triangles connected to non-convex node
            Tnew=zeros(length(Tnc)+1,3); % there is one element not connected to N(2)
            Tnew(1,:)=N(2:4);
            Tnew(2,:)=[N(2) N(4) N(1)];
            k=3;
            for i=1:length(Tnc) % Loop on elements
                if Tnc(i)==e1 ||Tnc(i)==e2
                    continue;
                end
                OtherNodes=T(Tnc(i),T(Tnc(i),:)~=N(2));
                Ti=Tnc(sum(ismember(T(Tnc,:),[N(2) OtherNodes(1:2)]),2)==3);
                if isempty(Ti)
                    [n1,n2,n3]=OrientNodes(N(2),OtherNodes(1), OtherNodes(2),X);
                else
                    n1=T(Ti,1);
                    n2=T(Ti,2);
                    n3=T(Ti,3);
                end
                Tnew(k,:)=[n1 n2 n3];
                k=k+1;
            end
            Te=[e1 e2];
            Tnc=[Tnc' Te(~ismember([e1 e2],Tnc))];
            T(Tnc,:)=Tnew;
            S=SidesElem(Tnc,:);
            S1=S(S(:,1)>0,1);
            S2=S(S(:,2)>0,2);
            S3=S(S(:,3)>0,3);
            SidesIntExt(S1,3)=-1; % Mark as processed sides
            SidesIntExt(S2,3)=-1; % Mark as processed sides
            SidesIntExt(S3,3)=-1; % Mark as processed sides
        end
    end
    if nchanges>0
        SidesL=Lengths(T,X);
        [Sides,SidesElem,SidesIntExt]=BuildSides(T,sx);
    end
end
end
%%
function  [Sides,SidesL,SidesElem,SidesIntExt,T]=RemeshNonConvexL(Sides,SidesL,SidesElem,SidesIntExt,sx,T,X)
% Remeshes locally when triangles with counterclock wise direction exists
% INPUT:
% T(i,:) nodes forming triangle i
% X(i,:) coordinates of node i
% SidesElem(i,:) : Global side number of element i. 1st Side is opposite to 1st node in element connectivity.
% SidesIntExt(i,:) : List of Internal Sides [elem1 elem2 side1 side2], external sides [elem1 0     side1 0    ]
% OUTPUT:
% T(i,:) nodes forming triangles i
Node1=[1 2 4
    1 2 3
    1 3 4];
Node2=[3 4 2];
nchanges=1;
while nchanges>0
    nchanges=0;
    for e=1:nelem 
        Te=T(e,:);
        [~,~,~,Orient]=OrientNodes(Te(1),Te(2),Te(3),X);
        if Orient
            nchanges=nchanges+1;
            Tnc=find(sum(ismember(T,Te),2)>0);
            TLocal=T(Tnc,:);
            SidesLocal=Sides(Tnc,:);
            [~,xExternal]=ExternalInternal(TLocal,SidesLocal);
            xExternal=SortExternal(xExternal,TLocal,SidesLocal);
            TLocalNew=delaunay(X(Tnc,1),X(Tnc,2));
            nelemNew=size(TLocaNew,1);
            for en=1:nelemNew
                Ten=TLocalNew(e,:);
                if IsExternal(xExternal,Ten,X)
                    TLocalNew(e,:)=0;
                end
            end
            TLocalNew(TLocalNew(:,1)<0,:)=[];
            nc=length(Tnc);
            nLocal=size(TLocalNew,1);
            if nLocal>nc
                T(Tnc,:)=TLocalNew(1:nc,:);
                T(end+1:end+nLocal-nc)=TLocalNew(nc+1:end,:);
            elseif size(TLocalNew,1)<nc
                T(Tnc(1:nLocal),:)=TLocalNew;
                T(Tnc(nLocal+1:nc),:)=[];
            else
                T(Tnc,:)=TLocalNew;
            end
        end
    end
    if nchanges>0
        SidesL=Lengths(T,X);
        [Sides,SidesElem,SidesIntExt]=BuildSides(T,sx);
    end
end
end
%%
function Is=IsExternal(xExternal,T,X)
% Determines if triangles in T are external. UNFINISHED.
% INPUT
% xExternal = ordered list of nodes forming boundary in clockwise direction
% T(i,:)    = nodes forming triangle i
% X(i,:)    = coordinates of node i
% OUTPUT
% Is(i)  = true if triangle i is external to boundary set by xExternal
nelem=size(T,1);
for e=1:nelem
    Te=T(e,:);
    Ti=ismember(xExternal,Te);
    if sum(Ti)==0
    elseif sum(Ti)==1
    elseif sum(Ti)==2
    else
    end
end
end
%%
function T=SetOptimal(delta,Sides,SidesL,SidesElem,SidesIntExt,sx,T,X,RemodelExternal)
% Swaps rectiangles formed by two triangles in order to get more optimal
% shapes regulated by delta
nelem=size(T,1);
ns=size(T,2); % Number of sides (faces). 2D=3, 3D=4
nchanges=1;
iter=0;
while nchanges
    iter=iter+1;
    nchanges=false;
    for e=1:nelem
        for s=1:ns
            f=Sides(e,s);
            if f>0 && (RemodelExternal || (min(Sides(e,:))>0 && min(Sides(f,:))>0)) % Not external
                % Check whether flip is favourable
                ks=SidesElem(e,s);
                if SidesIntExt(ks,1)==e % Not in 3D
                    se=SidesIntExt(ks,3); % =s
                    sf=SidesIntExt(ks,4);
                else
                    se=SidesIntExt(ks,4);
                    sf=SidesIntExt(ks,3); % =se
                end
                l1=SidesL(e,:); % Original triangles. Elem e
                l2=SidesL(f,:); % Original triangles. Elem f
                LN=norm(X(T(e,se),:)-X(T(f,sf),:)); % New length at crossing line
                L1=[SidesL(e,sx(se+1)) LN SidesL(f,sx(sf+2))]; % Flipped triangles (Not in 3D)
                L2=[SidesL(e,sx(se+2)) SidesL(f,sx(sf+1)) LN];
                T1=[T(e,sx(se+2)) T(e,se) T(f,sf)];
                T2=[T(e,sx(se+1)) T(f,sf) T(e,se)];
                w1=AspectRatioP(l1); % AspectRatioP(X(T(e,:),:)); %max(l1)/min(l1);
                w2=AspectRatioP(l2); % AspectRatioP(X(T(f,:),:)); %max(l2)/min(l2);
                W1=AspectRatioP(L1); % AspectRatioP(X(T(T1,:),:));
                W2=AspectRatioP(L2); % AspectRatioP(X(T(T2,:),:));
                [~,~,~,Ore]=OrientNodes(T(e,1),T(e,2),T(e,3),X);
                [~,~,~,Orf]=OrientNodes(T(f,1),T(f,2),T(f,3),X);
                if max(w1,w2)>max(W1,W2)*(1+delta) && Convex([T1(1:2) T2(1:2)], X) && ~Ore && ~Orf % FLIP ACCORDING TO AVOID MAXIMUM ASPECT RATIO (MODIFY IF NECESSARY)
                    % with entangled meshes, element with T1 and T2 may
                    % already exist. Avoid swapping in this case
                    if isempty(find(sum(ismember(T,T1),2)>=3,1)) && isempty(find(sum(ismember(T,T2),2)>=3,1))
                        T(e,:)=T1;
                        T(f,:)=T2;
                        % Update side lengths
                        auxe=SidesL(e,:);
                        SidesL(e,:)=[LN SidesL(f,sx(sf+2)) SidesL(e,sx(se+1))];
                        SidesL(f,:)=[LN auxe(1,sx(se+2))   SidesL(f,sx(sf+1))];
                        fi=Sides(e,sx(se+1));
                        if fi>0, SidesL(fi,Sides(fi,:)==e)=SidesL(e,3);end
                        fi=Sides(e,sx(se+2));
                        if fi>0, SidesL(fi,Sides(fi,:)==e)=SidesL(f,2);end
                        fi=Sides(f,sx(sf+1));
                        if fi>0, SidesL(fi,Sides(fi,:)==f)=SidesL(f,3);end
                        fi=Sides(f,sx(sf+2));
                        if fi>0, SidesL(fi,Sides(fi,:)==f)=SidesL(e,2);end
                        % Update List of sides
                        si=SidesElem(e,se);
                        SidesIntExt(si,:)=[e f 1 1];
                        si=SidesElem(e,sx(se+1));
                        fi=Sides(e,sx(se+1));
                        SidesIntExt(si,:)=[e fi 3 GetSideNum(e,fi,SidesIntExt(si,:))];
                        si=SidesElem(f,sx(sf+2));
                        fi=Sides(f,sx(sf+2));
                        SidesIntExt(si,:)=[e fi 2 GetSideNum(f,fi,SidesIntExt(si,:))];
                        si=SidesElem(e,sx(se+2));
                        fi=Sides(e,sx(se+2));
                        SidesIntExt(si,:)=[f fi 2 GetSideNum(e,fi,SidesIntExt(si,:))];
                        si=SidesElem(f,sx(sf+1));
                        fi=Sides(f,sx(sf+1));
                        SidesIntExt(si,:)=[f fi 3 GetSideNum(f,fi,SidesIntExt(si,:))];
                        auxe=SidesElem(e,:); % Update SidesElem (adjacent elements)
                        SidesElem(e,:)=[SidesElem(e,se) SidesElem(f,sx(sf+2)) SidesElem(e,sx(se+1))];
                        SidesElem(f,:)=[auxe(1,se)      auxe(1,sx(se+2))      SidesElem(f,sx(sf+1))];
                        Sides(e,:)=GetElemNum(e,SidesIntExt(SidesElem(e,:),1:2)); % Get adjacent elements
                        Sides(f,:)=GetElemNum(f,SidesIntExt(SidesElem(f,:),1:2));
                        fi=Sides(e,2);
                        if fi>0, Sides(fi,Sides(fi,:)==f)=e;end
                        fi=Sides(e,3);
                        if fi>0,Sides(fi,Sides(fi,:)==e)=e;end
                        fi=Sides(f,2);
                        if fi>0,Sides(fi,Sides(fi,:)==e)=f;end
                        fi=Sides(f,3);
                        if fi>0,Sides(fi,Sides(fi,:)==f)=f;end
                        nchanges=nchanges+1;
                    end
                end
            end
        end
    end
end
end
%%
function [I,Nodes,Knots]=Convex(T,X)
% Check if set of points in T form a convex polytope
% T      = list of nodes in X forming poilytop
% X(i,:) = coordinates of node i
% OUTPUT:
% I      = false if T,X forms convex polytope
% Nodes  = Nodes(1:3) nodes forming nont-convex side. N(2) is the centered one.
%          Nodes(4:end) the reminaing nodes
% knots  = true: at least one know exist (number of turns =0)
n=length(T);
T(end+1:end+2)=[T(1) T(2)]; % Force connected loop
N=[1:n 1 2];
Nodes=1:n;
I=true;
alpha=0;
for i=1:n
    x1=X(T(i),:);
    x2=X(T(i+1),:);
    x3=X(T(i+2),:);
    x12=[x2-x1 0];
    x23=[x3-x2 0];
    v=cross(x12,x23);
    if x12*x23'>0
        alpha=alpha+sign(v(3))*asin(norm(v)/norm(x12)/norm(x23));
    else
        alpha=alpha+sign(v(3))*(pi-asin(norm(v)/norm(x12)/norm(x23)));
    end
    if v(3)<0
        I=false;
        Nodes=T(N([i i+1 i+2])); % Three nodes forming non convex side
        if length(T(~ismember(T(1:n),Nodes)))~=4-n+1
            error('Could not assign 4th node to qudrilateral while checking convexity.');
        else
            Nodes(4:n)=T(~ismember(T(1:n),Nodes));
        end
    end
end
if abs(alpha-2*pi)<1e4*eps
    Knots=false;
else
    Knots=true;
end
end
