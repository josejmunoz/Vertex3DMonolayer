%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,Edges]=delaunayF(T,X,Y,TolF,Delta,xExternal,RemodelExternal,RemodelEntangled)
% Apply Delaunay tiranguilation, filtering boundary triangles with aspect
% ratio w<TolF
if nargin<8
    RemodelEntangled=true;
end 
if nargin<7
   RemodelExternal=false;
end
if nargin<6
    xExternal=[];
end
if abs(Delta)<eps || isempty(T)
    C=delaunay(X,Y);
else
    C=DelaunayR([X Y],T,Delta,xExternal,RemodelExternal,RemodelEntangled);
end
X=[X Y];
Edges=ConnectivityEdges(C);
ntri=size(C,1);
Sides=[2 3
       3 1
       1 2];
% Filter external triangles
w=zeros(ntri,1);
Removed=true;
while Removed
    Removed=false;
    for i=1:ntri
        if min(Edges(i,:))==0 % External edge
            nodes=C(i,:);
            l=[norm(X(nodes(2),:)-X(nodes(1),:))...
                norm(X(nodes(3),:)-X(nodes(2),:))...
                norm(X(nodes(3),:)-X(nodes(1),:))];
            w(i)=AspectRatioL(l);
            if w(i)>TolF % Remove triangle
                for s=1:3 % Loop on sides
                    if Edges(i,s)==0
                        for j=1:2 % Loop on remaining two other edges
                            iadj=Edges(i,Sides(s,j)); % Non-external edge
                            if iadj>0
                                Edges(iadj,Edges(iadj,:)==i)=0; % Set as external
                            end
                        end
                    end
                end
                Edges(i,1)=-1;
                C(i,1)=-1;
                Removed=true;
            end
        end
    end
end
% Remove filtered triangles
Edges(Edges(:,1)==-1,:)=[];
C(C(:,1)==-1,:)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w=AspectRatioL(l)
% l= vector of lengths
% w=aspect ratio as circumradius/inradius
dim=length(l)-1;
if dim==2
    w=0.5*(l(1)+l(2)-l(3))*(l(1)+l(3)-l(2))*(l(3)+l(2)-l(1))/l(1)/l(2)/l(3);
end
w=1/w;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
