function [L,Ln,L0,Stress]=RemodelLS(C,COld,Cv,CvOld,LnOld,LOld,StressOld,X,X0,Y,VNew2Old,L0Old)
% Sets to current length the  new rest length of new elements.
[COld,IOld]=sortrows(sort(COld,2));
[CvOld,IOldv]=sortrows(sort(CvOld,2));
nelen=size(C,1);
nelev=size(Cv,1);
L.V.L=zeros(nelev,1);
L.D.L=zeros(nelen,1);
L0.V.L=zeros(nelev,1);
L0.D.L=zeros(nelen,1);
Stress.V.s=zeros(nelev,1);
Stress.D.s=zeros(nelen,1);
for i=1:nelen
    nodes=C(i,:);
    e=IncludedD(nodes,COld);
    if e>0
        L.D.L(i)=LOld.D.L(IOld(e));
        Ln.D.L(i)=LnOld.D.L(IOld(e));
        Stress.D.s(i)=StressOld.D.s(IOld(e));
    else
        L.D.L(i)=norm(X(nodes(1),:)-X(nodes(2),:));
        Ln.D.L(i)=L.D.L(i);
        Stress.D.s(i)=0;
    end
    L0.D.L(i)=norm(X0(nodes(1),:)-X0(nodes(2),:));
end
for i=1:nelev
%     nodes=Cv(i,:);                                       %% Malik Commented 
%     e=IncludedD(nodes,CvOld);                            %% Malik Commented 
%     if e>0 && CheckL0                                    %% Malik Commented 
%-----------------------------------MALIK ADDED----------------------------
     vert=Cv(i,:);
     vertOld=VNew2Old(vert);
     e=IncludedV(vertOld,CvOld);
%      e=0;
    if e>0
        L.V.L(i)=LOld.V.L(IOldv(e));
        Ln.V.L(i)=LnOld.V.L(IOldv(e));
        L0.V.L(i)=L0Old.V.L(IOldv(e)); 
        Stress.V.L(i)=StressOld.V.s(IOldv(e));
    else
        L.V.L(i)=norm(Y(vert(1),:)-Y(vert(2),:));
        Ln.V.L(i)=L.V.L(i);
        L0.V.L(i)=L.V.L(i);
        Stress.V.s(i)=0;
    end
%     L0.V.L(i)=L.V.L(i); % Assume rest length = current length
%     Ln.V.L(i)=L.V.L(i); % Assume rest length = current length   
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function     e=IncludedD(nodes,C)
% Checks if pair of elements in nodes is included in matrix C along a row.  
% Assumes C is sorted along row and then along columns
% e= row where list nodes is found. If =0, not found.
nodes=sort(nodes);
nelem=size(C,1);
i=1;
e=0;
while i<=nelem
    if nodes==C(i,:)
        e=i;
        break;
    end
    i=i+1;
end
end

function  e=IncludedV(vertOld,CvOld)
% % Checks if pair vertOld in included in CvOld. 
e=0;
if isnan(vertOld(1)) || isnan(vertOld(2))
    % Old Vertices are unknown  
    return
end 
nelem=size(CvOld,1);
i=1;
while i<=nelem
    if any(vertOld(1)==CvOld(i,:)) && any(vertOld(2)==CvOld(i,:))
        
        e=i;
        break;
    end
    i=i+1;
end
end



