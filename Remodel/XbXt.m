function [Set,Tb,Tt,X,X0,Xb,Xt]=XbXt(Cell,Set,T,X,X0)
% Returns list of
% Xb = bottom nodes
% Xt = top nodes
Tb=T;
Tt=T;
ncell=length(Cell);
Xt=zeros(ncell,1);
Xb=zeros(ncell,1);
for i=1:ncell
    Xt(i)=Cell{i}.xTopp;
    Xb(i)=Cell{i}.xBott;
end
% Add boundary cells
nnodes=size(X,1);
XBt=zeros(nnodes,1);
XBb=zeros(nnodes,1);
kt=0;
kb=0;
IncludedT=zeros(nnodes,1);
IncludedB=zeros(nnodes,1);
for i=1:ncell
    for v=1:size(Cell{i}.xTop,1)
        for s=1:2
            if ~IncludedT(Cell{i}.xTop(v,s))
                if min(abs(Cell{i}.xTop(v,s)-Xt))>0 % Not in counted top nodes
                    kt=kt+1;
                    XBt(kt)=Cell{i}.xTop(v,s);
                    IncludedT(Cell{i}.xTop(v,s))=1;
                end
            end
        end                                           %%%   Added MALIK
    end                                               %%%   Added MALIK
    for v=1:size(Cell{i}.xBot,1)                      %%%   Added MALIK
        for s=1:2                                     %%%   Added MALIK
            if ~IncludedB(Cell{i}.xBot(v,s))
                if min(abs(Cell{i}.xBot(v,s)-Xb))>0 % Not in counted bottom nodes
                    kb=kb+1;
                    XBb(kb)=Cell{i}.xBot(v,s);
                    IncludedB(Cell{i}.xBot(v,s))=1;
                end
            end
        end
    end
end
% Add potential nodes at the corner (not connected to cells) that may be missing
TT=T;
TT(any(TT==0,2),:)=[];
Sides=[2 3
    1 3
    1 2];
for e=1:size(TT,1)
    for n=1:3
        node=TT(e,n);
        if min(IncludedT(node))==0 && min(IncludedB(node))==0 && min(abs(node-Xb))>0 && min(abs(node-Xt))>0
            na=TT(e,Sides(n)); % Search connecting node
            if IncludedB(na)>0 || min(abs(na-Xb))==0 % Is a bottom node
                kb=kb+1;
                XBb(kb)=node;
                IncludedB(node)=1;
            else
                kt=kt+1;
                XBt(kt)=node;
                IncludedT(node)=1;
            end
        end
    end
end
XBt(kb+1:end)=[];
XBb(kt+1:end)=[];
Xt=[Xt; XBt];
Xb=[Xb; XBb]; % Some nodes
%-------------------------------Malik Added----------------------------
% it was found that after remodelimg the bottom and top node IDs are not consistent
% such that (IXt(i) ~= IXb(i)+max(IXb)) which is essential for discretization,
% temporary solution. 
  Xt=Xb+max(Xb);
%----------------------------------------------------------------------

% Build Tt and Tt
for e=1:size(T,1)
    if min(abs(T(e,1)-Xb))==0 && min(abs(T(e,2)-Xb))==0 && min(abs(T(e,3)-Xb))==0
        Tt(e,1)=0;
    else
        Tb(e,1)=0;
    end
end
Tt(Tt(:,1)==0,:)=[];
Tb(Tb(:,1)==0,:)=[];
nnodes=size(Xb,1);
NodesI(Xb)=1:nnodes;
Tb=[NodesI(Tb(:,1))' NodesI(Tb(:,2))' NodesI(Tb(:,3))'];
NodesI=1:2*nnodes;
NodesI(Xt)=nnodes+1:2*nnodes;
Tt=[NodesI(Tt(:,1))' NodesI(Tt(:,2))' NodesI(Tt(:,3))'];
end