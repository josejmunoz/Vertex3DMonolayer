function sExt=SortExternal(xExt,C,Edges)
% Orders external nodes in the same sense as appearing in triangulation
% INPUT:
% xExt       : list of external nodes
% C(i,:)     : nodes forming triangle i
% Edges(i,s) : element connected to side s of triangle i. Side s is
% opposite to node in position s of T(i,:).
% OUTPUT:
% xExt       : ordered list (clockwise) of external nodes
sExt=0*xExt;
xExtE=sExt; % List of nodes that have been countd
sExt(1)=xExt(1);
% 
Next=[2 3 1];
Side=[3 1 2];
% Build external elements
eExt=sum(ismember(Edges,0),2)>=1;
Cext=C(eExt,:); % External triangles
EdgesExt=Edges(eExt,:);
n=xExt(1);
xExtE(xExt==n)=1;
i=1;
Finish=false;
while ~Finish
    Found=false;
    Ci=find(sum(ismember(Cext,n),2)>0); % External triangles that contain n
    for ex=1:length(Ci) % Loop on elements that contain external node n
        for s=1:3 % Loop on nodes of triangle
            if Cext(Ci(ex),s)==n
                n1=Cext(Ci(ex),Next(s));
                if sum(ismember(sExt,[n n1]))>=2 
                    ni=find(sExt==n1);
                    if ~isempty(ni)
                        if ni(1)>1
                            if sExt(ni(1)-1)==n
                                continue% Skip if this sequence has already been found (this may happen when domain connected through one single node)
                            end
                        end
                    end
                end
                if any(ismember(xExt,n1))&& EdgesExt(Ci(ex),Side(s))==0 % Found external node linked to n through an external edge
                    if sum(xExtE==0)==0
                        Finish=true;
                        break
                    else
                        i=i+1;
                        sExt(i)=n1;
                        xExtE(xExt==n1)=i;
                        n=n1;
                        Found=true;
                        break;
                    end
                end
            end
        end
        if Found
            break;
        end
    end
end
end
