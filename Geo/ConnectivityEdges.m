function Edges=ConnectivityEdges(Ct)
% Forms list of adjacent triangles
% INPUT:
% Ct(i,:) : List of nodes forming triangle i
% OUTPUT:
% Edges(i,:)= List of triangles (vertices) to which triangle i (vertex) is
% connected. If =0, is extrnal side
% nEdges(i) = number of triangles to which triangle i is connected
ntri=size(Ct,1); % = number of triangles
Edges=zeros(ntri,3);
nEdges=zeros(ntri,1);
Sides=[2 3
    3 1
    1 2];
for i=1:ntri
    lnodi=Ct(i,:);
    % Search neighbours per triangle
    if nEdges(i)<3
        for j=1:ntri
            lnodj=Ct(j,:);
            if sum(ismember(lnodi,lnodj))==2 && min(abs(Edges(i,:)-j))>0
                nEdges(i)=nEdges(i)+1;
                nEdges(j)=nEdges(j)+1;
                sj=0;
                for s=1:3
                    if sum(lnodi(Sides(s,:))==lnodj([2 1]))==2
                        sj=3;
                    elseif sum(lnodi(Sides(s,:))==lnodj([3 2]))==2
                        sj=1;
                    elseif sum(lnodi(Sides(s,:))==lnodj([1 3]))==2
                        sj=2;
                    end
                    if sj>0
                        Edges(i,s)=j;
                        Edges(j,sj)=i;
                        break;
                    end
                end
                if nEdges(i)==3
                    break;
                end
            end
        end
    end  
end
end