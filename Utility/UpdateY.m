function Y=UpdateY(T,X,N,Y,Yr,YrZ,Ymid)
% Update vertex positions
% Yr(j)=1 : vertex is relaxed
%      =0 : vertex is not relaxed (Attached to nodes)
dim=size(X,2);
nvert=size(Y,1);
if nargin==3
    YrZ=zeros(nvert,1);
    Yr=zeros(nvert,1);
    Y=zeros(nvert,dim);
end
for i=1:nvert
    if Yr(i)==0 && Ymid(i)==0
        Y(i,1:2)=0;
        if YrZ(i)==0
            Y(i,dim)=0;
        end
        x=X(T(i,:),:); % 3 nodes that interpolate vertex i
        for n=1:size(T,2)
            Y(i,1:2)=Y(i,1:2)+N(i,n)*x(n,1:2);
            if YrZ(i)==0
                Y(i,dim)=Y(i,dim)+N(i,n)*x(n,dim);
            end
        end
    end
end