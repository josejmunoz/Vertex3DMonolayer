function C=MatrixC(Cb,nTriB,TriB,Xb,Xbt,Xt,Ct,nTriT,TriT)
%          C=MatrixC(Cb,Ct,nTriNodeB,nTriNodeT,TriNodeB,TriNodeT,Xb,Xbt,Xt);
% Another use        C=MatrixC(Cb,nTriNodeB,TriNodeB,Xb,Xbt,Xt);
% Function that builds the connectivity of the bar eleemnts between nodes 
%  it works by doing three nested loop 
%  1) the first is over the nodes Xb/Xt(i)
%  2) the secound is over the tringules(t) connected to nodes(i)
%  3) the third is over the nodes(s) on triangle(t)  
    

%% Output
% Returns connectivity matrix C of nodes
% C(i,:)=nodes forming nodal element i

%% Initilize 
nnodes=length(nTriB);     % num of nodes (should be the same B/T)
nTrianglesB=size(Cb,1);    % num of Tringules Bottom
IncludedB=zeros(nnodes,1); % bottom Global indicator 
IncludedT=zeros(nnodes,1); % Top    Global indicator
C=zeros(4*nTrianglesB*3,2); % ???
k=0;

if nargin == 6 
    %% Same Top/Bottom
    for i=1:length(Xb)
        xb=Xb(i);
        xt=Xt(i);
        C(k+1,:)=[xb xt]; % vertical
        k=k+1;
        iIncluded=zeros(nnodes,1);
        iIncluded(xb)=1;
        for t=1:nTriB(xb) % Loop on connected triangles
            iTri=TriB(xb,t);
            for s=1:3 % Loop on nodes of connected triangles
                xbb=Cb(iTri,s); % Bottom node
                if ~IncludedB(xbb) && ~iIncluded(xbb)
                    xtt=Xbt(xbb); % Top node
                    C(k+1,:)=[xb xbb]; % horizontal bottom
                    C(k+2,:)=[xt xtt]; % horizontal top
                    C(k+3,:)=[xt xbb]; % diagonal
                    k=k+3;
                    iIncluded(xbb)=1;
                end
            end
        end
        IncludedB(xb)=1; % All nodal elements from xb and xt have been completed
    end
    
else
    %% Different Top/Botoom
    for i=1:length(Xb)                  
        xb=Xb(i);          % take top node 
        xt=Xt(i);          % take bottom node
        C(k+1,:)=[xb xt]; % vertical bar (CenterB-CenterT) 
        k=k+1;
        % Bottom ----------------------------------------------------------
        iIncluded=zeros(nnodes,1); %% Local indicator    
        iIncluded(xb)=1;
        for t=1:nTriB(xb)       % Loop over connected triangles (t)
            iTri=TriB(xb,t);    % triangle number (t)
            for s=1:3            % Loop on nodes (s) of connected triangles (t)
                xbb=Cb(iTri,s);  % Bottom node
                if ~IncludedB(xbb) && ~iIncluded(xbb)   % if the node is not included 
                    C(k+1,:)=[xb xbb]; % horizontal bottom
                    k=k+1;
                    iIncluded(xbb)=1;
                end
            end
        end
        % Top -------------------------------------------------------------
        iIncluded=zeros(nnodes,1);    
        iIncluded(xt-nnodes)=1;
        for t=1:nTriT(xt-nnodes) % Loop over connected triangles
            iTri=TriT(xt-nnodes,t)-nTrianglesB;   % triangle numver 
            for s=1:3 % Loop on nodes of connected triangles
                xtt=Ct(iTri,s); % Bottom node
                if ~IncludedT(xtt-nnodes) && ~iIncluded(xtt-nnodes)   % if the node is not included 
                    C(k+1,:)=[xt xtt]; % horizontal top
                    C(k+2,:)=[xt xtt-nnodes]; % diagonal
                    k=k+2;
                    iIncluded(xtt-nnodes)=1;
                end
            end
        end

        IncludedB(xb)=1; % All nodal elements from xb and xt have been completed
        IncludedT(xt-nnodes)=1; % All nodal elements from xb and xt have been completed 
    end 
end 
C(k+1:end,:)=[];



end



% older version 
% function C=MatrixC(Cb,nTri,Tri,Xb,Xbt,Xt)
% % Returns connectivity matrix C of nodes
% % C(i,:)=nodes forming nodal element i
% nnodes=length(nTri);
% nTriangles=size(Cb,1);
% Included=zeros(nnodes,1);
% C=zeros(4*nTriangles*3,2);
% k=0;
% for i=1:length(Xb)
%     xb=Xb(i);
%     xt=Xt(i);
%     C(k+1,:)=[xb xt]; % vertical
%     k=k+1;
%     iIncluded=zeros(nnodes,1);
%     iIncluded(xb)=1;
%     for t=1:nTri(xb) % Loop on connected triangles
%         iTri=Tri(xb,t);
%         for s=1:3 % Loop on nodes of connected triangles
%             xbb=Cb(iTri,s); % Bottom node
%             if ~Included(xbb) && ~iIncluded(xbb)
%                 xtt=Xbt(xbb); % Top node
%                 C(k+1,:)=[xb xbb]; % horizontal bottom
%                 C(k+2,:)=[xt xtt]; % horizontal top
%                 C(k+3,:)=[xt xbb]; % diagonal
%                 k=k+3;
%                 iIncluded(xbb)=1;
%             end
%         end
%     end
%     Included(xb)=1; % All nodal elements from xb and xt have been completed
% end
% C(k+1:end,:)=[];
% end





