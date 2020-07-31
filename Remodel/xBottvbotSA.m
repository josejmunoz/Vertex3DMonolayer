function [Cell]=xBottvbotSA(Cb,CellOld,Ablated,nTriNode,TriNode,Y)
% INPUT:
% Cb(i,:)  = nodes forming triangle i at bottom network
% nTriNode = number of triangles connected to node i
% TriNode(i,:) = triangles connected to node i
% Ablated.VRingTop(i,1) vertex number
% Ablated.VRingTop(i,2:3) nodes that vertex is connected to (if last is 0,
% triangle of vertex just has one node.
% Y(i,:)   = x,y,z coordinate of vertex
% OUTPUT:
% Builds structure Cell with:
% Cell{i}.xBott=bottom node of cell i
% Cell{i}.xTopp=top node of cell i
% Cell{i}.vbot=list of vertices at bottom of cell i
% Cell{i}.vtop=list of vertices at top of cell i
% Cell{i}.xBot=list of nodal bar elements connected to cell at bottom layer
% Cell{i}.xTop=list of nodal bar elements connected to cell at top layer
% ncell=length(xInternal);
ncell=length(CellOld);
Cell=cell(ncell,1);
nTriangles=size(Cb,1);
nodes=2*size(TriNode,1);
for i=1:ncell
%     nB=xInternal(i);
    % Find vertices around bottom and top node
    nB=CellOld{i}.xBott;
    Cell{i}.xBott=nB;
    Cell{i}.xTopp=nB+nodes/2;

    nvert=nTriNode(nB);
    nnodes=nvert;
    if Ablated.Exist 
        if any(any(ismember(Ablated.VRingBot(:,[2 3]),nB),2))
            nnodes=nnodes+1; % If cell at wound ring, nB connected to another node and 2 vertices
            nvert=nvert+2;
        end 
    end 
    vert=TriNode(nB,1:nvert); % List of vertices around cell-center nB
    vbot=zeros(1,nvert); % Ordered list found vertices around cell-center i
    % vbot=Sequential list of vertices forming ring
    vbot(1)=vert(1); % Starting vertex
    vtop=zeros(1,nvert); % Ordered list found vertices around cell-center i
    vtop(1)=vert(1)+nTriangles; % Starting vertex
    xbot=zeros(nnodes,2);
    xbot(:,2)=nB; % all  bars connected to center
    vert(1)=0; % To avoid checking against itself
    n=1;
    v=1;
    % Loop on all vertices aroundnode nB. Build vbot, vtop
    while (v<=nvert-1)
        FoundNext=0;
        vc=1; % Loop on found triangles (vertices) connected to nB
        while ~FoundNext && vc<=nvert
            if vert(vc)>0
                ismemb=Cb(vert(vc),ismember(Cb(vert(vc),:),Cb(vbot(v),:))); % List of common nodes
                FoundNext=length(ismemb)==2; % Next triangle must have 2 common nodes
            end
            vc=vc+1;
        end
        if ~FoundNext % Not found because cell at wound ring (wound not triangularised)
            % find bar elements at wound ring connected to node nB 
            RowBar=Ablated.VRingBot(any(ismember(Ablated.VRingBot(:,[2 3]),nB),2),:);
            RowBarT=Ablated.VRingTop(any(ismember(Ablated.VRingBot(:,[2 3]),nB),2),:);
            % find vertices (triangles) connected to current tringle Cb(vbot(v),:))
            firstlogical=all(ismember(RowBar(:,[2 3]),Cb(vbot(v),:)),2);
            v1=RowBar(firstlogical,1);
            v1T=RowBarT(firstlogical,1);
            % the next vertex: other triangle connected to nB and not in vbot(v)
            v2=RowBar(~firstlogical,1);
            v2T=RowBarT(~firstlogical,1);
            % find remaining triangles in contact with nB and at other side of wound ring
            Found=0;
            vc=1;
            while ~Found && vc<=nvert
                if vert(vc)>0
                    ismemb=Cb(vert(vc),ismember(Cb(vert(vc),:),RowBar(~firstlogical,[2 3]))); % List of common nodes
                    Found=length(ismemb)==2; 
                end
                vc=vc+1; 
            end
            vbot(v+1)=v1; % ERROR i=56?
            vbot(v+2)=v2;
       
            vtop(v+1)=v1T;
            vtop(v+2)=v2T;
            v=v+2;
            if Found % If still remaining triangles connected to nB (and found), complete vbot with found triangle
                 vbot(v+1)=vert(vc-1);
                 vtop(v+1)=vert(vc-1)+nTriangles;
                 v=v+1;
            end 
            vert(vc-1)=0;
            nC=RowBar(firstlogical,2);
            if (nC==nB)
                nC=RowBar(firstlogical,3);
            end
%             xbot=[xbot ;[0 nB]];
            xbot(n+1,1)=nC;
                n=n+1;
           if Found 
                nC=ismemb(1);
                if (nC==nB)
                    nC=ismemb(2);
                end
%                 xbot=[xbot ;[0 nB]];
                xbot(n+1,1)=nC;
                n=n+1;
           end 
%           error('Triangle nexto to triangle %i %i %i not found',Cb(vbot(v),:));
        else % found vertex belongs to existing triangel (not at wound)
            vbot(v+1)=vert(vc-1);
            vtop(v+1)=vert(vc-1)+nTriangles;
            vert(vc-1)=0;
            nC=ismemb(1);
            if (nC==nB)
                nC=ismemb(2);
            end
            if n==1
                xbot(1,1)=Cb(vbot(v),~ismember(Cb(vbot(v),:),ismemb));
            end
            xbot(n+1,1)=nC;
            n=n+1; % node counter
            v=v+1; % vertex counter
        end
    end
    % Check orientation
    v1=Y(vbot(2),:)-Y(vbot(1),:);
    v2=Y(vbot(3),:)-Y(vbot(2),:);
    if cross(v1,v2)*[0 0 1]'<0
        vbot=vbot(end:-1:1);
        vtop=vtop(end:-1:1);
        xbot(:,1)=xbot(end:-1:1,1);
    end
    Cell{i}.vbot=vbot';
%   Cell{i}.vtop=vbot'+size(Ablated.VRingBot,1)+2*nTriangles;
    Cell{i}.vtop=vtop';
    Cell{i}.xBot=xbot;
    Cell{i}.xTop=xbot+nodes/2;
    if length(CellOld)>=i
        if isfield(CellOld{i},'Vol0')
            Cell{i}.Vol0=CellOld{i}.Vol0;
        end
    end
end
end