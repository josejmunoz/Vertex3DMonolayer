function[Cell]=xBottvbotDA(Cb,CellOld,nTriNodeB,TriNodeB,Y,Ct,nTriNodeT,TriNodeT,Ablated)
% Another Use  [Cell]=xBottvbot(Cb,CellOld,xInternalB,nTriNodeB,TriNodeB,Y)  
% INPUT:
    % nTri = number of triangles for each node
    % Tri = list of triangles for each nodes
    % Builds structure Cell with:
    % Cell{i}.xBott=bottom node of cell i
    % Cell{i}.xTopp=top node of cell i
    % Cell{i}.vbot=list of vertices at bottom of cell i
    % Cell{i}.vtop=list of vertices at top of cell i
    % Cell{i}.xBot=list of nodal bar elements connected to cell at bottom layer
    % Cell{i}.xTop=list of nodal bar elements connected to cell at top layer


ncell=length(CellOld);   % number of Cells (it should be the same as the top) 
Cell=cell(ncell,1);         % Initialize 
nBTriangles=size(Cb,1);
nBnodes=size(TriNodeB,1);   % number of bottom nodes 

for i=1:ncell               % loop over cells  
    %% bottom
    nB=CellOld{i}.xBott;
    Cell{i}.xBott=nB;                % bottom Node
    nvertB=nTriNodeB(nB);           % Number of vertex/tri connected to nB/i
    nnodesB=nvertB;
    if Ablated.Exist 
        if any(any(ismember(Ablated.VRingBot(:,[2 3]),nB),2))
            nnodesB=nnodesB+1;
            nvertB=nvertB+2;
        end 
    end 
    vertB=TriNodeB(nB,1:nvertB);    % List of vertices round cell-center i (they are to be orderd)
    % order the list of found verstices (vf ordered) (vert List) 
    vfB=zeros(size(vertB));          % Ordered list found vertices around cell-center i
    vfB(1)=vertB(1);                 % Starting vertex
    xbot=zeros(nnodesB,2);     % List Nodel bars 
    xbot(:,2)=nB; % all  bars connected to center
    vertB(1)=0; % To avoid checking against itself (remove first vertex)
    n=1;
    v=1;
    while (v<=nvertB-1) % loop over vertex/tri connected to nB/i
        FoundNext=0;
        vc=1; % Current triangle inspected  (counter)
        % while the next vertics is not found && and the inspected vc < nvert-1
        while ~FoundNext && vc<=nvertB 
            if vertB(vc)>0  % 
                % List of common nodes
                ismemb=Cb(vertB(vc),ismember(Cb(vertB(vc),:),Cb(vfB(v),:))); 
                FoundNext=length(ismemb)==2; 
            end
            vc=vc+1;
        end
        if ~FoundNext
            % the next verstex is not found
            % find the elements connected to node nB 
            RowBar=Ablated.VRingBot(any(ismember(Ablated.VRingBot(:,[2 3]),nB),2),:);
            % find which one vertex connected to the current tringle Cb(vf(v),:))
            firstlogical=all(ismember(RowBar(:,[2 3]),Cb(vfB(v),:)),2);
            v1=RowBar(firstlogical,1);
            % the next vertex
            v2=RowBar(~firstlogical,1);

            % find the closing triangle from the other side
            Found=0;
            vc=1;
            while ~Found && vc<=nvertB
                if vertB(vc)>0
                    ismemb=Cb(vertB(vc),ismember(Cb(vertB(vc),:),RowBar(~firstlogical,[2 3]))); % List of common nodes
                    Found=length(ismemb)==2; 
                end
                vc=vc+1; 
            end
            if length(v1)~=1
                malik='debuge'
            end 
            vfB(v+1)=v1;
            vfB(v+2)=v2;
            v=v+2;

            if Found 
                 vfB(v+1)=vertB(vc-1);
                 v=v+1;
            end 
            vertB(vc-1)=0;
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
%             error('Triangle nexto to triangle %i %i %i not found',Cb(vfB(v),:));
        else
            % the next vertex is found 
            vfB(v+1)=vertB(vc-1); % put the next vertex in order 
            vertB(vc-1)=0;       % mark as (an ordered vertex)
            
            nC=ismemb(1);       % take the first node of the shared nodes    
            if (nC==nB)         % if it the same as the cell node nB/i
                nC=ismemb(2);   % take the other 
            end
            
            if v==1
                xbot(1,1)=Cb(vfB(v),~ismember(Cb(vfB(v),:),ismemb));
            end
            xbot(n+1,1)=nC;
            n=n+1;
            v=v+1;
        end
    end
    % Check orientation
    v1=Y(vfB(2),:)-Y(vfB(1),:);
    v2=Y(vfB(3),:)-Y(vfB(2),:);
    if cross(v1,v2)*[0 0 1]'<0
        vfB=vfB(end:-1:1);
        xbot(:,1)=xbot(end:-1:1,1);
    end
    Cell{i}.vbot=vfB';
    Cell{i}.xBot=xbot;
        
    %% Top
    nT=CellOld{i}.xTopp;
    nt=nT-nBnodes;                         % going from globl indices to partical 
    Cell{i}.xTopp=nT;                      % Top Node
    nvertT=nTriNodeT(nt);                  % Number of vertex/tri connected to nT/i
    nnodesT=nvertT;
    if Ablated.Exist 
        if any(any(ismember(Ablated.VRingTop(:,[2 3]),nT),2))
            nnodesT=nnodesT+1;
            nvertT=nvertT+2;
        end 
    end
    vertT=TriNodeT(nT-nBnodes,1:nvertT);   % List of vertices round cell-center i (they are to be orderd)
    vfT=zeros(size(vertT));       % Ordered list found vertices around cell-center i
    vfT(1)=vertT(1);              % Starting vertex
    xtop=zeros(nnodesT,2); % List Nodel bars 
    xtop(:,2)=nT; % all  bars connected to center
    vertT(1)=0; % To avoid checking against itself (remove first vertex)
    n=1;
    v=1;
    while (v<=nvertT-1) % loop over vertex/tri connected to nB/i
        FoundNext=0;
        vc=1; % Current triangle inspected  (counter)
        % while the next vertics is not found && and the inspected vc < nvert-1
        while ~FoundNext && vc<=nvertT 
            if vertT(vc)>0  % 
                vertt=vertT-size(Cb,1);
                vft=vfT-size(Cb,1);
                % List of common nodes
                ismemb=Ct(vertt(vc),ismember(Ct(vertt(vc),:),Ct(vft(v),:))); 
                FoundNext=length(ismemb)==2; 
            end
            vc=vc+1;
        end
        if ~FoundNext
            % the next verstex is not found
            % find the elements connected to node nB 
            RowBar=Ablated.VRingTop(any(ismember(Ablated.VRingTop(:,[2 3]),nT),2),:);
            % find which one vertex connected to the current tringle Cb(vf(v),:))
            firstlogical=all(ismember(RowBar(:,[2 3]),Ct(vfT(v)-size(Cb,1),:)),2);
            v1=RowBar(firstlogical,1);
            % the next vertex
            v2=RowBar(~firstlogical,1);

            % find the closing triangle from the other side
            Found=0;
            vc=1;
            while ~Found && vc<=nvertT
                if vertT(vc)>0
                    vertt=vertT-size(Cb,1);
                    ismemb=Ct(vertt(vc),ismember(Ct(vertt(vc),:),RowBar(~firstlogical,[2 3]))); % List of common nodes
                    Found=length(ismemb)==2; 
                end
                vc=vc+1; 
            end
            vfT(v+1)=v1;
            vfT(v+2)=v2;
            v=v+2;

            if Found 
                 vfT(v+1)=vertT(vc-1);
                 v=v+1;
            end 
            vertT(vc-1)=0;
            nC=RowBar(firstlogical,2);
            if (nC==nT)
                nC=RowBar(firstlogical,3);
            end
%             xbot=[xbot ;[0 nB]];
            xtop(n+1,1)=nC;
                n=n+1;
           if Found 
                nC=ismemb(1);
                if (nC==nT)
                    nC=ismemb(2);
                end
%                 xbot=[xbot ;[0 nB]];
                xtop(n+1,1)=nC;
                n=n+1;
           end 
            
            
            % the next verstex is not found
%             error('Triangle nexto to triangle %i %i %i not found',Ct(vf(v),:));
        else
            % the next vertex is found 
            vfT(v+1)=vertT(vc-1); % put the next vertex in order 
            vertT(vc-1)=0;        % mark as (an ordered vertex)

            nC=ismemb(1);         % take the first node of the shared nodes    
            if (nC==nT)           % if it the same as the cell node nB/i
                nC=ismemb(2);     % take the other 
            end

            if v==1
                xtop(1,1)=Ct(vft(v),~ismember(Ct(vft(v),:),ismemb));
            end
            xtop(n+1,1)=nC;
            n=n+1;
            v=v+1;
        end
     end

    % Check orientation
    v1=Y(vfT(2),:)-Y(vfT(1),:);
    v2=Y(vfT(3),:)-Y(vfT(2),:);
    if cross(v1,v2)*[0 0 1]'<0
        vfT=vfT(end:-1:1);
        xtop(:,1)=xtop(end:-1:1,1);
    end
    Cell{i}.vtop=vfT';
    Cell{i}.xTop=xtop;
 

    if length(CellOld)>=i
        if isfield(CellOld{i},'Vol0')
            Cell{i}.Vol0=CellOld{i}.Vol0;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%