function[Cell]=xBottvbot(Cb,CellOld,xInternalB,nTriNodeB,TriNodeB,Y,Ct,xInternalT,nTriNodeT,TriNodeT)
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

if isempty(CellOld)
    ncell=length(xInternalB);   % number of Cells (it should be the same as the top)
else 
    ncell=length(CellOld);
end 
Cell=cell(ncell,1);         % Initialize 
nBTriangles=size(Cb,1);
nBnodes=size(TriNodeB,1);   % number of bottom nodes 

for i=1:ncell               % loop over cells  
    %% bottom
    if isempty(CellOld)
        nB=xInternalB(i);   % number of Cells (it should be the same as the top)
    else 
        nB=CellOld{i}.xBott;
    end     
    Cell{i}.xBott=nB;
    
    nvertB=nTriNodeB(nB);           % Number of vertex/tri connected to nB/i
    vertB=TriNodeB(nB,1:nvertB);    % List of vertices round cell-center i (they are to be orderd)
    % order the list of found verstices (vf ordered) (vert List) 
    vfB=zeros(size(vertB));          % Ordered list found vertices around cell-center i
    vfB(1)=vertB(1);                 % Starting vertex
    xbot=zeros(length(vertB),2);     % List Nodel bars 
    xbot(:,2)=nB; % all  bars connected to center
    vertB(1)=0; % To avoid checking against itself (remove first vertex)
    
    for v=1:nvertB-1 % loop over vertex/tri connected to nB/i
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
            error('Triangle next to the triangle [%i %i %i] was not found \nMight be filtered, try to chenge the value of (Set.RemodelTolF)',Cb(vfB(v),1),Cb(vfB(v),2),Cb(vfB(v),3));
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
            xbot(v+1,1)=nC;
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
    if nargin>6     % Different Top/Bottom
        if isempty(CellOld)
           nT=xInternalT(i);   % number of Cells (it should be the same as the top)
        else 
            nT=CellOld{i}.xTopp;
        end  
        nt=nT-nBnodes;                         % going from globl indices to partical 
        Cell{i}.xTopp=nT;                      % Top Node
        nvertT=nTriNodeT(nt);                  % Number of vertex/tri connected to nT/i
        vertT=TriNodeT(nT-nBnodes,1:nvertT);   % List of vertices round cell-center i (they are to be orderd)
        vfT=zeros(size(vertT));       % Ordered list found vertices around cell-center i
        vfT(1)=vertT(1);              % Starting vertex
        xTop=zeros(length(vertT),2); % List Nodel bars 
        xTop(:,2)=nT; % all  bars connected to center
        vertT(1)=0; % To avoid checking against itself (remove first vertex)

        for v=1:nvertT-1 % loop over vertex/tri connected to nB/i
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
                % the next vertex is not found
                error('Triangle nexto to triangle %i %i %i not found',Ct(vf(v),:));
            else
                % the next vertex is found 
                vfT(v+1)=vertT(vc-1); % put the next vertex in order 
                vertT(vc-1)=0;        % mark as (an ordered vertex)

                nC=ismemb(1);         % take the first node of the shared nodes    
                if (nC==nT)           % if it the same as the cell node nB/i
                    nC=ismemb(2);     % take the other 
                end

                if v==1
                    xTop(1,1)=Ct(vft(v),~ismember(Ct(vft(v),:),ismemb));
                end
                xTop(v+1,1)=nC;
            end
        end
        % Check orientation
        v1=Y(vfT(2),:)-Y(vfT(1),:);
        v2=Y(vfT(3),:)-Y(vfT(2),:);
        if cross(v1,v2)*[0 0 1]'<0
            vfT=vfT(end:-1:1);
            xTop(:,1)=xTop(end:-1:1,1);
        end
        Cell{i}.vtop=vfT';
        Cell{i}.xTop=xTop;
    else    % same Top/Bottom 
        Cell{i}.xTopp=nB+nBnodes;
        Cell{i}.vtop=vfB'+nBTriangles;
        Cell{i}.xTop=xbot+nBnodes;
    end 

    if length(CellOld)>=i
        if isfield(CellOld{i},'Vol0')
            Cell{i}.Vol0=CellOld{i}.Vol0;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%