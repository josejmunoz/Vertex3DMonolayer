function [Ablated,Set]=Ablation(Ablated,Cell,C,Cv,X0,Set,Y)
% INPUT
% Cell = Cell database
% C    = connectivity
% X0   = Nodal positions
% Set  = Settings database
% OUTPUT
% Ablated:
%   Ablated.Nodal(e)=1, ablated nodal element e, =0: Non ablated nodal element
%   Ablated.Cell(e)=1, ablated cell e, =0: Non ablated cell
%   Ablated.Cv(e)=0, vertex bar element e is not ablated
%   Ablated.dofY   , list of dof corresponding to relaxed vertex nodes
%   (Startingg from 1)
% Cell{c}.vtopR(j)=j: top vertex j is relaxed
%                 =0: top vertex j is not relaxed (it is attached to nodes)
% Cell{c}.vbotR(j)=j: bottom vertex j is relaxed
%                 =0: bottom vertex j is not relaxed (it is attached to nodes)
% Set.Ablation = cell number being ablated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute total number of vertex elements:
nBar=Cell{1}.nElem;
nelem=size(C,1);
tol =min([max(X0(:,1))-min(X0(:,1)),...
    max(X0(:,2))-min(X0(:,2)),...
    max(X0(:,3))-min(X0(:,3))])/ceil(sqrt(nelem/3))*1e-10;
Ablated.Yr=zeros(nBar,1); % List of relaxed vertices
Ablated.YrZ=zeros(nBar,1); % List of relaxed vertices only on Z
% Define number of ablated cells=nodes ablated
Set.Ablation=SetAblation(Cell,Set,X0); % Give node numbers of cells being ablated
% loop on nodal network
ncell=length(Cell);
Ablated.Nodal=zeros(nelem,1);
Ablated.Exist=false;
for e=1:nelem
    x01=X0(C(e,1),:);
    x02=X0(C(e,2),:);
    for jj=1:length(Set.Ablation)
        % cell center model/find the coordinates that incuded to the cell that we
        % want increase  contractility
        if (norm(X0(Cell{Set.Ablation(jj)}.xTopp,:)-x02)<tol && norm(X0(Cell{Set.Ablation(jj)}.xBott,:)-x01)<tol)||...     % for lateral surface
                (norm(X0(Cell{Set.Ablation(jj)}.xTopp,:)-x01)<tol && norm(X0(Cell{Set.Ablation(jj)}.xBott,:)-x02)<tol)
            Ablated.Nodal(e)=3;
            Ablated.Exist=true;
        elseif (norm(X0(Cell{Set.Ablation(jj)}.xTopp,:)-x02)<tol || norm(X0(Cell{Set.Ablation(jj)}.xTopp,:)-x01)<tol)    % for top surface
            Ablated.Nodal(e)=1;
            Ablated.Exist=true;
        elseif (norm(X0(Cell{Set.Ablation(jj)}.xBott,:)-x02)<tol || norm(X0(Cell{Set.Ablation(jj)}.xBott,:)-x01)<tol)    % for bottom surface
            Ablated.Nodal(e)=2;
            Ablated.Exist=true;
        end
    end
end
% Loop on vertex network
Ablated.Cell=cell(ncell,1);
Ablated.Cv=zeros(nBar,1);
ibar=0;
for c=1:ncell
    if isempty(Set.Ablation) || min(abs(c-Set.Ablation))>0  % Non-Ablated Cell
%         nv=length(Cell{c}.vtop);
%         Ablated.Cell{c}.Top=zeros(nv,1);                 %% Malik commented
%         Ablated.Cell{c}.Bot=zeros(nv,1);                 %% Malik commented
%         Ablated.Cell{c}.LatD=zeros(nv,1);                %% Malik commented
%         Ablated.Cell{c}.LatV=zeros(nv,1);                %% Malik commented
         Ablated.Cell{c}.eType=zeros(1,length(Cell{c}.eType));  %% Malik Added
         Ablated.Cell{c}.Ablated=0; % Ablated cell
         Ablated.Cell{c}.wounded=0;
    else
%         nv=length(Cell{c}.vtop);
%         Ablated.Cell{c}.Top=4*ones(nv,1); % Vertex top                   %% Malik commented
%         Ablated.Cell{c}.Bot=5*ones(nv,1); % Vertex Bottom                %% Malik commented
%         Ablated.Cell{c}.LatD=6*ones(nv,1); % Vertex Lateral Diagonal     %% Malik commented
%         Ablated.Cell{c}.LatV=7*ones(nv,1); % Vertex Lateral Vertical     %% Malik commented
        Ablated.Cell{c}.eType=Cell{c}.eType+3;  
        Ablated.Cell{c}.Ablated=1; % Ablated cell
        Ablated.Cell{c}.wounded=0;
%         Ablated.Cv(ibar+1:ibar+length(Cell{c}.vtop)*4)=1; % Mark all bars in cell as ablated
        Ablated.Cv(Cell{c}.CCv)=1;
        WCCv=Cv(Cell{c}.CCv,:);
        Yr=unique(WCCv);
        Ablated.Yr(Yr)=Yr;
        Ablated.Exist=true;
    end
    ibar=ibar+length(Cell{c}.vtop)*4;
end

%---------------------------------- Malik Added begin----------------------
for i=1:length(Set.Ablation)      % loop over ablated cells 
    cw=Set.Ablation(i);           % wounded Cell 
%     WCCv=Cv(Cell{cw}.CCv,:);      % Cv of the ablated Cell 
%     for j=1:size(WCCv,1)          % loop over the elements of ablated cells 
        for c=1:ncell             % loop over unwounded cells (can be improved using loop over neighbouring cells )
            unWCCv=Cv(Cell{c}.CCv,:);
            for jj=1:size(unWCCv,1) % loop over the elements of the unwounded cell (c) 
                if min(abs(c-Set.Ablation))>0 % Non-Ablated Cell
%                     if isequal(WCCv(j,:),unWCCv(jj,:)) || isequal(WCCv(j,:),fliplr(unWCCv(jj,:)))
                    if all(ismember(unWCCv(jj,:),[Cell{cw}.vtop;  Cell{cw}.vbot]))    
                        Ablated.Cell{c}.eType(jj)=Cell{c}.eType(jj)+6;
                        Ablated.Cell{c}.wounded=1;
                    end
                end
            end
        end
%     end
end 
%---------------------------------- Malik Added end -----------------------


% Mark vertices to be relaxed (regardless of whether they belong to wounded
% cells or not)
%--------------------------------------------------------------------------
% Malik Comment: I can not find in which context the following part is
% used, not sure if need to be modified. 
nVert=0;
for c=1:ncell
    Cell{c}.vtopR=zeros(length(Cell{c}.vtop),1);
    for j=1:length(Cell{c}.vtop)
        nVert=max(nVert,Cell{c}.vtop(j));
        if Ablated.Yr(Cell{c}.vtop(j))==1
            Cell{c}.vtopR(j)=1;
        end
        if Set.yRelaxationZ
            Ablated.YrZ(Cell{c}.vtop)=Cell{c}.vtop;
        end
    end
    Cell{c}.vbotR=zeros(length(Cell{c}.vbot),1);
    for j=1:length(Cell{c}.vbot)
        nVert=max(nVert,Cell{c}.vbot(j));
        if Ablated.Yr(Cell{c}.vbot(j))==1
            Cell{c}.vbotR(j)=1;
        end
    end
end


% Set Ablated.dofY with list of relaxed dof
dim=size(X0,2);
if Set.yRelaxation
    ny=length(Ablated.Yr(Ablated.Yr>0));
    Ablated.dofY=zeros(dim*ny,1);
    k=0;
    for i=1:length(Ablated.Yr)                             %%%%% remove basal dofs 
        if Ablated.Yr(i)>0
            if Y(Ablated.Yr(i),3)==0 && Set.FixedBasal==1
                Ablated.dofY(k+1:k+dim-1)=(Ablated.Yr(i)-1)*dim+(1:dim-1);
            else 
                Ablated.dofY(k+1:k+dim)=(Ablated.Yr(i)-1)*dim+(1:dim);
            end 
             Ablated.dofY(Ablated.dofY==0)=[];   
            k=k+dim;
        end
    end
    Ablated.Yr(Set.nvert+1:end)=[];
else
    Ablated.Yr=zeros(nVert,1); % Set to 0 if no yRelaxation
    Ablated.dofY=[];
end
% Relax on z. Set Ablated.dofYZ with list of relaxed dof
if Set.yRelaxationZ
    ny=length(Ablated.YrZ(Ablated.YrZ>0));
    Ablated.dofYZ=zeros(ny,1);
    k=0;
    for i=1:length(Ablated.YrZ)
        if Ablated.YrZ(i)>0
            k=k+1;
            Ablated.dofYZ(k)=Ablated.YrZ(i)*dim;
        end
    end
else
    Ablated.dofYZ=[];
end
Ablated.YrZ(Set.nvert+1:end)=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ablation=SetAblation(Cell,Set,X0)
% if isfield(Set,'Ablation')
%     Ablation=Set.Ablation;
%     return
% end
% Give node numbers of cells being ablated
if isempty(Set.CellCentres) && Set.AblationN>0 % Choose according to Set.nx and Set.ny
    if Set.AblationN>=(Set.nx-1)*(Set.ny-1)/2 || Set.nx<=2 || Set.ny<=2
        warning('Number of ablated cells in Set.AblationN too big. Set to 1');
        Set.AblationN=1;
    end
    c=floor((Set.nx-1)/2)*(Set.ny-1)+floor((Set.ny-1)/2)+mod(Set.ny-1,2); % center cell
    if Set.AblationN==0
        Ablation=[];
    elseif Set.AblationN==1
        Ablation=c;
    elseif Set.AblationN==2
        Ablation=[c c+1];
    elseif Set.AblationN==6
        Ablation=[c-Set.nx+1 c-Set.nx+2 c-1 c c+1 c+Set.nx-2];
    elseif Set.AblationN==7
        Ablation=[c-Set.nx+1 c-Set.nx+2 c-1 c c+1 c+Set.nx-2 c+Set.nx-1];
    else % case =2
        Ablation=(1:N)+c-floor(N/2);
    end   
elseif Set.AblationN>0 % Choose from center of patch initial configuration
    nnodes=size(X0,1)/2; % 
    X0=X0(1:nnodes,:); % Take only bottom
    xc=(min(X0(:,1:2))+max(X0(:,1:2)))/2; % Center
    nAblated=0;
    L=max(X0(:,1:2))-min(X0(:,1:2));
    dx0=sqrt(L(1)*L(2)/size(X0,1));
    dx=dx0;
    nodes=[];
    if Set.AblationN>=size(X0,1)/6
        warning('Number of ablated cells in Set.AblationN too big. Set to 1');
        Set.AblationN=1;
    end
    while nAblated<Set.AblationN
        nodesp=nodes;
        nodes=find(sqrt((X0(:,1)-xc(1)).^2+(X0(:,2)-xc(2)).^2)<dx);
        dx=dx+dx0;
        nAblated=length(nodes);
    end
    if nAblated>Set.AblationN
        nodesI=setdiff(nodes,nodesp);
        nodesd=zeros(size(nodesI));
        for i=1:length(nodesI)
            nodesd(i)=sqrt((X0(nodesI(i),1)-xc(1)).^2+(X0(nodesI(i),2)-xc(2)).^2);
        end
        nodesIS=sortrows([nodesd nodesI]);
        nodes=[nodesp' nodesIS(1:(Set.AblationN-length(nodesp)),2)']';
    end
    Ablation=nodes;
    % Convert node number to cell number
    ncell=length(Cell);
    Converted=zeros(nnodes,1);
    nConverted=0;
    for c=1:ncell
        if min(abs(Ablation-Cell{c}.xBott))==0 && ~Converted(Cell{c}.xBott)
            Ablation(Ablation==Cell{c}.xBott)=c;
            Converted(Cell{c}.xBott)=1;
            nConverted=nConverted+1;
            if nConverted==length(Ablation)
                break;
            end
        end
    end
else
    Ablation=[];
end
end

