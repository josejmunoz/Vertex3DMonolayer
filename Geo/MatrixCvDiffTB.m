function  [Cv,Cell,Y,Set,Ymid,Faces]=MatrixCvDiffTB(Cell,Cb,Ct,Y,Set,Ablated)
% Build vertex connectivity with different top-bottom triangulation 
%% Output
% - Cv(i,:) = vertices connecting vertex bar element i
% - Cell  .CCv= Vector of the global IDs of the bar-elements in the matrix Cv
%         .eType= Vector of size(CCv) definse the type of the bar elements
%                = 1 bottom element , =2 top elements , =3 Lateral 
%         .Tri= matrix with triangles that defines the lateral surface
%         .nTri: Total number of surface triangles 
%         .nElem: Total number of bar-elements (size(Cv,1))
% - Ymid: vector for size(Y) 
%      - Ymid(i)=0 --> Top/bottom vertex 
%      - Ymid(i)=1 --> mid-plane vertex or vertex on the lateral side  

%% Initialize  
ncell=length(Cell);                      % number of Cells 
IncludeB=zeros(size(Cb,1),1);            % indicator for included  Bottom
IncludeT=zeros(size(Ct,1),1);            % indicator for included Top 
IntVerB=ones(size(Cb,1),1);                % indicator for vertices without cell (Bottom) 
IntVerT=ones(size(Ct,1),1);                % indicator for vertices without cell (Top)
NonHorVer=zeros(size(Y,1),2);            % vertices with non Horizontal connection 
nMidY=0;                                 % number of Add Vertices 
Top2BottomN=max(max(Cb));                % Num of nodes bottom
Top2BottomTri=size(Cb,1);                % Num of Tri bottom
cel=cell(ncell,1); 
Faces=zeros(4*ncell,5);                  % Face indicator 
FaceCount=0; 

if ~isempty(Ablated) && ~Ablated.Exist
    Ablated=[];
end 

%% Initial loop (The reason: is to define the type of the "bottom-top" connectivity.)
for n=1:ncell

    nTrsB=Cell{n}.vbot;                 % IDs of bottom Triangles
    nTrsT=Cell{n}.vtop;                 % IDs of Top  Tringles  

    if ~isempty(Ablated) && Ablated.Cell{n}.wounded
        IntVerB(Cell{n}.vbot(~ismember(Cell{n}.vbot,Ablated.VRingBot(1,:))))=0;                % mark these vertices as interior (connected)  
        IntVerT(Cell{n}.vtop(~ismember(Cell{n}.vtop,Ablated.VRingTop(1,:)))-Top2BottomTri)=0;  % mark these vertices as interior (connected)        
        TrsB=zeros(length(nTrsB),3);                   
        TrsT=zeros(length(nTrsT),3); 
        rVB=ismember(nTrsB,Ablated.VRingBot(:,1));
        rVT=ismember(nTrsT,Ablated.VRingTop(:,1));
        TrsB(~rVB,:)=Cb(nTrsB(~rVB),:);  TrsB(rVB,:)=nan(length(rVB(rVB)),3);
        TrsT(~rVT,:)=Ct(nTrsT(~rVT)-Top2BottomTri,:);  TrsT(rVT,:)=nan(length(rVT(rVT)),3);
    else 
        IntVerB(Cell{n}.vbot)=0;                % mark these vertices as interior (connected)  
        IntVerT(Cell{n}.vtop-Top2BottomTri)=0;  % mark these vertices as interior (connected)
        TrsB=Cb(nTrsB,:);                   % Nodes of bottom Triangles
        TrsT=Ct(nTrsT-Top2BottomTri,:);        % Nodes of Top Triangles
    end 
    b2tType=zeros(size(Cell{n}.vbot,1),2);  % 1 = vertical  2 = "Y"  3 ="Y'" Vector of connectivity type
    
    % Bottom-Top
    for iTb=1:size(TrsB,1)        % loop over bottom Triangles\Ver 
        Tb=TrsB(iTb,:);           % take the nodes of that Tri\ver (iTb)
        iTt=1;                    % counter fot top Tri\ver
        Found=0;                  
        % Search for the identical Tri on the top (while is not found) 
        while ~Found && iTt <= size(TrsT,1)        
            Tt=TrsT(iTt,:)-Top2BottomN;      % Take a top Tri\ver   (not sure if general)
            ismemb=Tb(:,ismember(Tb,Tt));    % check commen nodes 
            Found=length(ismemb)==3;         % if they are 3, than break
            iTt=iTt+1;
        end 
            
        %==================================================================    
        if Found  % It is found, -> vertical connectivity (|) Type=1
               IncludeB(nTrsB(iTb))=1;                   %include bottom vertex.  
               IncludeT(nTrsT(iTt-1)-Top2BottomTri)=1;   %include top vertex. 
               b2tType(iTb,:)=[1 iTt-1];            % define the type=1   
        %==================================================================
        elseif ~isempty(Ablated) && ismember(nTrsB(iTb),Ablated.VRingBot(:,1))
                %------- There are (4) cases for wounded vertices (not sure !!)
                RowEleV=ismember(Ablated.VRingBot(:,1),nTrsB(iTb));
                EleV=Ablated.VRingBot(RowEleV,[2 3]);
                RowEleVT=all(ismember(Ablated.VRingTop(:,[2 3]),EleV+Top2BottomN),2);
                if any(RowEleVT) 
                  %-- Case(1) vertical connectivity (|) Type=1 one-to-one
                    Index=1:length(nTrsT);
                    Indexlogic=ismember(nTrsT,Ablated.VRingTop(RowEleVT));
                    b2tType(iTb,:)=[1 Index(Indexlogic)];            % define the type=1   
                else 
                    NRelaxedBot=sum(ismember(nTrsB,Ablated.VRingBot(:,1)));
                    NRelaxedTop=sum(ismember(nTrsT,Ablated.VRingTop(:,1)));
                    if NRelaxedBot==2 && NRelaxedTop==0
                   %-- Case(2)  connectivity (Y') Type=3 (two relaxed to one unrelaxed)
                        b2tType(iTb,:)=[3 0];
                    else 
                    % find the traingle above the wound ele
                    TTriWithEleV=ismember(TrsT,EleV+Top2BottomN);
                    sTTriWithEleV=sum(TTriWithEleV,2)==2;
                    
                        if sum(sTTriWithEleV)==2 % Type=2 or Type=3
                    %-- Case(3) connectivity (Y) Type=2 (one relaxed to two (relaxed and unrelaxed))
                            % unrelaxed 
                             NnReV=sum(ismember(TrsT,Ablated.VRingTop(:,[2 3])),2)==2 & sTTriWithEleV;
                             NnReV=nTrsT(NnReV);
                             ReV=all(ismember(Ablated.VRingTop(:,[2 3]),TrsT(sTTriWithEleV,:)),2);
                             ReV=Ablated.VRingTop(ReV,1);
                             Index=1:length(nTrsT);
                             IIfle=Index(ismember(nTrsT,[NnReV ReV]));
                              if IIfle(1)==1 && IIfle(2)==length(nTrsT) % define the type=2
                                   b2tType(iTb,:)=[2 IIfle(1)]; % connectivity "Y"
                              else 
                                   b2tType(iTb,:)=[2 IIfle(2)]; % connectivity "Y"
                              end 
                              Ifle=nTrsT(IIfle);
                              if NonHorVer(Ifle(1),2)>0 &&  NonHorVer(Ifle(2),2)>0
                                  % modify The position of the added vertex
                                       Y(NonHorVer(Ifle(1),2),:)=(1/4).*Y(nTrsB(iTb),:)+(3/4).*Y(NonHorVer(Ifle(1),2),:);
%                                        IncludeB(nTrsB(iTb))=1;
                                       NonHorVer(nTrsB(iTb),:)=[1 NonHorVer(Ifle(1),2)];   
                              else 
                                       AddY=sum(Y([nTrsB(iTb) Ifle(1) Ifle(2)],:))/3;     % Position of the added vertex
                                       IAddY=1+size(Y,1);                                 % The ID of the added vertex
                                       nMidY=nMidY+1;                               % Update the Num of added vertices
                                       NonHorVer(Ifle,1)=1;                               % Update the vector of "non-horizental connectivity"
                                       NonHorVer(Ifle,2)=IAddY;                           % Update the vector of "non-horizental connectivity"
                                       NonHorVer(nTrsB(iTb),:)=[1 IAddY];                 % Update the vector of "non-horizental connectivity"
                                       Y=[Y;AddY];                                        % Add new Vertex
%                                        FaceCount=FaceCount+1;
%                                        Faces(FaceCount,:)=[nTrsB(iTb) 0 Ifle(1) Ifle(2) IAddY];
%                                        IncludeB(nTrsB(iTb))=1;
                                       IncludeT(NnReV-Top2BottomTri)=1;
                                      if IIfle(1)==1 && IIfle(2)==length(nTrsT)        % define the type=2
                                           b2tType(iTb,:)=[2 IIfle(1)]; % connectivity "Y"
                                      else 
                                           b2tType(iTb,:)=[2 IIfle(2)]; % connectivity "Y"
                                      end
                              end
                        else 
                            % Case(4) connectivity (Y) Type=2 (one relaxed to two, relaxed and unrelaxed)
                              b2tType(iTb,:)=[3 0];
                        end 
                    end 
                end 
        %==================================================================    
        else % It was not found. -> not vertical connectivity  Type=2 or Type=3
           Sim=ismember(TrsT-Top2BottomN,Tb); % check commen nodes between top-Tris and bottom tri
           I=sum(Sim,2);        
           [I,~]=find(I==2);    %  find the indices of top tir\ver with two commen nodes 
           TriW2com=TrsT(I,:);  %  Take nodes of top tir\vers with two commen nodes
           fle=zeros(size(TriW2com));   % indicator
           % Searsh for nodes which are commen between the top tringules
           % and excluded in bottom. if there is one, than we have "Y" shape.  
           for ii=1:size(TriW2com,1)   
                TriW2com2=TriW2com;
                TriW2com2(ii,:)=[];
                fle(ii,:)=ismember(TriW2com(ii,:),TriW2com2) & ~Sim(I(ii),:) ;
           end
           fle=sum(fle,2);
           [Ifle,~]=find(fle>0);
           IIfle=I(Ifle);      % the local indices of top two tri/ver connected to bottom tri/ver nTrsB(iTb) 
           Ifle=nTrsT(IIfle);  % the IDs of top two tri/ver connected to bottom tri/ver nTrsB(iTb)

           if ~isempty(Ifle)  % connectivity of type "Y" =2
               if IncludeT(Ifle(1)-Top2BottomTri) ==1 &&  IncludeT(Ifle(2)-Top2BottomTri) ==1
                  % Case(1) these two top vertices are already connected before 
                  % modify The position of the added vertex
                   Y(NonHorVer(Ifle(1),2),:)=(1/4).*Y(nTrsB(iTb),:)+(3/4).*Y(NonHorVer(Ifle(1),2),:);
                   IncludeB(nTrsB(iTb))=1;
                   NonHorVer(nTrsB(iTb),:)=[1 NonHorVer(Ifle(1),2)];   
                   
                  if IIfle(1)==1 && IIfle(2)==length(nTrsT) % define the type=2
                       b2tType(iTb,:)=[2 IIfle(1)]; % connectivity "Y"
                  else 
                       b2tType(iTb,:)=[2 IIfle(2)]; % connectivity "Y"
                  end 
               elseif IncludeT(Ifle(1)-Top2BottomTri) ==0 &&  IncludeT(Ifle(2)-Top2BottomTri) ==0
                   % Case(2) these two top vertices are not connected, (Add vertex in the middle)
                   AddY=sum(Y([nTrsB(iTb) Ifle(1) Ifle(2)],:))/3;     % Position of the added vertex
                   IAddY=1+size(Y,1);                                 % The ID of the added vertex
                   nMidY=nMidY+1;                               % Update the Num of added vertices
                   NonHorVer(Ifle,1)=1;                               % Update the vector of "non-horizental connectivity"
                   NonHorVer(Ifle,2)=IAddY;                           % Update the vector of "non-horizental connectivity"
                   NonHorVer(nTrsB(iTb),:)=[1 IAddY];                 % Update the vector of "non-horizental connectivity"
                   Y=[Y;AddY];                                        % Add new Vertex
%                    FaceCount=FaceCount+1;
%                    Faces(FaceCount,:)=[nTrsB(iTb) 0 Ifle(1) Ifle(2) IAddY];
                     
                   IncludeB(nTrsB(iTb))=1;
                   IncludeT(Ifle-Top2BottomTri)=1;
                  if IIfle(1)==1 && IIfle(2)==length(nTrsT)        % define the type=2
                       b2tType(iTb,:)=[2 IIfle(1)]; % connectivity "Y"
                  else 
                       b2tType(iTb,:)=[2 IIfle(2)]; % connectivity "Y"
                  end
               end 
           else 
        %==================================================================
                if ~isempty(Ablated) && Ablated.Cell{n}.wounded
                    % CASE (5)
                    RowEle=all(ismember(Ablated.VRingTop(:,[2 3])-Top2BottomN,Tb),2);
                    Index=1:length(nTrsT);
                    IIfle=Index(ismember(nTrsT,Ablated.VRingTop(RowEle,1)));
                    if length(IIfle)==2 
                       if IIfle(1)==1 && IIfle(2)==length(nTrsT) % define the type=2
                           b2tType(iTb,:)=[2 IIfle(1)]; % connectivity "Y"
                      else 
                           b2tType(iTb,:)=[2 IIfle(2)]; % connectivity "Y"
                      end 
                      Ifle=nTrsT(IIfle);
                      % find the opposite vertex 
                      OppVer=all(ismember(Ablated.VRingBot(:,[2 3]),TrsB(iTb,:)),2);
                        
                       AddY=sum(Y([nTrsB(iTb) Ifle(1) Ifle(2) Ablated.VRingBot(OppVer,1)],:))/4;     % Position of the added vertex
                       IAddY=1+size(Y,1);                                 % The ID of the added vertex
                       nMidY=nMidY+1;                               % Update the Num of added vertices
                       NonHorVer(Ifle,1)=1;                               % Update the vector of "non-horizental connectivity"
                       NonHorVer(Ifle,2)=IAddY;                           % Update the vector of "non-horizental connectivity"
                       NonHorVer(nTrsB(iTb),:)=[1 IAddY];                 % Update the vector of "non-horizental connectivity"
                       NonHorVer(Ablated.VRingBot(OppVer,1),:)=[1 IAddY];                 % Update the vector of "non-horizental connectivity"
                       Y=[Y;AddY];                                        % Add new Vertex
%                                        FaceCount=FaceCount+1;
%                                        Faces(FaceCount,:)=[nTrsB(iTb) 0 Ifle(1) Ifle(2) IAddY];
                       IncludeB(nTrsB(iTb))=1;
%                                IncludeT(Ifle-Top2BottomTri)=1;
                      if IIfle(1)==1 && IIfle(2)==length(nTrsT)        % define the type=2
                           b2tType(iTb,:)=[2 IIfle(1)]; % connectivity "Y"
                      else 
                           b2tType(iTb,:)=[2 IIfle(2)]; % connectivity "Y"
                      end
                    else 
                       b2tType(iTb,:)=[3 0]; % define the connectivity type "Y'" =3

                    end 
                    % Define NonHorVer for the opposite triangle 
                    
                    
                    
                    
                else 
                  %==================================================================
                   b2tType(iTb,:)=[3 0]; % define the connectivity type "Y'" =3
                end 
           end  
        end 
    cel{n}.b2tType=b2tType;
    end
end 
    
%%  Doubl-check 
% bottm-top (the reason is to define connectivity for missing boundary vertices)
IncludeB(IntVerB==1)=1;  % remove not needed vertices 
[Id,~]=find(IncludeB==0); % not included 
if ~isempty(Id)
    for j=1:length(Id)          %loop over unconnected bottom vertices
        ID=Id(j);               % first Tri/ver Id 
        Tb=Cb(ID,:);            % take the corresponding nodes of that Tri\ver
        Sim=ismember(Ct-Top2BottomN,Tb); % check commen nodes between the Tri (Tb) and all top tri
        I=sum(Sim,2);        
        [I,~]=find(I==2);    %  find the Id of top tri\ver with two commen nodes 
        TriW2com=Ct(I,:);    %  Taka nodes of top tir\ver with two commen nodes
        fle=zeros(size(TriW2com));   % indicator
        % Searsh for nodes which are commen between top tringules
        % and excluded in bottom. if there is one, than we have "Y" shape.  
        for ii=1:size(TriW2com,1)   
               TriW2com2=TriW2com;
               TriW2com2(ii,:)=[];
               fle(ii,:)=ismember(TriW2com(ii,:),TriW2com2) & ~Sim(I(ii),:) ;
        end
       fle=sum(fle,2);
       [Ifle,~]=find(fle>0);
       Ifle=I(Ifle);      % the IDs of top two tri/ver connected to bottom tri/ver Id(j)
       if ~isempty(Ifle)    
           if IncludeT(Ifle(1)) ==1 &&  IncludeT(Ifle(2)) ==1
               Y(NonHorVer(Ifle(1)+Top2BottomTri,2),:)=(1/4).*Y(ID,:)...
                                               +(3/4).*Y(NonHorVer(Ifle(1)+Top2BottomTri,2),:);
               IncludeB(ID)=1;
               NonHorVer(ID,:)=[1 NonHorVer(Ifle(1)+Top2BottomTri,2)];
%                TheRowIndiFace=any(Faces==NonHorVer(Ifle(1)+Top2BottomTri,2),2);
%                Faces(TheRowIndiFace,2)=ID;
               
           end 
       end 
    end 
end 

%% top-bottom  (the reason connect "Y'" in the corners) 
IncludeT(IntVerT==1)=1;  % remove not needed vertices 
[Id,~]=find(IncludeT==0);
for j=1:length(Id)          %loop over unconnected top vertices
    ID=Id(j);               % first Tri/ver Id 
    Tt=Ct(ID,:);         % take the corresponding nodes of that Tri\ver
    Sim=ismember(Cb+Top2BottomN,Tt); % check commen nodes between top-Tri and all bottom tris
    I=sum(Sim,2);        
    [I,~]=find(I==2);    %  find the Id of bottom tris\vers with two commen nodes 
    TriW2com=Cb(I,:)+Top2BottomN;    %  Taka nodes of bottom tir\vers with two commen nodes
    fle=zeros(size(TriW2com));   % indicator
    % Searsh for nodes which are commen between bottom tringules
    % and excluded in top. if there is one, than we have "Y'" shape.  
    for ii=1:size(TriW2com,1)   
       TriW2com2=TriW2com;
       TriW2com2(ii,:)=[];
       fle(ii,:)=ismember(TriW2com(ii,:),TriW2com2) & ~Sim(I(ii),:) ;
    end
    fle=sum(fle,2);
    [Ifle,~]=find(fle>0);
    Ifle=I(Ifle);     % the IDs of bottom two tri/ver connected to bottom tri/ver Id(j)
    if ~isempty(Ifle)   % just to make sure  
       AddY=sum(Y([ID+Top2BottomTri Ifle(1) Ifle(2)],:))/3;     % Position of the added vertex
       IAddY=1+size(Y,1);                                       % The ID of the added vertex
       nMidY=nMidY+1;                                     % Update the Num of added vertices
       NonHorVer(Ifle,1)=1;                                     % Update the vector of "non-horizental connectivity"
       NonHorVer(Ifle,2)=IAddY;                                 % Update the vector of "non-horizental connectivity"
       NonHorVer(ID+Top2BottomTri,:)=[1 IAddY];                 % Update the vector of "non-horizental connectivity"
       Y=[Y;AddY];                                       % Add new Vertex
       IncludeT(ID)=1;
       IncludeB(Ifle)=1;
    end 
end
    
    
    

%% Secound loop to bulid connectiviy 
Cv=zeros(ncell*24,3);             
k=0;                  % Bar element ID 
kk=0;                 % Tringules ID
for n=1:ncell         % loop over Cells 
    cellCv=zeros(80);      % IDs cellular bar elements
    nk=0;                  % local  counter ID 
    
    % Horizontal
    vtop=[Cell{n}.vtop' Cell{n}.vtop(1)]; 
    vbot=[Cell{n}.vbot' Cell{n}.vbot(1)];
    b2tType=[cel{n}.b2tType; cel{n}.b2tType(1,:)];
    
    % Top
    for v=1:length(vtop)-1               % loop around all vertices
        Cv(k+1,:)=[vtop(v) vtop(v+1) 2];   % horiz top
        k=k+1;                           % Global counter ID 
        cellCv(nk+1)=k;
        nk=nk+1;                         % local  counter ID 
        kk=kk+1;                         % number of triangles 
    end
    
    % Bottom
    for v=1:length(vbot)-1               % loop around all vertices
        Cv(k+1,:)=[vbot(v) vbot(v+1) 1];   % horiz top
        k=k+1;                           % Global counter ID 
        cellCv(nk+1)=k;
        nk=nk+1;                         % local  counter ID 
        kk=kk+1;                         % number of triangles 
    end
    

    
    
    %% Lateral
    % loop over bottom bar elements
    % initiate 
    iType=b2tType(1,1);    % vertex (i)
    ib2t=b2tType(1,2);     % vertex (j)
    i=Cell{n}.vbot(1);
    Cell{n}.Tri=[];         % matrix with triangles of the lateral surafce 
    
    for jj=2:length(vbot)
        jType=b2tType(jj,1);
        jb2t=b2tType(jj,2);
        j=vbot(jj);
        if ~isempty(Ablated) && ismember(i,Ablated.VRingBot(:,1)) && ismember(j,Ablated.VRingBot(:,1)) && Set.YmidWound==0 
            %% CASE(0)  (Wound Surface)
            [Cv,FaceTri,k,nk,cellCv]=WoundLateralFace(i,j,iType,ib2t,jType,jb2t,Cv,k,nk,vtop,cellCv,NonHorVer);
        elseif iType==1 && jType==1   
            %% CASE(1)  (||)
            [VFace,FaceSt]=CheckFace(i,j,Faces);
            if FaceSt==1 % the face already triangluized
                % build Lateral bars 
                Cv(k+1,:)=[j vtop(jb2t) 3];
                Cv(k+2,:)=[i VFace 3];
                Cv(k+3,:)=[j VFace 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(jb2t) VFace 3];
                cellCv(nk+1:nk+5)=k+1:k+5;
                k=k+5;
                nk=nk+5;
                % build tringlues 
                FaceTri=[i          j          VFace;
                         vtop(ib2t) i          VFace;
                         vtop(jb2t) vtop(ib2t) VFace;
                         j          vtop(jb2t) VFace];
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(jb2t)],:))/4;
                IAddY=1+size(Y,1);
                nMidY=nMidY+1;
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(jb2t) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[j vtop(jb2t) 3];
                Cv(k+2,:)=[i IAddY 3];
                Cv(k+3,:)=[j IAddY 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(jb2t) IAddY 3];
                cellCv(nk+1:nk+5)=k+1:k+5;
                k=k+5;
                nk=nk+5;
                FaceTri=[i          j          IAddY;
                         vtop(ib2t) i          IAddY;
                         vtop(jb2t) vtop(ib2t) IAddY;
                         j          vtop(jb2t) IAddY];
            end
        elseif iType==1 && jType==2   
            %% CASE(2) (| Y)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % build Lateral bars
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+6,:)=[j C 3];
                Cv(k+7,:)=[C vtop(ib2t+1) 3];
                Cv(k+8,:)=[C vtop(jb2t) 3];
                cellCv(nk+1:nk+8)=k+1:k+8;
                k=k+8;
                nk=nk+8;
                % build triangules 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         vtop(ib2t+1)  j            C;
                         vtop(jb2t)    vtop(ib2t+1) C];
                
                
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % build Lateral bars
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+6,:)=[j C 3];
                Cv(k+7,:)=[C vtop(ib2t+1) 3];
                Cv(k+8,:)=[C vtop(jb2t) 3];
                cellCv(nk+1:nk+8)=k+1:k+8;
                k=k+8;
                nk=nk+8;
                % build triangles 
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         vtop(ib2t+1)  j            C;
                         vtop(jb2t)    vtop(ib2t+1) C];
            end
        elseif iType==1 && jType==3   
         %% CASE(3) (| Y')
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+6,:)=[C vtop(ib2t+1) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build triangules
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         vtop(ib2t+1)  j            C];
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;  %??
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+6,:)=[C vtop(ib2t+1) 3];                
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build triangules
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         vtop(ib2t+1)  j            C];
            end
        elseif iType==2 && jType==1   
           %% CASE(4) (Y |)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C=NonHorVer(i,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(jb2t) VFace 3];
                Cv(k+6,:)=[j vtop(jb2t) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build Triangles 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         i              vtop(ib2t)  C];
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(jb2t)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(jb2t) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(jb2t) IAddY 3];
                Cv(k+6,:)=[j vtop(jb2t) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build Triangles 
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         i              vtop(ib2t)  C];
            end
       elseif iType==2 && jType==2   
           %% CASE(5) (Y Y)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C1=NonHorVer(i,2);
            C2=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) VFace 3];
                Cv(k+6,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(ib2t+1) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9;
                % build triangles 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2;
                         vtop(jb2t)    vtop(ib2t+1) C2];
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) IAddY 3];
                Cv(k+6,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(ib2t+1) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9;
                % build triangles 
                FaceTri=[i             j        IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2;
                         vtop(jb2t)    vtop(ib2t+1) C2];
            end
        elseif iType==2 && jType==3   
           %% CASE(6) (Y Y')
            [VFace,FaceSt]=CheckFace(i,j,Faces);
            C1=NonHorVer(i,2);
            C2=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) VFace 3];
                Cv(k+6,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+7,:)=[C2 vtop(ib2t+1) 3];
                cellCv(nk+1:nk+7)=k+1:k+7;
                k=k+7;
                nk=nk+7;               
                % build triangles 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2];
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) IAddY 3];
                Cv(k+6,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+7,:)=[C2 vtop(ib2t+1) 3];
                cellCv(nk+1:nk+7)=k+1:k+7;
                k=k+7;
                nk=nk+7;               
                % build triangles 
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2];
            end
         elseif iType==3 && jType==1   
           %% CASE(7) (Y' |)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
             if jb2t ==1
                 jb2tt=size(vtop,2);
             else
                 jb2tt=jb2t;
             end 
             C1=NonHorVer(i,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(jb2tt-1) 3];
                Cv(k+4,:)=[vtop(jb2tt) VFace 3];
                Cv(k+5,:)=[vtop(jb2tt-1) VFace 3];
                Cv(k+6,:)=[j vtop(jb2tt) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6; 
                %build triangles
                FaceTri=[i               j          VFace;
                        vtop(jb2tt-1)    i          VFace;
                        vtop(jb2tt)  vtop(jb2tt-1)   VFace;
                        j           vtop(jb2tt)     VFace;
                        i           vtop(jb2tt-1)     C1]; 
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(jb2tt) vtop(jb2tt-1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(jb2tt) vtop(jb2tt-1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(jb2tt-1) 3];
                Cv(k+4,:)=[vtop(jb2tt) IAddY 3];
                Cv(k+5,:)=[vtop(jb2tt-1) IAddY 3];
                Cv(k+6,:)=[j vtop(jb2tt) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6; 
                %build triangles 
                FaceTri=[i               j          IAddY;
                        vtop(jb2tt-1)    i          IAddY;
                        vtop(jb2tt)   vtop(jb2tt-1)   IAddY;
                        j            vtop(jb2tt)     IAddY;
                        i           vtop(jb2tt-1)     C1]; 
            end
        elseif iType==3 && jType==2   
         %% CASE(8) (Y' Y)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C1=NonHorVer(i,2);
            C2=NonHorVer(j,2);
             if jb2t ==1
                 jb2tt=size(vtop,2)-1;
                 ib2tt=size(vtop,2)-2;
             elseif jb2t ==2
                 jb2tt=1;
                 ib2tt=size(vtop,2)-1;
             else 
                 ib2tt=jb2t-2;
                 jb2tt=jb2t-1;
             end 
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2tt) 3];
                Cv(k+4,:)=[j vtop(jb2tt) 3];
                Cv(k+5,:)=[vtop(jb2tt) VFace 3];
                Cv(k+6,:)=[vtop(ib2tt) VFace 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(jb2tt) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9; 
                % build triangles 
               FaceTri=[i              j            VFace;
                        vtop(ib2tt)    i            VFace;
                        vtop(jb2tt)  vtop(ib2tt)    VFace;
                        j            vtop(jb2tt)    VFace;
                        i            vtop(ib2tt)     C1;
                        vtop(jb2tt)  j               C2;
                        vtop(jb2t)   vtop(jb2tt)     C2]; 
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(jb2tt) vtop(ib2tt)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(jb2tt) vtop(ib2tt) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2tt) 3];
                Cv(k+4,:)=[j vtop(jb2tt) 3];
                Cv(k+5,:)=[vtop(jb2tt) IAddY 3];
                Cv(k+6,:)=[vtop(ib2tt) IAddY 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(jb2tt) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9; 
                % build triangles 
                FaceTri=[i              j            IAddY;
                         vtop(ib2tt)    i            IAddY;
                         vtop(jb2tt)  vtop(ib2tt)    IAddY;
                         j            vtop(jb2tt)    IAddY;
                         i            vtop(ib2tt)     C1;
                         vtop(jb2tt)  j               C2;
                         vtop(jb2t)   vtop(jb2tt)     C2]; 
            end
        elseif iType==3 && jType==3  
         %% CASE(9) (Y' Y')
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
             k1=NonHorVer(i,2);
             k2=NonHorVer(j,2);
            if NonHorVer(i,2)~=NonHorVer(j,2) % Y'_Y'
                [r1,~]=find(NonHorVer(:,2)==k1);
                [r2,~]=find(NonHorVer(:,2)==k2);
                [ib2tt,~]=ismember(r1,vtop(1:end-1));  %??
                [jb2tt,~]=ismember(r2,vtop(1:end-1));  %??
                if FaceSt==1 % the face already triangluized
                    % bulid bar elements 
                    Cv(k+1,:)=[i VFace 3];
                    Cv(k+2,:)=[j VFace 3];
                    if length([i r1(ib2tt) 3])~=3
                        malik='debuge';
                    end 
                    Cv(k+3,:)=[i r1(ib2tt) 3];
                    Cv(k+4,:)=[j r2(jb2tt) 3];
                    Cv(k+5,:)=[r2(jb2tt) VFace 3];
                    Cv(k+6,:)=[r1(ib2tt) VFace 3];
                    Cv(k+7,:)=[k2 r2(jb2tt) 3];
                    cellCv(nk+1:nk+7)=k+1:k+7;
                    k=k+7;
                    nk=nk+7;
                    %build triangles
                    FaceTri=[i             j           VFace;
                             r1(ib2tt)    i           VFace;
                             r2(jb2tt)  r1(ib2tt)    VFace;
                             j           r2(jb2tt)    VFace
                             i           r1(ib2tt)     k1;
                             r2(jb2tt)  j              k2]; 
                elseif FaceSt==0 % the face was not already triangluized
                       % Add Vertex
                    AddY=sum(Y([i j r2(jb2tt) r1(ib2tt)],:))/4;
                    IAddY=1+size(Y,1);
                    Y=[Y;AddY];
                    FaceCount=FaceCount+1;
                    nMidY=nMidY+1;
                    if ~any(jb2tt)
                    jb2tt=ib2tt;
                    end 
                    Faces(FaceCount,:)=[i j r2(jb2tt) r1(ib2tt) IAddY];
                    % bulid bar elements 
                    Cv(k+1,:)=[i IAddY 3];
                    Cv(k+2,:)=[j IAddY 3];
                    Cv(k+3,:)=[i r1(ib2tt) 3];
                    Cv(k+4,:)=[j r2(jb2tt) 3];
                    Cv(k+5,:)=[r2(jb2tt) IAddY 3];
                    Cv(k+6,:)=[r1(ib2tt) IAddY 3];
                    Cv(k+7,:)=[k2 r2(jb2tt) 3];
                    cellCv(nk+1:nk+7)=k+1:k+7;
                    k=k+7;
                    nk=nk+7;
                    %build triangles
                    FaceTri=[i             j           IAddY;
                             r1(ib2tt)    i           IAddY;
                             r2(jb2tt)  r1(ib2tt)    IAddY;
                             j         r2(jb2tt)    IAddY
                             i           r1(ib2tt)     k1;
                             r2(jb2tt)  j              k2];
                end
            else    % Y'
                   % build bar elemets 
                   Cv(k+1,:)=[i k1 3];
                   Cv(k+2,:)=[j k1 3];
                   cellCv(nk+1:nk+2)=k+1:k+2;
                    k=k+2;
                    nk=nk+2;
                   FaceTri=[i j k1];
                   
            end
                
        end
        iType=jType;
        ib2t=jb2t;
        i=j;
        Cell{n}.Tri=[Cell{n}.Tri ;FaceTri];
        kk=kk+size(FaceTri,1);
    end
    cellCv(nk+1:end)=[];
    Cell{n}.CCv=cellCv;
    Cell{n}.eType=Cv(cellCv,3);
end 
%%

for n=1:ncell                            % loop over Cells 
Cell{n}.nTri=kk;
Cell{n}.nElem=k;
end 
Set.nMidY=nMidY;           
Set.nvert=size(Y,1);
Ymid=zeros(size(Y,1)-nMidY,1);  
YYmid=length(Ymid)+1:size(Y,1);
Ymid=[Ymid; YYmid'];  

Cv(k+1:end,:)=[];
% eType=Cv(:,3);
Cv=Cv(:,[1 2]);

Faces(FaceCount+1:end,:)=[];

% if Set.AblationN>0 
%     for n=1:length(Cell)
%         Ablated.Cell{n}.b2tType=cel{n}.b2tType;
%     end 
%     Ablated.NonHorVer=NonHorVer;
% end 

end 
%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%end 

function [VFace,FaceSt]=CheckFace(i,j,IndicatorFaces)
         TEST=ismember(IndicatorFaces(:,1:2),[i j]); 
         RowTEST=sum(TEST,2);
         [r,~]=find(RowTEST==2);
         if length(r)>1
             error('go to CheckFace')
         end 
         if isempty(r)
             % faces is not trianglized
             VFace=[];
             FaceSt=0;
         else  % faces is already trianglized
             VFace=IndicatorFaces(r,5);
             FaceSt=1;
         end 
end 


function [Cv,FaceTri,k,nk,cellCv]=WoundLateralFace(i,j,iType,ib2t,jType,jb2t,Cv,k,nk,vtop,cellCv,NonHorVer)
if iType==1 && jType==1   
    %% CASE(1)  (||)
        % build Lateral bars 
    Cv(k+1,:)=[j vtop(jb2t) 3];
    Cv(k+2,:)=[i vtop(jb2t) 3];
    cellCv(nk+1:nk+2)=k+1:k+2;
    k=k+2;
    nk=nk+2;
    % build tringlues 
    FaceTri=[i          j          vtop(jb2t);
             vtop(jb2t) vtop(ib2t)          i];

elseif iType==1 && jType==2   
    %% CASE(2) (| Y)
    C=NonHorVer(j,2);
        % build Lateral bars
        Cv(k+1,:)=[i vtop(ib2t+1) 3];
        Cv(k+2,:)=[j vtop(ib2t+1) 3];
        Cv(k+3,:)=[j C 3];
        Cv(k+4,:)=[C vtop(ib2t+1) 3];
        Cv(k+5,:)=[C vtop(jb2t) 3];
        cellCv(nk+1:nk+5)=k+1:k+5;
        k=k+5;
        nk=nk+5;
        % build triangules 
        FaceTri=[i             j            vtop(ib2t+1);
                 i        vtop(ib2t+1)      vtop(ib2t)
                 vtop(ib2t+1)  j            C;
                 vtop(jb2t)    vtop(ib2t+1) C];

elseif iType==1 && jType==3   
 %% CASE(3) (| Y')
    C=NonHorVer(j,2);
        % bulid bar elements 
        Cv(k+1,:)=[i vtop(ib2t+1) 3];
        Cv(k+2,:)=[j vtop(ib2t+1) 3];
        Cv(k+3,:)=[C vtop(ib2t+1) 3];
        cellCv(nk+1:nk+3)=k+1:k+3;
        k=k+3;
        nk=nk+3;
        % build triangules
        FaceTri=[i             j            vtop(ib2t+1);
                 i        vtop(ib2t+1)      vtop(ib2t);
                 vtop(ib2t+1)  j            C];
elseif iType==2 && jType==1   
   %% CASE(4) (Y |)
    C=NonHorVer(i,2);
        % bulid bar elements 
        Cv(k+1,:)=[i vtop(ib2t) 3];
        Cv(k+2,:)=[i vtop(jb2t) 3];
        Cv(k+3,:)=[j vtop(jb2t) 3];
        cellCv(nk+1:nk+3)=k+1:k+3;
        k=k+3;
        nk=nk+3;
        % build Triangles 
        FaceTri=[i             j            vtop(jb2t);
                 i             vtop(jb2t)   vtop(ib2t);
                 i              vtop(ib2t)  C];
elseif iType==2 && jType==2   
   %% CASE(5) (Y Y)
    C1=NonHorVer(i,2);
    C2=NonHorVer(j,2);
        % bulid bar elements 
        Cv(k+1,:)=[i vtop(ib2t) 3];
        Cv(k+2,:)=[i vtop(ib2t+1) 3];
        Cv(k+3,:)=[j vtop(ib2t+1) 3];
        Cv(k+4,:)=[j C2 3];
        Cv(k+5,:)=[C2 vtop(ib2t+1) 3];
        Cv(k+6,:)=[C2 vtop(jb2t) 3];
        cellCv(nk+1:nk+6)=k+1:k+6;
        k=k+6;
        nk=nk+6;
        % build triangles 
        FaceTri=[i             j            vtop(ib2t+1);
                 i     vtop(ib2t+1)         vtop(ib2t);
                 i             vtop(ib2t)   C1;
                 vtop(ib2t+1)  j            C2;
                 vtop(jb2t)    vtop(ib2t+1) C2];
elseif iType==2 && jType==3   
   %% CASE(6) (Y Y')
    C1=NonHorVer(i,2);
    C2=NonHorVer(j,2);
        % bulid bar elements 
        Cv(k+1,:)=[i vtop(ib2t) 3];
        Cv(k+2,:)=[i vtop(ib2t+1) 3];
        Cv(k+3,:)=[j vtop(ib2t+1) 3];
        Cv(k+4,:)=[C2 vtop(ib2t+1) 3];
        cellCv(nk+1:nk+4)=k+1:k+4;
        k=k+4;
        nk=nk+4;               
        % build triangles 
        FaceTri=[i             j            vtop(ib2t+1);
                 i     vtop(ib2t+1)         vtop(ib2t);
                 i             vtop(ib2t)   C1;
                 vtop(ib2t+1)  j            C2];
 elseif iType==3 && jType==1   
   %% CASE(7) (Y' |)
     if jb2t ==1
         jb2tt=size(vtop,2);
     else
         jb2tt=jb2t;
     end 
     C1=NonHorVer(i,2);
        % bulid bar elements 
        Cv(k+1,:)=[i vtop(jb2tt-1) 3];
        Cv(k+2,:)=[i vtop(jb2tt) 3];
        Cv(k+3,:)=[j vtop(jb2tt) 3];
        cellCv(nk+1:nk+3)=k+1:k+3;
        k=k+3;
        nk=nk+3; 
        %build triangles
        FaceTri=[i             j            vtop(jb2tt);
                 i     vtop(jb2tt)         vtop(jb2tt-1);
                 i           vtop(jb2tt-1)     C1]; 
elseif iType==3 && jType==2   
 %% CASE(8) (Y' Y)
    C1=NonHorVer(i,2);
    C2=NonHorVer(j,2);
     if jb2t ==1
         jb2tt=size(vtop,2)-1;
         ib2tt=size(vtop,2)-2;
     elseif jb2t ==2
         jb2tt=1;
         ib2tt=size(vtop,2)-1;
     else 
         ib2tt=jb2t-2;
         jb2tt=jb2t-1;
     end 
        % bulid bar elements 
        Cv(k+1,:)=[i vtop(ib2tt) 3];
        Cv(k+2,:)=[i vtop(jb2tt) 3];
        Cv(k+3,:)=[j vtop(jb2tt) 3];
        Cv(k+4,:)=[j C2 3];
        Cv(k+5,:)=[C2 vtop(jb2tt) 3];
        Cv(k+6,:)=[C2 vtop(jb2t) 3];
        cellCv(nk+1:nk+6)=k+1:k+6;
        k=k+6;
        nk=nk+6; 
        % build triangles 
       FaceTri=[i             j            vtop(jb2tt);
                i     vtop(jb2tt)          vtop(ib2tt);
                i            vtop(ib2tt)     C1;
                vtop(jb2tt)  j               C2;
                vtop(jb2t)   vtop(jb2tt)     C2]; 
elseif iType==3 && jType==3  
 %% CASE(9) (Y' Y')
     k1=NonHorVer(i,2);
     k2=NonHorVer(j,2);
    if NonHorVer(i,2)~=NonHorVer(j,2) % Y'_Y'
        [r1,~]=find(NonHorVer(:,2)==k1);
        [r2,~]=find(NonHorVer(:,2)==k2);
        [ib2tt,~]=ismember(r1,vtop(1:end-1));  %??
        [jb2tt,~]=ismember(r2,vtop(1:end-1));  %??
            % bulid bar elements 
            Cv(k+1,:)=[i r1(ib2tt) 3];
            Cv(k+2,:)=[i r2(jb2tt) 3];
            Cv(k+3,:)=[j r2(jb2tt) 3];
            Cv(k+4,:)=[k2 r2(jb2tt) 3];
            cellCv(nk+1:nk+4)=k+1:k+4;
            k=k+4;
            nk=nk+4;
            %build triangles
            FaceTri=[i             j          r2(jb2tt);
                     i        r2(jb2tt)       r1(ib2tt);
                     i           r1(ib2tt)     k1;
                     r2(jb2tt)  j              k2]; 
    else    % Y'
           % build bar elemets 
           Cv(k+1,:)=[i k1 3];
           Cv(k+2,:)=[j k1 3];
           cellCv(nk+1:nk+2)=k+1:k+2;
            k=k+2;
            nk=nk+2;
           FaceTri=[i j k1];

    end        
end
end 






