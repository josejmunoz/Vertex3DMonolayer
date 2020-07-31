function [VNew2Old,Y,Y0]=MapVertices(Y,T,TOld,Cv,CvOld,YOld,Ymid,Set,VNew2OldRing,Ablated,YmidOld,Y0)
%function - Relate the current vertices with one in the previuos mesh (checking triangles).
%         - Correct the position of mid-plane vertices in case they were found. 
%% Output: 
%     -VNew2Old: array of size(Y,1), such that VNew2Old(i)= the ID of the
%                  old vertex corresponding to vertex i  
%                  VNew2Old(i)=NaN --> vertex i not found 
VNew2Old=zeros(size(Y,1),1);
nTBV=size(Ymid(Ymid==0),1);
Y0Old=Y0;
Y0=zeros(size(Y));

%% Map only bottom and top vertices by comparing the nodal triangles 
% loop over the top and bottom vertices/Triangles
for i=1:nTBV 
    if Ablated.Yr(i)>0
        Row=VNew2OldRing(:,1)==i;
        if length(VNew2OldRing(Row,1))~=1
            error('Recheck !!!!')
        end 
        VNew2Old(i,:)=VNew2OldRing(Row,2);
%         VNew2Old(i,:)=nan;
    else
         ToldInTi=ismember(TOld,T(i,:)); % check commen nodes 
         ToldInTi=sum(ToldInTi,2);       
         [vID,~]=find(ToldInTi==3);      % find the one with 3 commen nodes 
        if ~isempty(vID) && length(vID) ==1 
            % This vertex/triangle is still there 
            VNew2Old(i,:)=vID;
            Y0(i,:)=Y0Old(vID,:);
        else 
            % it is not threr 
            VNew2Old(i,:)=NaN; 
            Y0(i,:)=Y(i,:);

        end 
    end 
end 


%% Map mid-plane vertices to the  previous  mesh
if Set.ModelTop==2 || Set.ModelTop==3 
    % loop over mid-plane Y
    % Describe the old mid-plane vertices in term of the connected top and bottom vertices (Faces)
    nYOld=size(YOld,1);         % total number of old vertives 
    nTBVOld=size(YmidOld(YmidOld==0),1);       % number of top and bottom vertices 
    IOld=nTBVOld+1:nYOld;   
    Faces=zeros(nYOld-nTBVOld,4);   
    for i=nTBVOld+1:nYOld
        RowsVi=any(CvOld==i,2);   % rows with i vertex (logical) 
        RowsVi=CvOld(RowsVi,:);   % rows with i vertex 
        VV=RowsVi(RowsVi~=i);  % vector of vertices connected to ith vertex 
        VV=unique(VV);
        if length(VV)==3
            Faces(i-nTBVOld,:)=[VV; NaN];
        else 
           Faces(i-nTBVOld,:)=VV;
        end 
    end 


    for i=nTBV+1:size(Y,1)
        % find vertices connected the (i) vertex 
        RowsVi=any(Cv==i,2);       % rows with i vertex (logical) 
        RowsVi=Cv(RowsVi,:);       % rows of Cv with i vertex 
        Vi=RowsVi(RowsVi~=i);      % vector of vertices connected to ith vertex (should be size 4) or three in the corners 3
        Vi=unique(Vi);
        ViNew2Old=VNew2Old(Vi);   % Transform Vi in terms of Old vertices
        Indicator=true(length(IOld),1);
        for ii=1:size(Vi,1)
            Indicator=any(Faces==ViNew2Old(ii),2) & Indicator;
        end 
        if length(IOld(Indicator))~=1
            VNew2Old(i,:)=NaN;
            Y0(i,:)=Y(i,:);

        else  
            VNew2Old(i,:)=IOld(Indicator);
            Y(i,:)=YOld(IOld(Indicator),:);
            Y0(i,:)=Y0Old(IOld(Indicator),:);

        end 
    end 
end 
end 