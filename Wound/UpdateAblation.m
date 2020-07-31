function Ablated=UpdateAblation(Ablated,Set,Cell,Cv,C,Ymid)
%% Updates information in struct Ablated:
%     Ablated.Yr
%     Ablated.dofY
%     Ablated.Nodal
%     Ablated.Cv
%     Ablated.YrZ
%     Ablated.dofYZ
%     Ablated.Cell=cell(ncell,1);
%     Ablated.WTri
%
if ~Ablated.Exist
    Ablated.Yr=zeros(Set.nvert,1);
    Ablated.dofY=[];
    Ablated.Nodal=zeros(size(C,1),1);
    Ablated.Cv=zeros(size(Cv,1),1);
    Ablated.YrZ=zeros(size(Ablated.Yr)); %% ????
    Ablated.dofYZ=[]; 
    ncell=length(Cell);
    Ablated.Cell=cell(ncell,1);
    for c=1:ncell             % loop over unwounded cells (can be improved using loop over neighbouring cells )
        Ablated.Cell{c}.eType=zeros(size(Cell{c}.CCv));
        Ablated.Cell{c}.Ablated=0; % Ablated cell
        Ablated.Cell{c}.wounded=0;
    end
    return
end

% Update Wounded verices 
Yr=[1:Set.nvert]';
VRingBot=Ablated.VRingBot(:,1);
VRingTop=Ablated.VRingTop(:,1);
Yr(~ismember(Yr,VRingBot) & ~ismember(Yr,VRingTop))=0;

if Set.FixedBasal==0
    dofY=3.*(kron(Yr(Yr>0),[1;1;1])-1);
    dofY=kron(ones(length(Yr(Yr>0)),1),[1;2;3])+dofY;
else 
    dofYBot=3.*(kron(VRingBot(:,1),[1;1])-1);
    dofYBot=kron(ones(length(VRingBot(:,1)),1),[1;2])+dofYBot;
    dofYTop=3.*(kron(VRingTop(:,1),[1;1;1])-1);
    dofYTop=kron(ones(length(VRingTop(:,1)),1),[1;2;3])+dofYTop;
    dofY=[dofYBot; dofYTop];
    Ablated.dofY=unique(dofY);
end

ncell=length(Cell);
if ~isfield(Ablated,'Cell')
    Ablated.Cell=cell(ncell,1);
end

YYmid=Ymid(Ymid>0);
wYmid=zeros(size(YYmid));
% find mid-plane wounded vertices
for c=1:ncell 
    wTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1);YYmid]),2) |...
            all(ismember(Cell{c}.Tri,[Ablated.VRingTop(:,1);YYmid]),2);
    wTri=Cell{c}.Tri(wTriRow,:); 
    wTri=reshape(wTri,[],1);
    wYmidRow=ismember(YYmid,wTri);
    wYmid(wYmidRow)=YYmid(wYmidRow);
end

WTri=zeros(3*length(Ablated.VRingBot(:,1)),3);
Tnt=0;

for c=1:ncell     
    unWCCv=Cv(Cell{c}.CCv,:);
    Ablated.Cell{c}.eType=zeros(size(Cell{c}.CCv));
    Ablated.Cell{c}.wounded=0;
    Ablated.Cell{c}.Ablated=0; % Ablated cell
    for jj=1:size(unWCCv,1) % loop over the elements of the unwounded cell (c) 
        if all(ismember(unWCCv(jj,:),[Ablated.VRingBot(:,1); Ablated.VRingTop(:,1); wYmid]));     
            Ablated.Cell{c}.eType(jj)=Cell{c}.eType(jj)+6;
            Ablated.Cell{c}.wounded=1;
        end
    end
    
    if Ablated.Cell{c}.wounded==1
        wwTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1); Ablated.VRingTop(:,1); wYmid]),2);
        wwTri=Cell{c}.Tri(wwTriRow,:); wwTri=flip(wwTri,2);
        nt=size(wwTri,1);
        WTri(Tnt+1:Tnt+nt,:)=wwTri;
        Tnt=Tnt+nt;
    end
end

WTri(Tnt+1:end,:)=[];

Ablated.WTri=WTri;
Ablated.Yr=Yr;
Ablated.dofY=dofY;
Ablated.Nodal=zeros(size(C,1),1);
Ablated.Cv=zeros(size(Cv,1),1);
Ablated.YrZ=zeros(size(Yr)); %% ????
Ablated.dofYZ=[];            %% ????

% % find mid-plane vertices on the wound-edge
% if Set.ModelTop==3
%      wYmid=zeros(5*length(Ablated.VRingBot(:,1)),1);
%      nv=0;
%     for c=1:ncell
%         if Ablated.Cell{c}.wounded==0
%             continue
%         else 
%             wTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1);Ymid(Ymid>0)]),2) |...
%                     all(ismember(Cell{c}.Tri,[Ablated.VRingTop(:,1);Ymid(Ymid>0)]),2);
% 
%             wYmid=[wYmid; unique(Cell{c}.Tri(wTriRow,:))];
%         end 
%         wYmid=wYmid(ismember(wYmid,Ymid(Ymid>0)));
%     end
% else 
%     wYmid=[];
% end 
% 
% 
% 
% 
% 
% WTri=zeros(3*length(Ablated.VRingBot(:,1)),3);
% Tnt=0;
% for c=1:ncell
%     if Ablated.Cell{c}.wounded==0
%         continue
%     else 
%         if Set.ModelTop==3
%             wTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1);Ablated.VRingTop(:,1)]),2) |... % Lateral 
%                     all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1);Ymid(Ymid>0)]),2) |...
%                     all(ismember(Cell{c}.Tri,[Ablated.VRingTop(:,1);Ymid(Ymid>0)]),2);
%                 
%             wTriRow=all(ismember(Cell{c}.Tri,Cell{c}.Tri(wTriRow,:)));
%             
%         else
%             wTriRow=all(ismember(Cell{c}.Tri,[Ablated.VRingBot(:,1);Ablated.VRingTop(:,1)]),2);    
%         end 
%         wTri=Cell{c}.Tri(wTriRow,:); wTri=flip(wTri,2);
%         nt=size(wTri,1);
%         WTri(Tnt+1:Tnt+nt,:)=wTri;
%         Tnt=Tnt+nt;
%     end 
% end 
% WTri(Tnt+1:end,:)=[];

