function [NRingBot,NRingTop]=NodalWoundRing(cNBot,cNTop,AbNBot,AbNTop,T)
%% Builds wound ring at top and bottom surface
% INPUT:
% AbNTop = Ablated nodes at top surface
% AbNBot = Ablated nodes at bottom surface
% cNTop = nodes at wound ring of top surface
% cNBot = nodes at wound ring of bottom surface

%% bottom
NRingBot=zeros(length(cNBot)+1,1); 

% Strat from the first  
NRingBot(1)=cNBot(1);
% find the two neighbouring cells 
Ismem=ismember(T,cNBot(1));  % Triangles connected to 1st node of bottom ring  
TrsI=T(sum(Ismem,2)>0,:); % Triangles connected to 1st node of bottom ring  
NeighNodes=unique(TrsI); % Set of nodes connected to 1st node of bottom ring
NeighNodes=NeighNodes(~ismember(NeighNodes,AbNBot) & ...
                    ~ismember(NeighNodes,cNBot(1)) & ismember(NeighNodes,cNBot)); % nodes at wound ring connected to 1st node of bottom ring
                
% take the first as the next 
if length(NeighNodes)~=2 
    error('Error building wound ring')
else 
    NRingBot(2)=NeighNodes(1);
end

% find the other node 
NPrev=[3 1 2];
for i=2:length(NRingBot)-1
    Ismem=ismember(T,NRingBot(i));    
    TrsI=T(sum(Ismem,2)>0,:);
    NeighNodes=unique(TrsI);
    Next=NeighNodes(~ismember(NeighNodes,AbNBot) & ...
                    ~ismember(NeighNodes,[NRingBot(i) NRingBot(i-1)]) & ismember(NeighNodes,cNBot));
    if length(Next)>1 % discard also not in contact with 2 previous nodes
        if i==2
            Next=Next(~ismember(Next,NRingBot(1)));
        else
            Next=Next(~ismember(Next,NRingBot(i-2)));
        end
    end
    if length(Next)>1 % Choose correct one from 2 possible nodes at wound ring
      Te= TrsI(sum(ismember(TrsI,[Next;NRingBot(i)]),2)==3,:);
      Next(1)=Te(NPrev(Te==NRingBot(i)));
    end
    NRingBot(i+1)=Next(1);
end 


%% Top
NRingTop=zeros(length(cNTop)+1,1); 

% Strat from the first  
NRingTop(1)=cNTop(1);
% find the two neighbouring cells 
Ismem=ismember(T,cNTop(1));    
TrsI=T(sum(Ismem,2)>0,:);
NeighNodes=unique(TrsI);
NeighNodes=NeighNodes(~ismember(NeighNodes,AbNTop) & ...
                      ~ismember(NeighNodes,cNTop(1)) & ismember(NeighNodes,cNTop));
                
% take the first as the next 
if length(NeighNodes)~=2 
    error('Error building wound ring')
else 
    NRingTop(2)=NeighNodes(1);
end

% find the other node 
for i=2:length(NRingTop)-1
    Ismem=ismember(T,NRingTop(i));    
    TrsI=T(sum(Ismem,2)>0,:);
    NeighNodes=unique(TrsI);
    Next=NeighNodes(~ismember(NeighNodes,AbNTop) & ...
                    ~ismember(NeighNodes,[NRingTop(i) NRingTop(i-1)]) & ismember(NeighNodes,cNTop));
    if length(Next)>1 % discard also not in contact with 2 previous nodes
        if i==2
            Next=Next(~ismember(Next,NRingTop(1)));
        else
            Next=Next(~ismember(Next,NRingTop(i-2)));
        end
    end
    if length(Next)>1 % Choose correct one from 2 possible nodes at wound ring
      Te= TrsI(sum(ismember(TrsI,[Next;NRingTop(i)]),2)==3,:);
      Next(1)=Te(NPrev(Te==NRingTop(i)));
    end
    NRingTop(i+1)=Next(1);
end 
% Ensure direction of ring is counterclock wise 
Te=T(sum(ismember(T,NRingBot(1:2)),2)==2,:);
Te=Te(sum(ismember(Te,AbNBot),2)==0,:); % Triangle at wound ring not containing ablted nodes
Pos=find(Te==NRingBot(1) | Te==NRingBot(2));
if (Te(Pos(1))==NRingBot(1) && Pos(2)-Pos(1)==1) || (Te(Pos(1))==NRingBot(2) && Pos(2)-Pos(1)==2)
    NRingBot=NRingBot(end:-1:1);
end
Te=T(sum(ismember(T,NRingTop(1:2)),2)==2,:);
Te=Te(sum(ismember(Te,AbNTop),2)==0,:); % Triangle at wound ring not containing ablted nodes
Pos=find(Te==NRingTop(1) | Te==NRingTop(2));
if (Te(Pos(1))==NRingTop(1) && Pos(2)-Pos(1)==1) || (Te(Pos(1))==NRingTop(2) && Pos(2)-Pos(1)==2)
    NRingTop=NRingTop(end:-1:1);
end
