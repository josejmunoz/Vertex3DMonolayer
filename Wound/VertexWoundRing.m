function [VRingBot,VRingTop]=VertexWoundRing(NRingBot,NRingTop,AbNBot,AbNTop,T)
% Builds vertex wound ring 
% OUTPUT:
% VRingBot(i,:) = (v1, node1, node2) vertex and 2 nodes associated to
%                  vertex i along Bottom of wound rin. If node2=0, it means the vertex
%                  belongs to a triangle that only has one node at the wound edge.
IndexTri=1:size(T,1);  
%% bottom
VRingBot=zeros(length(NRingBot)*2,3);
k=0;    
% loop over nodal wound ring
for i=1:length(NRingBot)-1
    Ismem=ismember(T,NRingBot(i)); 
    TrsI=T(sum(Ismem,2)>0,:);                         % Tris connected to node i 
    IndexTriI=IndexTri(sum(Ismem,2)>0);               % Tris\ver Ids  connected to node i
    ConAb=any(ismember(TrsI,AbNBot),2);             % logical of Tris the with ablated nodes  
    RestOfTheRing=NRingBot; RestOfTheRing(i)=[];  % nodes ring wto node i 
    ConOthers=any(ismember(TrsI,RestOfTheRing),2);      % logical of Tris with RestOf..
    Vsolo=IndexTriI(ConAb & ~ConOthers);              % Tris\ver which are connected to node i and not to nodal wound ring 
    if ~isempty(Vsolo)
        for j=1:length(Vsolo)  % they might be more then one    
            VRingBot(k+1,:)=[Vsolo(j) NRingBot(i) 0];
            k=k+1;
        end 
    end
    ConJ=any(ismember(TrsI,NRingBot(i+1)),2);       % logical of Tris node i+1
    Vpair=IndexTriI(ConAb & ConJ);                      % Tris\ver which are connected to node i and j and ablated nodes
    
    if length(Vpair)~=1
        error('???????')
    end 
    VRingBot(k+1,:)=[Vpair NRingBot(i)  NRingBot(i+1)];
    k=k+1;
end     
VRingBot(k+1:end,:)=[];


%% Top
VRingTop=zeros(length(NRingTop)*2,3);
k=0;
% loop over nodal wound ring
for i=1:length(NRingTop)-1
    Ismem=ismember(T,NRingTop(i)); 
    TrsI=T(sum(Ismem,2)>0,:);                         % Tris connected to node i 
    IndexTriI=IndexTri(sum(Ismem,2)>0);               % Tris\ver Ids  connected to node i
    ConAb=any(ismember(TrsI,AbNTop),2);             % logical of Tris the with ablated nodes  
    RestOfTheRing=NRingTop; RestOfTheRing(i)=[];  % nodes ring wto node i 
    ConOthers=any(ismember(TrsI,RestOfTheRing),2);      % logical of Tris with RestOf..
    Vsolo=IndexTriI(ConAb & ~ConOthers);              % Tris\ver which are connected to node i and not to nodal wound ring 
    if ~isempty(Vsolo)
        for j=1:length(Vsolo)  % they might be more then one    
            VRingTop(k+1,:)=[Vsolo(j) NRingTop(i) 0];
            k=k+1;
        end 
    end
    ConJ=any(ismember(TrsI,NRingTop(i+1)),2);       % logical of Tris node i+1
    Vpair=IndexTriI(ConAb & ConJ);                      % Tris\ver which are connected to node i and j and ablated nodes
    
    if length(Vpair)~=1
        error('???????')
    end 
    VRingTop(k+1,:)=[Vpair NRingTop(i)  NRingTop(i+1)];
    k=k+1;
end
VRingTop(k+1:end,:)=[];

end 

