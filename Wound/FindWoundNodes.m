function [AbNBot,AbNTop,cNBot,cNTop]=FindWoundNodes(Cell,SetAblation)
%% Returns Data related to Ablated nodes:
% AbNBot = Ablated node numbers at bottom surface
% AbNTop = Ablated node numbers at top surface
% cNBot = Nodes forming the bottom ring. 
% cNTop = Nodes forming the top ring
AbNBot=zeros(size(SetAblation));
cNBot=zeros(length(SetAblation)*8,2);
AbNTop=zeros(size(SetAblation));
cNTop=zeros(length(SetAblation)*8,2);
k1=0;
k2=0;

for i=1:length(SetAblation)
    AbNBot(i)=Cell{SetAblation(i)}.xBott;
    AbNTop(i)=Cell{SetAblation(i)}.xTopp;
    
    xBot=Cell{SetAblation(i)}.xBot;
    cNBot(k1+1:k1+size(xBot,1),:)=xBot;
    k1=k1+size(xBot,1);

    xTop=Cell{SetAblation(i)}.xTop;
    cNTop(k2+1:k2+size(xTop,1),:)=xTop;
    k2=k2+size(xTop,1);
end 

cNBot(k1+1:end,:)=[];
cNBot=unique(cNBot);
cNBot=cNBot(~ismember(cNBot,AbNBot));

cNTop(k2+1:end,:)=[];
cNTop=unique(cNTop);
cNTop=cNTop(~ismember(cNTop,AbNTop));

end 