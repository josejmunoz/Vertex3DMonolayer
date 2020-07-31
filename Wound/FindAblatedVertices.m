function [AbVBot,AbVTop,AbVmid]=FindAblatedVertices(Ablated,Cv,AbNBot,AbNTop,T,SetModelTop)
% Returns list of vertices that are removed due to wound creation
IndexTri=1:size(T,1);
AbVBot=IndexTri(all(ismember(T,AbNBot),2));
AbVTop=IndexTri(all(ismember(T,AbNTop),2));
% Find mid-plane ablated Vertices
if SetModelTop==2 || SetModelTop==3 
    HV=Cv(Ablated.Cv==0,:);
    HV=unique(HV);
    Yr=Ablated.Yr(Ablated.Yr>0);
    AbVmid=Yr(~ismember(Yr,HV)& ~ismember(Yr,[AbVBot AbVTop]));
else 
    AbVmid=[];
end 