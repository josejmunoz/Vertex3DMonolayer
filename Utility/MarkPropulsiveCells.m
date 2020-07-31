function [Set]=MarkPropulsiveCells(X,Cell,Set)

if ~Set.Propulsion
    return
end 

Region1=zeros(length(Cell),1);
Region2=zeros(length(Cell),1);

k1=1;
k2=1;
for i=1:length(Cell)
    if  X(Cell{i}.xBott,1)>=Set.PropulsiveRegionX1(1) && X(Cell{i}.xBott,1)<=Set.PropulsiveRegionX1(2) &&...
        X(Cell{i}.xBott,2)>=Set.PropulsiveRegionY1(1) && X(Cell{i}.xBott,2)<=Set.PropulsiveRegionY1(2)
        Region1(k1+1)=i;
        k1=k1+1;
    end  
    if  X(Cell{i}.xBott,1)>=Set.PropulsiveRegionX2(1) && X(Cell{i}.xBott,1)<=Set.PropulsiveRegionX2(2) &&...
        X(Cell{i}.xBott,2)>=Set.PropulsiveRegionY2(1) && X(Cell{i}.xBott,2)<=Set.PropulsiveRegionY2(2)
        Region2(k2+1)=i;
        k2=k2+1;
    end 
end
Set.PropulsiveCells.Region1=Region1(1:k1);
Set.PropulsiveCells.Region2=Region2(1:k2);

    
end 