function [n1,n2,n3,Flag]=OrientNodes(n1,n2,n3,X)
% Make sure that triangle [n1 n2 n3] is counterclockwise
% OUTPUT:
% Flag = true : triangle was oritnted clockwise, and changed to counterclockwise
x1=X(n1,:);
x2=X(n2,:);
x3=X(n3,:);
x12=x2-x1;
x13=x3-x1;
Flag=false;
if x12(1)*x13(2)-x12(2)*x13(1)<0
    aux=n2;
    n2=n3;
    n3=aux;
    Flag=true;
end
end
