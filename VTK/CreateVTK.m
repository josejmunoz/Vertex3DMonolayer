function []=CreateVTK(X,X0,C,Cv,Cell,Ablated,L,L0,Stress,Set,Y,Y0,NameFile)
t=Set.iIncr;
if Set.OutputVTK
    if Set.iIncr==0
        InitVTK();
    end
    Ablated.Exist=0;
    CreateVtkNod(t,X0,X,C,Ablated,L.D.L,L0.D.L,Stress.D.s,Set,NameFile)
    CreateVtkVert(t,Y,Cell,Cv,L.V.L,L0.V.L,Stress.V.s,Ablated,Y0,Set,NameFile)
    CreateVtkTri(t,Y,Cell,Ablated,Y0,X,X0,NameFile)
end
