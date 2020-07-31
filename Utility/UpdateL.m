function L = UpdateL( C,Cv,L,Ln,Ld,ld,Mat,X,Xn,Y,Yn,Set,Ablated)
%UpdateL updates resting length as a function of iterative change of positions
% INPUT:
% C  = nodal connectivity
% Cv = vertex connectivity
% L  = rest length at previous iterations
% X  = updated nodal positions
% Y  = updated vertex positions
% Yp = vertex positions at previous iteration
% OUTPUT:
% L  = updated rest lengths
% gL = residual of evolution law dot L=gamma*(l-L)
%
% Nodal positions
for e=1:size(C,1)
    x1=X(C(e,1),:);
    x2=X(C(e,2),:);
    x1n=Xn(C(e,1),:);
    x2n=Xn(C(e,2),:);
    l=norm(x1-x2);
    ln=norm(x1n-x2n);
    lt=(1-Set.theta)*ln+Set.theta*l;
    wound=Ablated.Nodal(e);
    gamma=Mat.D.gamma;
    [~,gamma]=CheckAblated(wound,Set,gamma);
    if Mat.nDelay>0
        L.D.L(e)=Ln.D.L(e)+Set.dt*gamma*(ld.D.l(e)-Ld.D.L(e));
    else
        L.D.L(e)=(Ln.D.L(e)+Set.dt*gamma*(lt-(1-Set.theta)*Ln.D.L(e)))/(1+Set.dt*gamma*Set.theta);
    end
end
% Vertex positions
for e=1:size(Cv,1)
    y1=Y(Cv(e,1),:);
    y2=Y(Cv(e,2),:);
    y1n=Yn(Cv(e,1),:);
    y2n=Yn(Cv(e,2),:);
    l=norm(y1-y2);
    ln=norm(y1n-y2n);
    lt=(1-Set.theta)*ln+Set.theta*l;
    wound=Ablated.Cv(e);
    gamma=Mat.V.gamma;
    [~,gamma]=CheckAblated(wound,Set,gamma);
    if Mat.nDelay>0
        L.V.L(e)=Ln.V.L(e)+Set.dt*gamma*(ld.V.l(e)-Ld.V.L(e));
    else
%         [e gamma]
        L.V.L(e)=(Ln.V.L(e)+Set.dt*gamma*(lt-(1-Set.theta)*Ln.V.L(e)))/(1+Set.dt*gamma*Set.theta);
    end
end

end

