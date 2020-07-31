function []=RemoveZigZag()


% if Set.smooth 
CvOld=Cv;
CellOld=Cell;
VRingBot1Old=VRingBot1;
VRingTop1Old=VRingTop1;
AblatedYrOld=Ablated.Yr;
[rowVBot,~]=find(VRingBot23(:,2)==0);
[rowVTop,~]=find(VRingTop23(:,2)==0);
RVBot=VRingBot1(rowVBot);
RVTop=VRingTop1(rowVTop);


% Update
for i=1:length(RVBot)
    Cv(CvOld>RVBot(i))=Cv(CvOld>RVBot(i))-1;
    Cv(CvOld>RVTop(i))=Cv(CvOld>RVTop(i))-1;
    Ablated.Yr(AblatedYrOld>RVBot(i))=Ablated.Yr(AblatedYrOld>RVBot(i))-1;
    Ablated.Yr(AblatedYrOld>RVTop(i))=Ablated.Yr(AblatedYrOld>RVTop(i))-1;
    VRingBot1(VRingBot1Old>RVBot(i))=VRingBot1(VRingBot1Old>RVBot(i))-1;
    VRingBot1(VRingBot1Old>RVTop(i))=VRingBot1(VRingBot1Old>RVTop(i))-1;
    VRingTop1(VRingTop1Old>RVBot(i))=VRingTop1(VRingTop1Old>RVBot(i))-1;
    VRingTop1(VRingTop1Old>RVTop(i))=VRingTop1(VRingTop1Old>RVTop(i))-1;
    for ncell=1:length(Cell)
        Cell{ncell}.vbot(CellOld{ncell}.vbot>RVBot(i))=Cell{ncell}.vbot(CellOld{ncell}.vbot>RVBot(i))-1;
        Cell{ncell}.vbot(CellOld{ncell}.vbot>RVTop(i))=Cell{ncell}.vbot(CellOld{ncell}.vbot>RVTop(i))-1;
        Cell{ncell}.vtop(CellOld{ncell}.vtop>RVBot(i))=Cell{ncell}.vtop(CellOld{ncell}.vtop>RVBot(i))-1;
        Cell{ncell}.vtop(CellOld{ncell}.vtop>RVTop(i))=Cell{ncell}.vtop(CellOld{ncell}.vtop>RVTop(i))-1;
    end
end 

%Remove zig-zag entities  
InsideBars=any(ismember(CvOld,[RVBot ;RVTop]),2);
IDBars=1:size(Cv,1);
IDBars=IDBars(InsideBars);
Cv(InsideBars,:)=[];
L.V.L(InsideBars)=[]; 
L0.V.L(InsideBars)=[]; 
Ln.V.L(InsideBars)=[]; 
Stress.V.s(InsideBars)=[];
Ablated.Cv(InsideBars,:)=[];

Y([VRingBot1Old(rowVBot); VRingTop1Old(rowVTop)],:)=[];
Y0([VRingBot1Old(rowVBot); VRingTop1Old(rowVTop)],:)=[];
Yn([VRingBot1Old(rowVBot); VRingTop1Old(rowVTop)],:)=[];

Ablated.Yr([VRingBot1Old(rowVBot); VRingTop1Old(rowVTop)])=[];   %??
Ablated.YrZ([VRingBot1Old(rowVBot) VRingTop1Old(rowVTop)])=[];  %??


VRingBot=[VRingBot1 VRingBot23];
VRingTop=[VRingTop1 VRingTop23];

Included=zeros(length(VRingBot1),1);
for ncell=1:length(Cell)
%         CellOld{ncell}.vbot
%         CellOld{ncell}.vtop
    RCellVBot=ismember(CellOld{ncell}.vbot,RVBot);
    RCellVTop=ismember(CellOld{ncell}.vtop,RVTop);    
    % Delet vertices 
    Cell{ncell}.vbot(RCellVBot)=[];
    Cell{ncell}.vtop(RCellVTop)=[];
    for j=1:length(IDBars)
        Cell{ncell}.CCv(CellOld{ncell}.CCv>IDBars(j))=Cell{ncell}.CCv(CellOld{ncell}.CCv>IDBars(j))-1;
    end 
    Cell{ncell}.CCv(ismember(CellOld{ncell}.CCv,IDBars))=[];
    Cell{ncell}.eType(ismember(CellOld{ncell}.CCv,IDBars))=[];
    Ablated.Cell{ncell}.eType(ismember(CellOld{ncell}.CCv,IDBars))=[];



    % create new elements
    for i=1:length(VRingBot1)
        % --------------------bot 
        if (VRingBot23(i,2)==0 && VRingBot23(i+1,2)==0) || ~ismember(VRingBot1(i),Cell{ncell}.vbot) || Included(i)>0
            continue 
        end 
        Included(i)=1;

        if VRingBot23(i,2)==0
            im1=i-1;
            nm1=0;
            while nm1==0;
                if im1<=0
                    im1=length(VRingBot1);
                end 
                cb1=VRingBot1(im1);
                nm1=VRingBot23(im1,2);
                im1=im1-1;
            end 

            ip1=i+1;
            np1=0;
            while np1==0;
                if ip1>length(VRingBot1)
                    ip1=1;
                end 
                cb2=VRingBot1(ip1);
                np1=VRingBot23(ip1,2);
                ip1=ip1+1;
            end 

            Cv=[Cv; [cb1 cb2]];
            Stress.V.s=[Stress.V.s; 0];
            Ablated.Cv=[Ablated.Cv; 0];
            L.V.L=[L.V.L; norm(Y(cb1,:)-Y(cb2,:))]; 
            L0.V.L=[L0.V.L; norm(Y(cb1,:)-Y(cb2,:))]; 
            Ln.V.L=[Ln.V.L; norm(Y(cb1,:)-Y(cb2,:))]; 
            Cell{ncell}.CCv=[Cell{ncell}.CCv size(Cv,1)];
            Cell{ncell}.eType=[Cell{ncell}.eType 1];
            Ablated.Cell{ncell}.eType=[Ablated.Cell{ncell}.eType 7];

        % --------------Top
            im1=i-1;
            nm1=0;
            while nm1==0;
                if im1<=0
                    im1=length(VRingTop1);
                end 
                ct1=VRingTop1(im1);
                nm1=VRingTop23(im1,2);
                im1=im1-1;
            end 

            ip1=i+1;
            np1=0;
            while np1==0;
                if ip1>length(VRingTop1)
                    ip1=1;
                end 
                ct2=VRingTop1(ip1);
                np1=VRingTop23(ip1,2);
                ip1=ip1+1;
            end 

            Cv=[Cv; [ct1 ct2]];
            Stress.V.s=[Stress.V.s; 0];
            Ablated.Cv=[Ablated.Cv; 0];
            L.V.L=[L.V.L; norm(Y(ct1,:)-Y(ct2,:))]; 
            L0.V.L=[L0.V.L; norm(Y(ct1,:)-Y(ct2,:))]; 
            Ln.V.L=[Ln.V.L; norm(Y(ct1,:)-Y(ct2,:))]; 
            Cell{ncell}.CCv=[Cell{ncell}.CCv size(Cv,1)];
            Cell{ncell}.eType=[Cell{ncell}.eType 2];
            Ablated.Cell{ncell}.eType=[Ablated.Cell{ncell}.eType 8];

            % ------------------ Vertical
    %         Cv=[Cv; [cb1 ct1]];
    %         Stress.V.s=[Stress.V.s; 0];
    %         Ablated.Cv=[Ablated.Cv; 0];
    %         L.V.L=[L.V.L; norm(Y(cb1,:)-Y(ct1,:))]; 
    %         L0.V.L=[L0.V.L; norm(Y(cb1,:)-Y(ct1,:))]; 
    %         Ln.V.L=[Ln.V.L; norm(Y(cb1,:)-Y(ct1,:))]; 
    %         
    %         Cv=[Cv; [cb2 ct2]];
    %         Stress.V.s=[Stress.V.s; 0];
    %         Ablated.Cv=[Ablated.Cv; 0];
    %         L.V.L=[L.V.L; norm(Y(cb2,:)-Y(ct2,:))]; 
    %         L0.V.L=[L0.V.L; norm(Y(cb2,:)-Y(ct2,:))]; 
    %         Ln.V.L=[Ln.V.L; norm(Y(cb2,:)-Y(ct2,:))]; 

            % ------------------Diagonal
            Cv=[Cv; [cb1 ct2]];
            Stress.V.s=[Stress.V.s; 0];
            Ablated.Cv=[Ablated.Cv; 0];
            L.V.L=[L.V.L; norm(Y(cb1,:)-Y(ct2,:))]; 
            L0.V.L=[L0.V.L; norm(Y(cb1,:)-Y(ct2,:))]; 
            Ln.V.L=[Ln.V.L; norm(Y(cb1,:)-Y(ct2,:))]; 
            Cell{ncell}.CCv=[Cell{ncell}.CCv size(Cv,1)];
            Cell{ncell}.eType=[Cell{ncell}.eType 3];
            Ablated.Cell{ncell}.eType=[Ablated.Cell{ncell}.eType 9];

        end 
    end 
end 
VRingBot(rowVBot,:)=[];
VRingTop(rowVTop,:)=[];
% end 