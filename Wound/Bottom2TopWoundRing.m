function [Cv,Cell,Y,Set,Ymid,Faces]=Bottom2TopWoundRing(Cell,Cb,Ct,Y,Set)


for n=1:ncell         % loop over Cells 
     if Ablated.Cell{n}.wounded ==0
         continue
     end 
     
for jj=2:length(vbot)
        jType=b2tType(jj,1);
        jb2t=b2tType(jj,2);
        j=vbot(jj);     
        if iType==1 && jType==1   
            %% CASE(1)  (||)
            [VFace,FaceSt]=CheckFace(i,j,Faces);
            if FaceSt==1 % the face already triangluized
                % build Lateral bars 
                Cv(k+1,:)=[j vtop(jb2t) 3];
                Cv(k+2,:)=[i VFace 3];
                Cv(k+3,:)=[j VFace 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(jb2t) VFace 3];
                cellCv(nk+1:nk+5)=k+1:k+5;
                k=k+5;
                nk=nk+5;
                % build tringlues 
                FaceTri=[i          j          VFace;
                         vtop(ib2t) i          VFace;
                         vtop(jb2t) vtop(ib2t) VFace;
                         j          vtop(jb2t) VFace];
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(jb2t)],:))/4;
                IAddY=1+size(Y,1);
                nMidY=nMidY+1;
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(jb2t) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[j vtop(jb2t) 3];
                Cv(k+2,:)=[i IAddY 3];
                Cv(k+3,:)=[j IAddY 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(jb2t) IAddY 3];
                cellCv(nk+1:nk+5)=k+1:k+5;
                k=k+5;
                nk=nk+5;
                FaceTri=[i          j          IAddY;
                         vtop(ib2t) i          IAddY;
                         vtop(jb2t) vtop(ib2t) IAddY;
                         j          vtop(jb2t) IAddY];
            end
        elseif iType==1 && jType==2   
            %% CASE(2) (| Y)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % build Lateral bars
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+6,:)=[j C 3];
                Cv(k+7,:)=[C vtop(ib2t+1) 3];
                Cv(k+8,:)=[C vtop(jb2t) 3];
                cellCv(nk+1:nk+8)=k+1:k+8;
                k=k+8;
                nk=nk+8;
                % build triangules 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         vtop(ib2t+1)  j            C;
                         vtop(jb2t)    vtop(ib2t+1) C];
                
                
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % build Lateral bars
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+6,:)=[j C 3];
                Cv(k+7,:)=[C vtop(ib2t+1) 3];
                Cv(k+8,:)=[C vtop(jb2t) 3];
                cellCv(nk+1:nk+8)=k+1:k+8;
                k=k+8;
                nk=nk+8;
                % build triangles 
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         vtop(ib2t+1)  j            C;
                         vtop(jb2t)    vtop(ib2t+1) C];
            end
        elseif iType==1 && jType==3   
         %% CASE(3) (| Y')
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+6,:)=[C vtop(ib2t+1) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build triangules
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         vtop(ib2t+1)  j            C];
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;  %??
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[j vtop(ib2t+1) 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+6,:)=[C vtop(ib2t+1) 3];                
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build triangules
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         vtop(ib2t+1)  j            C];
            end
        elseif iType==2 && jType==1   
           %% CASE(4) (Y |)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C=NonHorVer(i,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[vtop(ib2t) VFace 3];
                Cv(k+5,:)=[vtop(jb2t) VFace 3];
                Cv(k+6,:)=[j vtop(jb2t) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build Triangles 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         i              vtop(ib2t)  C];
            elseif FaceSt==0 % the face was not already triangluized
                % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(jb2t)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                nMidY=nMidY+1;
                FaceCount=FaceCount+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(jb2t) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[vtop(ib2t) IAddY 3];
                Cv(k+5,:)=[vtop(jb2t) IAddY 3];
                Cv(k+6,:)=[j vtop(jb2t) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6;
                % build Triangles 
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         i              vtop(ib2t)  C];
            end
       elseif iType==2 && jType==2   
           %% CASE(5) (Y Y)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C1=NonHorVer(i,2);
            C2=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) VFace 3];
                Cv(k+6,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(ib2t+1) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9;
                % build triangles 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2;
                         vtop(jb2t)    vtop(ib2t+1) C2];
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) IAddY 3];
                Cv(k+6,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(ib2t+1) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9;
                % build triangles 
                FaceTri=[i             j        IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2;
                         vtop(jb2t)    vtop(ib2t+1) C2];
            end
        elseif iType==2 && jType==3   
           %% CASE(6) (Y Y')
            [VFace,FaceSt]=CheckFace(i,j,Faces);
            C1=NonHorVer(i,2);
            C2=NonHorVer(j,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) VFace 3];
                Cv(k+6,:)=[vtop(ib2t+1) VFace 3];
                Cv(k+7,:)=[C2 vtop(ib2t+1) 3];
                cellCv(nk+1:nk+7)=k+1:k+7;
                k=k+7;
                nk=nk+7;               
                % build triangles 
                FaceTri=[i             j            VFace;
                         vtop(ib2t)    i            VFace;
                         vtop(ib2t+1)  vtop(ib2t)   VFace;
                         j             vtop(ib2t+1) VFace;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2];
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(ib2t) vtop(ib2t+1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(ib2t) vtop(ib2t+1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2t) 3];
                Cv(k+4,:)=[j vtop(ib2t+1) 3];
                Cv(k+5,:)=[vtop(ib2t) IAddY 3];
                Cv(k+6,:)=[vtop(ib2t+1) IAddY 3];
                Cv(k+7,:)=[C2 vtop(ib2t+1) 3];
                cellCv(nk+1:nk+7)=k+1:k+7;
                k=k+7;
                nk=nk+7;               
                % build triangles 
                FaceTri=[i             j            IAddY;
                         vtop(ib2t)    i            IAddY;
                         vtop(ib2t+1)  vtop(ib2t)   IAddY;
                         j             vtop(ib2t+1) IAddY;
                         i             vtop(ib2t)   C1;
                         vtop(ib2t+1)  j            C2];
            end
         elseif iType==3 && jType==1   
           %% CASE(7) (Y' |)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
             if jb2t ==1
                 jb2tt=size(vtop,2);
             else
                 jb2tt=jb2t;
             end 
             C1=NonHorVer(i,2);
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(jb2tt-1) 3];
                Cv(k+4,:)=[vtop(jb2tt) VFace 3];
                Cv(k+5,:)=[vtop(jb2tt-1) VFace 3];
                Cv(k+6,:)=[j vtop(jb2tt) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6; 
                %build triangles
                FaceTri=[i               j          VFace;
                        vtop(jb2tt-1)    i          VFace;
                        vtop(jb2tt)  vtop(jb2tt-1)   VFace;
                        j           vtop(jb2tt)     VFace;
                        i           vtop(jb2tt-1)     C1]; 
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(jb2tt) vtop(jb2tt-1)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(jb2tt) vtop(jb2tt-1) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(jb2tt-1) 3];
                Cv(k+4,:)=[vtop(jb2tt) IAddY 3];
                Cv(k+5,:)=[vtop(jb2tt-1) IAddY 3];
                Cv(k+6,:)=[j vtop(jb2tt) 3];
                cellCv(nk+1:nk+6)=k+1:k+6;
                k=k+6;
                nk=nk+6; 
                %build triangles 
                FaceTri=[i               j          IAddY;
                        vtop(jb2tt-1)    i          IAddY;
                        vtop(jb2tt)   vtop(jb2tt-1)   IAddY;
                        j            vtop(jb2tt)     IAddY;
                        i           vtop(jb2tt-1)     C1]; 
            end
        elseif iType==3 && jType==2   
         %% CASE(8) (Y' Y)
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
            C1=NonHorVer(i,2);
            C2=NonHorVer(j,2);
             if jb2t ==1
                 jb2tt=size(vtop,2)-1;
                 ib2tt=size(vtop,2)-2;
             elseif jb2t ==2
                 jb2tt=1;
                 ib2tt=size(vtop,2)-1;
             else 
                 ib2tt=jb2t-2;
                 jb2tt=jb2t-1;
             end 
            if FaceSt==1 % the face already triangluized
                % bulid bar elements 
                Cv(k+1,:)=[i VFace 3];
                Cv(k+2,:)=[j VFace 3];
                Cv(k+3,:)=[i vtop(ib2tt) 3];
                Cv(k+4,:)=[j vtop(jb2tt) 3];
                Cv(k+5,:)=[vtop(jb2tt) VFace 3];
                Cv(k+6,:)=[vtop(ib2tt) VFace 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(jb2tt) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9; 
                % build triangles 
               FaceTri=[i              j            VFace;
                        vtop(ib2tt)    i            VFace;
                        vtop(jb2tt)  vtop(ib2tt)    VFace;
                        j            vtop(jb2tt)    VFace;
                        i            vtop(ib2tt)     C1;
                        vtop(jb2tt)  j               C2;
                        vtop(jb2t)   vtop(jb2tt)     C2]; 
            elseif FaceSt==0 % the face was not already triangluized
                   % Add Vertex
                AddY=sum(Y([i j vtop(jb2tt) vtop(ib2tt)],:))/4;
                IAddY=1+size(Y,1);
                Y=[Y;AddY];
                FaceCount=FaceCount+1;
                nMidY=nMidY+1;
                Faces(FaceCount,:)=[i j vtop(jb2tt) vtop(ib2tt) IAddY];
                % bulid bar elements 
                Cv(k+1,:)=[i IAddY 3];
                Cv(k+2,:)=[j IAddY 3];
                Cv(k+3,:)=[i vtop(ib2tt) 3];
                Cv(k+4,:)=[j vtop(jb2tt) 3];
                Cv(k+5,:)=[vtop(jb2tt) IAddY 3];
                Cv(k+6,:)=[vtop(ib2tt) IAddY 3];
                Cv(k+7,:)=[j C2 3];
                Cv(k+8,:)=[C2 vtop(jb2tt) 3];
                Cv(k+9,:)=[C2 vtop(jb2t) 3];
                cellCv(nk+1:nk+9)=k+1:k+9;
                k=k+9;
                nk=nk+9; 
                % build triangles 
                FaceTri=[i              j            IAddY;
                         vtop(ib2tt)    i            IAddY;
                         vtop(jb2tt)  vtop(ib2tt)    IAddY;
                         j            vtop(jb2tt)    IAddY;
                         i            vtop(ib2tt)     C1;
                         vtop(jb2tt)  j               C2;
                         vtop(jb2t)   vtop(jb2tt)     C2]; 
            end
        elseif iType==3 && jType==3  
         %% CASE(9) (Y' Y')
            [VFace,FaceSt]=CheckFace(i,j,Faces); 
             k1=NonHorVer(i,2);
             k2=NonHorVer(j,2);
            if NonHorVer(i,2)~=NonHorVer(j,2) % Y'_Y'
                [r1,~]=find(NonHorVer(:,2)==k1);
                [r2,~]=find(NonHorVer(:,2)==k2);
                [ib2tt,~]=ismember(r1,vtop(1:end-1));  %??
                [jb2tt,~]=ismember(r2,vtop(1:end-1));  %??
                if FaceSt==1 % the face already triangluized
                    % bulid bar elements 
                    Cv(k+1,:)=[i VFace 3];
                    Cv(k+2,:)=[j VFace 3];
                    Cv(k+3,:)=[i r1(ib2tt) 3];
                    Cv(k+4,:)=[j r2(jb2tt) 3];
                    Cv(k+5,:)=[r2(jb2tt) VFace 3];
                    Cv(k+6,:)=[r1(ib2tt) VFace 3];
                    Cv(k+7,:)=[k2 r2(jb2tt) 3];
                    cellCv(nk+1:nk+7)=k+1:k+7;
                    k=k+7;
                    nk=nk+7;
                    %build triangles
                    FaceTri=[i             j           VFace;
                             r1(ib2tt)    i           VFace;
                             r2(jb2tt)  r1(ib2tt)    VFace;
                             j           r2(jb2tt)    VFace
                             i           r1(ib2tt)     k1;
                             r2(jb2tt)  j              k2]; 
                elseif FaceSt==0 % the face was not already triangluized
                       % Add Vertex
                    AddY=sum(Y([i j r2(jb2tt) r1(ib2tt)],:))/4;
                    IAddY=1+size(Y,1);
                    Y=[Y;AddY];
                    FaceCount=FaceCount+1;
                    nMidY=nMidY+1;
                    Faces(FaceCount,:)=[i j r2(jb2tt) r1(ib2tt) IAddY];
                    % bulid bar elements 
                    Cv(k+1,:)=[i IAddY 3];
                    Cv(k+2,:)=[j IAddY 3];
                    Cv(k+3,:)=[i r1(ib2tt) 3];
                    Cv(k+4,:)=[j r2(jb2tt) 3];
                    Cv(k+5,:)=[r2(jb2tt) IAddY 3];
                    Cv(k+6,:)=[r1(ib2tt) IAddY 3];
                    Cv(k+7,:)=[k2 r2(jb2tt) 3];
                    cellCv(nk+1:nk+7)=k+1:k+7;
                    k=k+7;
                    nk=nk+7;
                    %build triangles
                    FaceTri=[i             j           IAddY;
                             r1(ib2tt)    i           IAddY;
                             r2(jb2tt)  r1(ib2tt)    IAddY;
                             j         r2(jb2tt)    IAddY
                             i           r1(ib2tt)     k1;
                             r2(jb2tt)  j              k2];
                end
            else    % Y'
                   % build bar elemets 
                   Cv(k+1,:)=[i k1 3];
                   Cv(k+2,:)=[j k1 3];
                   cellCv(nk+1:nk+2)=k+1:k+2;
                    k=k+2;
                    nk=nk+2;
                   FaceTri=[i j k1];
                   
            end
                
        end
end
end 







end 