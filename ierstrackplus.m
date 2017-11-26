Xold=Xold-xhat;
L0=L0-ehat;

Xnew=Xold;
XnewG=Xnew(1:8);
XnewL=Xnew(9:end);
k=0;
for i=1:3:length(XnewL)-2
    k=k+1;
    Xtemp(k,:)=XnewL(i:i+2);
end


Lnew=L0;
k=0;
for i=1:5:length(Lnew)-4
    k=k+1;
    Ltemp(k,:)=Lnew(i:i+4);
end


clear Fc X Y Z AZ EL Aconstc Atargetc Bc
%从新的X及L变量中提取出迭代初值，生成新的雅可比矩阵A,
k=0;
for i=1:5:length(L0)-4
   k=k+1;
    L0XYZAZEL(k,:)=L0(i:i+4);
end


Hseq=[9,9,9,9,8,9,9,9,9,9,8,8];   %<=====================扣除了第五个点的第二个观测
k=0;
for i=1:length(Hseq)%size(L0XYZAZEL,1)
   % for j=1:length(Hseq)
    for j=1:Hseq(i)
        k=k+1;
        %disp(num2str(k));
        Fc(k).Fc=getFtc(Ltemp(k,1),Ltemp(k,2),Ltemp(k,3),XnewG(1),XnewG(2),XnewG(3),XnewG(5),...
        XnewG(6),XnewG(7),XnewG(4),Ltemp(k,4),Ltemp(k,5),XnewG(8),Xtemp(i,3),Xtemp(i,1),Xtemp(i,2));
       
        [pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA,pFpa,pFpb,pFpOE,pFpXYZ,pFpAZ,pFpEL]=...
        ipdtc(Ltemp(k,1),Ltemp(k,2),Ltemp(k,3),Ltemp(k,4),Ltemp(k,5),XnewG(1),XnewG(2),XnewG(3),...
        XnewG(4),XnewG(5),XnewG(6),XnewG(7),XnewG(8),Xtemp(i,1),Xtemp(i,2),Xtemp(i,3));
        %将L0XYZZAEL中的量依次赋给X/Y/Z/AZ/EL.
        X=L0XYZAZEL(k,1);Y=L0XYZAZEL(k,2);Z=L0XYZAZEL(k,3);
        AZ=L0XYZAZEL(k,4);EL=L0XYZAZEL(k,5);
        %组建矩阵Aconstc
        Aconstc(k).Aconstc=[pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA];
        %组建矩阵Atargetc
        Atargetc(i).Atargetc(j).At=[pFpa,pFpb,pFpOE];
        %组建矩阵Bc,表示pFpL
        Bc(k).Bc=[pFpXYZ,pFpAZ,pFpEL]; 
    end
end
Hend=i;
%add S
Sseq=[23,23,23,23,23];   %<=====================扣除了各方位的第首末观测点
%k=0;
for i=Hend+1:Hend+length(Sseq)%size(L0XYZAZEL,1)
   % for j=1:length(Hseq)
    for j=1:Sseq(i-Hend)
        k=k+1;
        %disp(num2str(k));
        Fc(k).Fc=getFtc(Ltemp(k,1),Ltemp(k,2),Ltemp(k,3),XnewG(1),XnewG(2),XnewG(3),XnewG(5),...
        XnewG(6),XnewG(7),XnewG(4),Ltemp(k,4),Ltemp(k,5),XnewG(8),Xtemp(i,3),Xtemp(i,1),Xtemp(i,2));
       
        [pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA,pFpa,pFpb,pFpOE,pFpXYZ,pFpAZ,pFpEL]=...
        ipdtc(Ltemp(k,1),Ltemp(k,2),Ltemp(k,3),Ltemp(k,4),Ltemp(k,5),XnewG(1),XnewG(2),XnewG(3),...
        XnewG(4),XnewG(5),XnewG(6),XnewG(7),XnewG(8),Xtemp(i,1),Xtemp(i,2),Xtemp(i,3));
        %将L0XYZZAEL中的量依次赋给X/Y/Z/AZ/EL.
        X=L0XYZAZEL(k,1);Y=L0XYZAZEL(k,2);Z=L0XYZAZEL(k,3);
        AZ=L0XYZAZEL(k,4);EL=L0XYZAZEL(k,5);
        %组建矩阵Aconstc
        Aconstc(k).Aconstc=[pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA];
        %组建矩阵Atargetc
        Atargetc(i).Atargetc(j).At=[pFpa,pFpb,pFpOE];
        %组建矩阵Bc,表示pFpL
        Bc(k).Bc=[pFpXYZ,pFpAZ,pFpEL]; 
    end
end


%利用Fc求出闭合差Omg
Omg=[];
for i=1:length(Fc)
   Omg=[Omg;Fc(i).Fc];
end
Omg=-Omg;
%组建矩阵B,将其Bc所有子阵组合起来
B=[];
for i=1:length(Bc)
    B=blkdiag(B,Bc(i).Bc);
end
%组建矩阵A,将Aconstc和Atargetc所有子阵结合起来
Ac=[];
for i=1:length(Aconstc)
    Ac=[Ac;Aconstc(i).Aconstc];
end
At=[];
for i=1:length(Atargetc)
    At00=[];
   for j=1:length(Atargetc(i).Atargetc)       
       At00=vertcat(At00,Atargetc(i).Atargetc(j).At) ; 
   end
   At0(i).At=At00;
end

At=[];
for i=1:length(At0)
    At=blkdiag(At,At0(i).At);
end
A=[Ac,At];

Acurve=A'*inv(B*QLL*B')*A;
bcurve=-A'*inv(B*QLL*B')*Omg;
xhat=Acurve\bcurve;
khat=inv(B*QLL*B')*(-Omg-A*xhat);
ehat=QLL*B'*inv(B*QLL*B')*(-Omg-A*xhat);
sgmhat2=-khat'*Omg/(size(A,1)-size(A,2));

%XVC=sgmhat2*inv(Acurve);


vx=Acurve*xhat-bcurve;


load  LL0.mat 
V=LL0-L0;
sg2=V'*inv(QLL)*V/(size(A,1)-size(A,2));

XVC=sg2*inv(Acurve);

disp('对圆心整体拟合的结果为：')
disp(num2str([Xold(1:3)']))
disp(num2str(sqrt(diag((XVC(1:3,1:3))*1e6))'))
disp('sigma of Post-fit residuals：')
disp(num2str(sqrt(V'*inv(QLL)*V/sum(diag(inv(QLL))))));

for i=1:8
    if i>4 && i<8
fprintf('%f ± %f\n',[Xold(i),sqrt(XVC(i,i))]*180/pi*3600); % mas
    elseif i==8
fprintf('%f ± %f\n',[Xold(i),sqrt(XVC(i,i))]*180/pi); % mas
    else
fprintf('%f ± %f\n',[Xold(i),sqrt(XVC(i,i))]);
    end
end
disp('Cor XY ')

disp(num2str(XVC(1,2)/sqrt(XVC(1,1)*XVC(2,2))))

disp(' ')

for i=1:8
    if i>4 && i<8
fprintf('%.1f(%.1f) &',[Xold(i),sqrt(XVC(i,i))]*180/pi*3600); % mas
    elseif i==8
fprintf('%.1f(%.1f) &',[Xold(i),sqrt(XVC(i,i))]*180*60/pi); % min
    else
fprintf('%.1f(%.1f) & ',[Xold(i),sqrt(XVC(i,i))]*1e3);
    end
end
%appended error of TP  
fprintf(' - &')
%pre-fit resdual(closing error)  
fprintf('%.1f &', presgm*1e3);
%post-fit resdual
fprintf('%.1f & ',sqrt(V'*inv(QLL)*V/sum(diag(inv(QLL))))*1e3);
%chi2
fprintf('%.1f & ',sg2);
fprintf('\\\\\n')


%求后验证残差Nabla


%sqrt(diag(XVC))

%norm(xhat,2)


% P=inv(QLL);
% Qnn=inv()
% P*ehat
% 
%XS0=inv(A'*inv(B*QLL*B')*A)*(A'*inv(B*QLL*B'))*QLL*(inv(A'*inv(B*QLL*B')*A)*(A'*inv(B*QLL*B')))';
%%

% 
% 
% 
% 
% Xnew=Xold;
% XnewG=Xnew(1:8);
% XnewL=Xnew(9:end);

% 
% 
% %for i=1:3:length(XnewL)-2
% k=0;
% for i=1:5:length(Lnew)-4
%     k=k+1;
%     Ltemp(k,:)=Lnew(i:i+4);
% end
% k=0;
% for i=1:3:length(XnewL)-2
%     k=k+1;
%     Xtemp(k,:)=XnewL(i:i+2);
% end
% 
% %Hseq=[]
% %%
% %从新的X及L变量中提取出迭代初值，生成新的雅可比矩阵B。
% % for i=1:size(Ltemp,1)
% %     %Fc(i).Fc=getF(X,Y,Z,XPR,YPR,ZPR,af,bt,gm,e,AZ,EL,OA,OE,a,b);
% %     %计算Omg
% %     Fc(i).Fc=getF(Ltemp(i,1),Ltemp(i,2),Ltemp(i,3),XnewG(1),XnewG(2),XnewG(3),XnewG(5),...
% %         XnewG(6),XnewG(7),XnewG(4),Ltemp(i,4),Ltemp(i,5),XnewG(8),Xtemp(i,3),Xtemp(i,1),Xtemp(i,2));
% %     %计算各项偏导数
% %     [pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA,pFpa,pFpb,pFpOE,pFpXYZ,pFpAZ,pFpEL]=...
% %         ipd(Ltemp(i,1),Ltemp(i,2),Ltemp(i,3),Ltemp(i,4),Ltemp(i,5),XnewG(1),XnewG(2),XnewG(3),...
% %         XnewG(4),XnewG(5),XnewG(6),XnewG(7),XnewG(8),Xtemp(i,1),Xtemp(i,2),Xtemp(i,3));
% %     %组建矩阵Bc,表示pFpL
% %     Bc(i).Bc=[pFpXYZ,pFpAZ,pFpEL];    
% % end
% 
% %从新的X及L变量中提取出迭代初值，生成新的雅可比矩阵A,
% k=0;
% for i=1:5:length(L0)-4
%    k=k+1;
%     L0XYZAZEL(k,:)=L0(i:i+4);
% end
% Hseq=[9,9,9,9,8,9,9,9,9,9,9,9];   %<=====================扣除了第五个点的第二个观测
% k=0;
% for i=1:length(Hseq)%size(L0XYZAZEL,1)
%    % for j=1:length(Hseq)
%     for j=1:Hseq(i)
%         k=k+1;
%         disp(num2str(k));
%         Fc(k).Fc=getF(Ltemp(i,1),Ltemp(i,2),Ltemp(i,3),XnewG(1),XnewG(2),XnewG(3),XnewG(5),...
%         XnewG(6),XnewG(7),XnewG(4),Ltemp(i,4),Ltemp(i,5),XnewG(8),Xtemp(i,3),Xtemp(i,1),Xtemp(i,2));
%         [pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA,pFpa,pFpb,pFpOE,pFpXYZ,pFpAZ,pFpEL]=...
%         ipd(Ltemp(i,1),Ltemp(i,2),Ltemp(i,3),Ltemp(i,4),Ltemp(i,5),XnewG(1),XnewG(2),XnewG(3),...
%         XnewG(4),XnewG(5),XnewG(6),XnewG(7),XnewG(8),Xtemp(i,1),Xtemp(i,2),Xtemp(i,3));
%         %将L0XYZZAEL中的量依次赋给X/Y/Z/AZ/EL.
%         X=L0XYZAZEL(k,1);Y=L0XYZAZEL(k,2);Z=L0XYZAZEL(k,3);
%         AZ=L0XYZAZEL(k,4);EL=L0XYZAZEL(k,5);
%         %组建矩阵Aconstc
%         Aconstc(k).Aconstc=[pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA];
%         %组建矩阵Atargetc
%         Atargetc(i).Atargetc(j).At=[pFpa,pFpb,pFpOE];
%         %组建矩阵Bc,表示pFpL
%         Bc(k).Bc=[pFpXYZ,pFpAZ,pFpEL]; 
%         
%         
%     end
% end
% %利用Fc求出闭合差Omg
% Omg=[];
% for i=1:length(Fc)
%    Omg=[Omg;Fc(i).Fc];
% end
% Omg=-Omg;

% 
% % %求QLL
% % k=0;
% % for i=1:length(H)
% %     for j=1:size(H(i).sgmXYZ,1)
% %         k=k+1;
% %         QLLc(k).QLLc=blkdiag(H(i).VNEU(j).VNEU,H(i).VAZEL(j).VAZEL);
% %     end
% % end
% % QLL=[];
% % for i=1:length(QLLc)
% %     QLL=blkdiag(QLL,QLLc(i).QLLc);
% % end