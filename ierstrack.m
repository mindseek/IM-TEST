nsi=0;%+0.00046204+0.00020898+0.00005898;%+0.00062142+0.00025334;%+0.000667+0.00025942;%+0.00032963+0.00016902+4.9568e-05+0.00041775; %
clearvars -except nsi;clc;%close all
flnm=dir('H/');
flnm(1:2)=[];
ks=0;
%依H中坐标点的顺序给出的EL和AZ的初值。
EL=[10,20,30,40,50,60,70,80,88];
AZ=[70,70,130,130,190,190,250,250,310,310,10,10];  
%从H文件夹中的文件中读取各个靶标点的坐标方位及其误差到变量H中。

%nsi=0+0.00032963+0.00016902+4.9568e-05; %<=====================
for i=1:length(flnm)
     if strcmp(flnm(i).name(1),'H')
         ks=ks+1;
         TEMP=load(['H/' flnm(i).name]);
         H(ks).XYZ=TEMP(:,1:3);%+nsi*randn(size(TEMP,1),3
         %H(ks).XYZ=TEMP(:,1:3)-nsi*ones(size(TEMP,1),3);
         H(ks).sgmXYZ=TEMP(:,4:6);
         %------------------
         H(ks).NEU=H(ks).XYZ;
         for j=1:size(H(ks).sgmXYZ,1)
         H(ks).VXYZ(j).VXYZ=diag(H(ks).sgmXYZ(j,:).^2+nsi.^2);
         H(ks).VNEU(j).VNEU=H(ks).VXYZ(j).VXYZ;
         H(ks).VAZEL(j).VAZEL=diag([((5e-3)*pi/180).^2+(nsi/norm([15.2,2.5])).^2,((5e-3)*pi/180).^2+(nsi/norm([15.2,2.5])).^2]);%<====== 天线方位俯仰指向设计精度0.005度
         end
         H(ks).AZEL=[AZ(i)*ones(9,1) EL'];
     end
end

%% 剔除第三次方位 20度位置的一次观测

i=[5 11 12];j=[2 1 1];
for ii=1:length(i)
H(i(ii)).XYZ(j(ii),:)=[];H(i(ii)).sgmXYZ(j(ii),:)=[];
H(i(ii)).NEU(j(ii),:)=[];H(i(ii)).AZEL(j(ii),:)=[];H(i(ii)).VXYZ(j(ii))=[];
H(i(ii)).VNEU(j(ii))=[];H(i(ii)).VAZEL(j(ii))=[];
end
%
%i=[5 11 12];j=[2 1 1];
%i=[11];j=[1];

% H(i).XYZ(j,:)=[];H(i).sgmXYZ(j,:)=[];
% H(i).NEU(j,:)=[];H(i).AZEL(j,:)=[];H(i).VXYZ(j)=[];
% H(i).VNEU(j)=[];H(i).VAZEL(j)=[];
%% 加载碗形分布的散点数据
flnms=dir('S/');
flnms(1:2)=[];
ks=0;
ELs=[88];
%AZs=[45:-15:-315];%+7.5;
AZs=[45:-15:-315];
%a0s=[429.5*1.31 426.4*1.1 388.1*1.15 344.0*1.1 312.0*1.23]*10/340;
%b0s=[0.1 0 0 0 0.2]*1.0;%*ones(1,5);
a0s=[429.5*1.21 426.4*1.12 388.1*1.09 344.0*1.05 312.0*0.96]*10/340;
b0s=[0.1 0 0 0 0.2]*1.0;%*ones(1,5);
ks=0;
for i=1:length(flnms)
    if strcmp(flnms(i).name(1),'S')
        ks=ks+1;
        TEMP=load(['S/' flnms(i).name]);
         S(ks).XYZ=TEMP(:,1:3);%+nsi*randn(size(TEMP,1),3
         %H(ks).XYZ=TEMP(:,1:3)-nsi*ones(size(TEMP,1),3);
         S(ks).sgmXYZ=TEMP(:,4:6);
         %------------------
         S(ks).NEU=S(ks).XYZ;
         for j=1:size(S(ks).sgmXYZ,1)
         S(ks).VXYZ(j).VXYZ=diag(S(ks).sgmXYZ(j,:).^2+nsi.^2);
         S(ks).VNEU(j).VNEU=S(ks).VXYZ(j).VXYZ;
         S(ks).VAZEL(j).VAZEL=diag([((5e-3)*pi/180).^2+(nsi/norm([a0s(ks),b0s(ks)])).^2,((5e-3)*pi/180).^2+(nsi/norm([a0s(ks),b0s(ks)])).^2]);%<====== 天线方位俯仰指向设计精度0.005度
         end
         S(ks).AZEL=[AZs' ELs*ones(25,1)];
    end
end
i=[1 1 2 2 3 3 4 4 5 5]  ;j=[1 24 1 24 1 24 1 24 1 24];
for ii=1:length(i)
S(i(ii)).XYZ(j(ii),:)=[];S(i(ii)).sgmXYZ(j(ii),:)=[];
S(i(ii)).NEU(j(ii),:)=[];S(i(ii)).AZEL(j(ii),:)=[];S(i(ii)).VXYZ(j(ii))=[];
S(i(ii)).VNEU(j(ii))=[];S(i(ii)).VAZEL(j(ii))=[];
end




%% 首次解算用到的X的初值
OA=38.15*pi/180;OE=-56.49*pi/180;
a=15.23;%a=15.30;
e=0.001;
XPR=63.3648;YPR=-12.394;ZPR=8.67867;
af=0.005*pi/180;bt=0.005*pi/180;gm=0.005*pi/180;
b0=2.51; %2.14
XoldG=[XPR,YPR,ZPR,e,af,bt,gm,OA];
XoldL=[];
for i=1:length(H)
    if mod(i,2) == 1;b=-b0;%-2.14;
    elseif mod(i,2) == 0; b=b0;%2.14;
    end
    XoldL=[XoldL,a,b,OE];
end
Xold=[XoldG,XoldL]';
%%
%H(20)=[];
k=0;
L0=[];F=[];Fc=[];
for i=1:length(H)
    for j=1:size(H(i).sgmXYZ,1)
        X=H(i).XYZ(j,1);Y=H(i).XYZ(j,2);Z=H(i).XYZ(j,3);
        AZ=H(i).AZEL(j,1)*pi/180;EL=H(i).AZEL(j,2)*pi/180;
        %
        k=k+1;
        if mod(i,2) == 1
        b=-b0;%-2.14;
        elseif mod(i,2) == 0
        b=b0;%2.3;%b=2.14; 
        end        
        %F的方程结果，计算闭合差(注意加负号)，直接代入到F中，此处加一个子程序getF
        Fc(k).Fc=getFtc(X,Y,Z,XPR,YPR,ZPR,af,bt,gm,e,AZ,EL,OA,OE,a,b);
        %%xia Atargetc
        %pFpa=[];pFpb=[];pFpOE=[];
        [pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA,pFpa,pFpb,pFpOE,pFpXYZ,pFpAZ,pFpEL]=...
            ipdtc(X,Y,Z,AZ,EL,XPR,YPR,ZPR,e,af,bt,gm,OA,a,b,OE);      
        %组建矩阵Bc,表示pFpL
        %组建矩阵Bc,表示pFpL
        Bc(k).Bc=[pFpXYZ,pFpAZ,pFpEL];
        %组建矩阵Aconstc
        Aconstc(k).Aconstc=[pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA];
        %组建矩阵Atargetc
        Atargetc(i).Atargetc(j).At=[pFpa,pFpb,pFpOE];
        %组建XYZ,AZ,EL的初值矢量
        L0=[L0,X,Y,Z,AZ,EL];
    end
end

%%


% OAs=38.15*pi/180;
% OEs=-asin([5.6/6.5  5.2/6 4.8/5.3 4.2/4.5 3.6/3.7]);
% % OEs=asin(6/6.5);
% % OEs=asin(5.3/6.5);
% % OEs=asin(4.5/6.5);
% % OEs=asin(3.7/6.5);
% %a0s=[429.5 426.4 388.1 344.0 312.0]*10/340;
% a0s=[429.5 426.4 388.1 344.0 312.0]*10/340;
% b0s=-[0.2 0.2 0.2 0.2 0.2];%*ones(1,5);




%modified values

OAs=38.15*pi/180;
OEs=-asin([5.6/6.5  5.2/6 4.8/5.3 4.2/4.5 3.6/3.7]);
% OEs=asin(6/6.5);
% OEs=asin(5.3/6.5);
% OEs=asin(4.5/6.5);
% OEs=asin(3.7/6.5);
%a0s=[429.5 426.4 388.1 344.0 312.0]*10/340;


XoldLs=[];
for i=1:length(S)
    XoldLs=[XoldLs,a0s(i),b0s(i),OEs(i)];

end
Xold=[Xold' XoldLs]';
%a0s(1)=a0s(1)-1;
%
%k=0;
%L0=[];F=[];Fc=[];
for i=1:length(S)
    for j=1:size(S(i).sgmXYZ,1)
        X=S(i).XYZ(j,1);Y=S(i).XYZ(j,2);Z=S(i).XYZ(j,3);
        AZ=S(i).AZEL(j,1)*pi/180;EL=S(i).AZEL(j,2)*pi/180;
        %       
         k=k+1;
%         if mod(i,2) == 1
%         b=-b0;%-2.14;
%         elseif mod(i,2) == 0
        %b=b0s;%2.3;%b=2.14; 
       % end        
        %F的方程结果，计算闭合差(注意加负号)，直接代入到F中，此处加一个子程序getF
        Fc(k).Fc=getFtc(X,Y,Z,XPR,YPR,ZPR,af,bt,gm,e,AZ,EL,OA,OEs(i),a0s(i),b0s(i));
        %%xia Atargetc
        %pFpa=[];pFpb=[];pFpOE=[];
        [pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA,pFpa,pFpb,pFpOE,pFpXYZ,pFpAZ,pFpEL]=...
            ipdtc(X,Y,Z,AZ,EL,XPR,YPR,ZPR,e,af,bt,gm,OA,a0s(i),b0s(i),OEs(i));      
        %组建矩阵Bc,表示pFpL
        %组建矩阵Bc,表示pFpL
        Bc(k).Bc=[pFpXYZ,pFpAZ,pFpEL];
        %组建矩阵Aconstc
        Aconstc(k).Aconstc=[pFpXYZPR,pFpe,pFpaf,pFpbt,pFpgm,pFpOA];
        %组建矩阵Atargetc
        Atargetcs(i).Atargetc(j).At=[pFpa,pFpb,pFpOE];
        %组建XYZ,AZ,EL的初值矢量
        L0=[L0,X,Y,Z,AZ,EL];
    end
end
%利用Fc求出闭合差Omg
Omg=[];
for i=1:length(Fc)
   Omg=[Omg;Fc(i).Fc]; 
end
Omg=-Omg;%已为Omg加了符号

%从Omg中找靶标点在XYZ方向上的变化规律，用以修正初值
k1=0;k2=0;k3=0;
for i=1:length(Omg)
    if mod(i,3)==1
        k1=k1+1;
        Xf(k1)=Omg(i);
    elseif  mod(i,3)==2
        k2=k2+1;
        Yf(k2)=Omg(i);
    elseif  mod(i,3)==0
        k3=k3+1;
        Zf(k3)=Omg(i);    
    end
end
% figure
% plot(Xf,'ro');hold on;
% plot(Yf,'go');hold on;
% plot(Zf,'bo');hold off;
% xlabel('观测靶标位置数目')
% ylabel('单位：米')
% title('靶标点拟前残差')
% %disp(num2str(sqrt(Omg.^2/length(Omg))));
% %bufuplt=norm(Omg,2);
 bufuplt=sqrt(sum(Omg(:).^2)/length(Omg));
 presgm=bufuplt;
% disp(num2str(bufuplt));
%求QLL
k=0;
for i=1:length(H)
    for j=1:size(H(i).sgmXYZ,1)
        k=k+1;
        QLLc(k).QLLc=blkdiag(H(i).VNEU(j).VNEU,H(i).VAZEL(j).VAZEL);
    end
end
%add S
for i=1:length(S)
    for j=1:size(S(i).sgmXYZ,1)
        k=k+1;
        QLLc(k).QLLc=blkdiag(S(i).VNEU(j).VNEU,S(i).VAZEL(j).VAZEL);
    end
end



QLL=[];
for i=1:length(QLLc)
    QLL=blkdiag(QLL,QLLc(i).QLLc);
end

qLLS=diag(QLL);

%for i=1:length(qLLS)
% qLLS0=find(qLLS~=qLLS(4));
% QLLS=diag(qLLS0);
% 
% disp('aaaaa')
% persgm1=sqrt(Omg'*inv(QLLS)*Omg/sum(diag(inv(QLLS))));
% disp(num2str(persgm1));
% fprintf('%.1f  \n',persgm1*1e3)

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


Ats=[];
for i=1:length(Atargetcs)
    At00=[];
   for j=1:length(Atargetcs(i).Atargetc)       
       At00=vertcat(At00,Atargetcs(i).Atargetc(j).At) ; 
   end
   At0s(i).At=At00;
end

Ats=[];
for i=1:length(At0s)
    Ats=blkdiag(Ats,At0s(i).At);
end

Att=blkdiag(At,Ats);
A=[Ac,Att];

%组建矩阵B,将其Bc所有子阵组合起来
B=[];
for i=1:length(Bc)
    B=blkdiag(B,Bc(i).Bc);
end

Acurve=A'*inv(B*QLL*B')*A;
bcurve=-A'*inv(B*QLL*B')*Omg;
xhat=Acurve\bcurve;
khat=inv(B*QLL*B')*(-Omg-A*xhat);
ehat=QLL*B'*inv(B*QLL*B')*(-Omg-A*xhat);

sgmhat2=-khat'*Omg/(size(A,1)-size(A,2));  %<=======================
L0=L0';


%从Omg中找靶标点在XYZ方向上的变化规律，用以修正初值
clear Xf Yf Zf
k1=0;k2=0;k3=0;
for i=1:length(ehat)
    if mod(i,5)==1
        k1=k1+1;
        Xf(k1)=ehat(i);
    elseif  mod(i,5)==2
         k2=k2+1;
        Yf(k2)=ehat(i);
    elseif  mod(i,5)==3
        k3=k3+1;
        Zf(k3)=ehat(i);   
        
    end
end
% figure
% plot(Xf,'ro');hold on;
% plot(Yf,'go');hold on;
% plot(Zf,'bo');hold off;
% xlabel('观测靶标位置数目')
% ylabel('单位：米')
% title('靶标点一次迭代拟后残差')
LL0=L0;
bufuplt=sqrt(sum(ehat(:).^2)/length(ehat));
disp(num2str(bufuplt));
save LL0 LL0 QLL








