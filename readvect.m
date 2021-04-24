function [order,A,C,D,E,frcv]=readvect(filepath)
%================================================
% Order %阶数
% Re_An %极点A-实数部分-奇数组
% Im_An %极点A-虚数部分-奇数组
% Re_Cn %留数C-实数部分-奇数组
% Im_Cn %留数C-虚数部分-奇数组
% D     %D-实数部分
% E     %E-实数部分
%================================================
way=2;%写入方式：1-fprint；2-csvwrite
Nsmp=3600;
offset=8e4;
S=1j*((1:Nsmp)+offset);
if way==1
    ifnan=1;%尾列为NaN
elseif way==2
    ifnan=0;
end

RecoverData=readmatrix(filepath,'Delimiter',',');
[nrow,ncol]=size(RecoverData);

Ndata=nrow;
order=RecoverData(:,1);
D=RecoverData(:,end-1-ifnan);
E=RecoverData(:,end-ifnan)*1e-6;

MaxOrder=(ncol-3-ifnan)/2;
RemainData=RecoverData(:,2:end-2-ifnan);

Re_An=RemainData(:,1:MaxOrder/2);
Im_An=RemainData(:,1+MaxOrder/2:MaxOrder);
Re_Cn=RemainData(:,1+MaxOrder:MaxOrder*3/2);
Im_Cn=RemainData(:,1+MaxOrder*3/2:MaxOrder*2);
A(:,1:2:MaxOrder)=Re_An+1j*Im_An;
A(:,2:2:MaxOrder)=Re_An-1j*Im_An;
C(:,1:2:MaxOrder)=Re_Cn+1j*Im_Cn;
C(:,2:2:MaxOrder)=Re_Cn-1j*Im_Cn;

for i=1:Ndata
    for ii=1:Nsmp
        sub=S(ii)-A(i,:);%广播
        frcv(i,ii)=sum(C(i,:)./(sub))+D(i)+S(ii)*E(i);%由output恢复的拟合函数,
    end
    size(frcv);
end
%由于忽略了B和浮点储存的原因，恢复后会增加1e-26级别的误差。
frcv;
end