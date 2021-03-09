function [frcv,smp,A,C,D,E]=readvect(k0, k1, k2, k3, k4, k5, k6, k7, k8);
% k0=Norder;%阶数
% k1=r_an(:,1:2:end); %极点A-实数部分-奇数组
% k2=i_an1(:,1:2:end); %极点A-（虚数部分绝对值-80000）-奇数组
% k3=r_c(:,1:2:end); %留数C-实数部分-奇数组
% k4=i_c(:,1:2:end);% 留数C-虚数部分-奇数组
% k5=r_cn(:,Norder+1); %D-实数部分
% k6=(1e6)*r_cn(:,Norder+2);%1e6*E-实数部分
% k7=Nsmp;%采样点数
% k8=offset;%S的位移量
smp=k7;
offset=k8;
A(1:2:k0)=k1+1j*(k2+offset);
A(2:2:k0)=k1-1j*(k2+offset);
C(1:2:k0)=k3+1j*k4;
C(2:2:k0)=k3-1j*k4;
D=k5;
E=k6*1e-6;
S=1j*((1:k7)+offset);
for i=1:smp
    sub=S(i)-A;%广播
    frcv(i)=sum(C./(sub))+D+S(i)*E;%由output恢复的拟合函数,
end
%由于忽略了B和浮点储存的原因，恢复后会增加1e-26级别的误差。
end