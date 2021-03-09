clc;clear all;close all;
%% 对方向图公式进行变量初始化
F=12.5e9;   %工作频率
lambda=physconst('lightspeed')/F;   %工作波长
e=2.9;  %微带结构的等效介电常数εe
k=2*pi/lambda;   %相位常数
ks=2*pi*sqrt(e)/lambda;  %介质中的相位常数
l=0.004;    %贴片宽度
p=0.005;    %单元周期间隔
N=5;%阵列单元的个数

%% 对方向图进行采样
Fs=3600;
S1=1:1:Fs;%-pi-pi一共采了W个点
thta=-pi+pi/(Fs/2):pi/(Fs/2):pi; % -pi to pi

%% 对阵列单元赋初值并产生方向图
M=5;%157;%训练数据和测试数据一共M组
rng(0);
m=round(rand(M,N)*1);%产生M*N的0/1的随机数作为阵列的输入
% ni=0:M-1;
% for i=1:N
%    m(:,N+1-i)=mod(bitshift(ni(:),-i+1),2);
% end
% m(1,:)=1

f1=cos(ks*l*cos(thta)/2);   %阵元方向性函数
f2=zeros(M,Fs);f3=zeros(M,Fs);F=zeros(M,Fs);   
for i=1:M
    n=zeros(1,Fs);
    for h=1:N
        n=n+m(i,h)*exp((-1i*(h-1)*(k*p*sin(thta)-ks*p)));
    end
    f2(i,:)=abs(n);  %阵因子
    f3(i,:)=(f1.*f2(i,:));  %总方向图
    max_f=max(f3(i,:));   
    F(i,:)=20*log10(f3(i,:)/max_f);
end
%% 频率
Fs = 1024*4;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
t = 0:1/Fs:2*pi;        % Time vector
FF=fft(F(1,:));

N=40; %order of approximation
poles=-2*pi*logspace(0,4,N); %Initial poles
weight=ones(1,Fs); 
figure(2);
plot(S1,F(1,:),S1, FF);

%
%对向量拟合技术设置迭代条件
def.relax=1;      %Use vector fitting with relaxed non-triviality constraint 
def.stable=0;     %Enforce stable poles 
def.asymp=3;      %Include only D in fitting (not E)    
def.skip_pole=0;  %Do NOT skip pole identification 
def.skip_res=0;   %Do NOT skip identification of residues (C,D,E)  
def.cmplx_ss=1;   %Create complex state space model 
def.spy1=0;       %No plotting for first stage of vector fitting 
def.spy2=1;       %Create magnitude plot for fitting of f(s)  
def.logx=1;       %Use logarithmic abscissa axis 
def.logy=0;       %Use logarithmic ordinate axis  
def.errplot=1;    %Include deviation in magnitude plot 
def.phaseplot=0;  %exclude plot of phase angle (in addition to magnitiude) 
def.legend=1;     %Do include legends in plots 
opts=def;

S2=80001:1:83600;%F2为方向图新的频率区间
%}

s=S1*1i; %复数频率

%% 向量拟合

%[SER,poles,rmserr,fit,opts]=vectfit3_1(F(1,:),s,poles,weight,opts);
[SER,poles,rmserr,fit,opts]=vectfit3_1(FF,s,poles,weight,opts);

%% 恢复
A=SER.A;
B=SER.B;
C=SER.C;
D=SER.D;
E=SER.E;
I=eye(N,N);
for i=1:Fs
    fit2(i) = C*(s(i)*I-A)^-1*B +D +s(i)*E;
end
figure(3);
%plot(S1,F(1,:),S1,real(fit),S1,real(fit2));
plot(S1,F(1,:),S1,real(ifft(fit)),S1,real(ifft(fit-fit2)));


%figure(2);
%plot(s,F(1,:));
%hold on;


