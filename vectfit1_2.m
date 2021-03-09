%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 在（1_2）版本中，对于阵元数为Nunit的阵列天线，随机生成Ndata组幅值序列
% 和Ndata组样本数据，每组进行Nsmp的采样，进行了阶数为Norder的vecterfit拟合
% 拟合进行了P次迭代，得到了[A,C,D,E,ff]等参数，并通过out储存
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;clear all;close all;
%% 变量初始化
Fw=12.5e9;   %工作频率
lambda=physconst('lightspeed')/Fw;   %工作波长
e=2.9;  %微带结构的等效介电常数εe
k=2*pi/lambda;   %相位常数
ks=2*pi*sqrt(e)/lambda;  %介质中的相位常数
l=0.004;    %贴片宽度
p=0.005;    %单元周期间隔
Nunit=1;    %阵列单元的个数

%% 采样方向图
Nsmp=3600*1;  %采样点数
S1=1:1:Nsmp;   %采样序列
thta=-pi+2*pi/Nsmp:2*pi/Nsmp:pi; % -pi to pi

%% 产生阵列方向图
Ndata=157;%训练数据和测试数据一共Ndata组
rng(0);
m=round(rand(Ndata,Nunit)*1);%产生Ndata*Nunit的0/1的随机数作为阵列的输入
m(sum(m,2)==0,1)=1;%确保每组幅值至少一个不为0
    
Funit=cos(ks*l*cos(thta)/2);    %阵元方向性函数 Funit
Farray=zeros(Ndata,Nsmp);       %阵因子 Farry
Ftotal=zeros(Ndata,Nsmp);       %总方向图 Ftotal=Funit*Farry
%Fnorm=zeros(Ndata,Nsmp);        %归一化方向图 Ftotal/max{Ftotal}
for i=1:Ndata
    n=zeros(Nunit,Nsmp);
    for h=1:Nunit
        n(h,:)=m(i,h)*exp((-1i*(h-1)*(k*p*sin(thta)-ks*p)));
    end
    Farray(i,:)=sum(n);  %阵因子？？
    Ftotal(i,:)=(Funit.*Farray(i,:));  %总方向图
    %max_f=max(Funit.*abs(Farray(i,:)));   
    %Fnorm(i,:)=20*log10(Funit.*abs(Farray(i,:))/max_f);
end
%% Vectfit参数设置 
%取测试数据
Ntest=157;    % Ntest < Ndata
Ftest=Ftotal(1:Ntest,:);    

%由于vectfit是在复频域进行拟合，需要把坐标轴挪到远离低频的复数域
offset=8e4;
S2=1i*(S1+offset);

%
%对向量拟合技术设置条件
def.relax=1;      %Use vector fitting with relaxed non-triviality constraint 
def.stable=0;     %Enforce stable poles 
def.asymp=3;      % Include only D in fitting (not E)    
def.skip_pole=0;  % Do NOT skip pole identification 
def.skip_res=0;   % Do NOT skip identification of residues (C,D,E)  
def.cmplx_ss=1;   %Create complex state space model 
def.spy1=0;       % No plotting for first stage of vector fitting 
def.spy2=0;       % Create magnitude plot for fitting of f(s)  
def.logx=1;       % Use logarithmic abscissa axis 
def.logy=0;       % Use logarithmic ordinate axis  
def.errplot=1;    % Include deviation in magnitude plot 
def.phaseplot=0;  % exclude plot of phase angle (in addition to magnitiude) 
def.legend=1;     % Do include legends in plots 
opts=def;
%}

weight=ones(1,Nsmp);    %设置权重
Norder=20;              %设置阶数，因极点和留数成共轭对出现，需为偶数
%设置最佳初始极点
initpoles=zeros(1,Norder);  %initpoles为初始极点矩阵
for k=1:2:Norder
    beta=offset+1+Nsmp*(k-1)/Norder;
    alpha=(beta-1)/100;
    initpoles(k)=-alpha+1i*beta;
    initpoles(k+1)=-alpha-1i*beta;
end
%% 向量拟合与迭代计算
P=1;%P为迭代次数
An=zeros(Ntest,Norder);
Cn=zeros(Ntest,Norder+2);%An和Cn分别为极点和留数矩阵,'+2'为D，E
err=zeros(Ntest,1); %均方误差
FF=zeros(Ntest,Nsmp);%拟合函数
for i=1:Ntest
    a2=zeros(P,Norder);
    c2=zeros(P,Norder+2);
    [a2(1,:),c2(1,:),rmserr,fit]=vectfit3(Ftest(i,:),S2,initpoles,weight,opts);
    for j=1:P-1
        [a2(j+1,:),c2(j+1,:),rmserr,fit]=vectfit3(Ftest(i,:),S2,a2(j,:),weight,opts);
    end
    An(i,:)=a2(P,:);%??P-1
    Cn(i,:)=c2(P,:);
    err(i)=rmserr;
    FF(i,:)=fit;
end
%% 极点和留数分组储存
r_an=real(An);i_an=imag(An);i_an1=abs(i_an)-offset; 
r_cn=real(Cn);i_cn=imag(Cn);
r_c=r_cn(:,1:Norder);
i_c=i_cn(:,1:Norder);

k=Norder;%阶数
k1=r_an(:,1:2:end); %极点A-实数部分-奇数组
k2=i_an1(:,1:2:end); %极点A-（虚数部分绝对值-80000）-奇数组
k3=r_c(:,1:2:end); %留数C-实数部分-奇数组
k4=i_c(:,1:2:end);% 留数C-虚数部分-奇数组
k5=r_cn(:,Norder+1); %D-实数部分
k6=(1e6)*r_cn(:,Norder+2);%1e6*E-实数部分
k7=Nsmp;%采样点数
k8=offset;%S的位移量
for i=1:Ndata
    
    out=[k k1(i) k2(i) k3(i) k4(i) k5(i) k6(i) k7 k8];
    %out=[k k1 k2 k3 k4 k5 k6 k7 k8];
end
%{
A=k1±j*(k2+offset)
C=k3±j*k4
D=k5
E=k6*1e-6
%}

%% 极点和留数复原
% 对k需要前置处理
% [ff,smp,A,C,D,E]=readvect(k, k1, k2, k3, k4, k5, k6, k7, k8)
% err=sqrterr(fit,ff)

%% 求均方误差函数
function err=sqrterr(f1,f2);
err=sum(abs(f1-f2).^2);
end


