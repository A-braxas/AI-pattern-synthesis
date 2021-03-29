%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 在（1_2）版本中，对于阵元数为Nunit的阵列天线，随机生成Ndata组幅值序列
% 和Ndata组样本数据，每组进行Nsmp的采样，进行了阶数为Norder的vecterfit拟合
% 拟合进行了P次迭代，得到了[A,C,D,E,ff]等参数，并通过readvect复原
% ---------------------------
% 在（1_3）版本中，将完成以下目标：
% 1.阶数Norder变成Ndata*1维度（✓）
% 2.对out进行文件读写(✓)
% 3.对阶数从10：2：50进行扫描，每次误差用数组储存，对误差画图分析，找出误差最小值对应的最小阶数（✓）
% 4.建立一个Ndata*1的[阶数]数组，其后再用此阶数对函数重新拟合，并记录极点留数（✓）
% 5。对迭代次数也可以进行扫描分析（×）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 流程：
% 1.天线变量初始化
% 2.生成方向图
% 3.vector fitting参数设置
% 4.扫描阶数
% 5.作阶数-误差关系图，找最小误差阶数，按最小阶数重新生成极点留数
% 6.复数按虚实和奇偶拆分，并输出极点留数储存矩阵
% 7.读取储存矩阵并复原测试
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
Nunit=5;    %阵列单元的个数

%% 采样方向图
Nsmp=3600*1;  %采样点数
S1=1:1:Nsmp;   %采样序列
thta=-pi+2*pi/Nsmp:2*pi/Nsmp:pi; % -pi to pi

%% 产生阵列方向图
Ndata=200;%训练数据和测试数据一共Ndata组
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
Ntest=200;    % Ntest < Ndata
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

%设置初始阶数，因极点和留数成共轭对出现，需为偶数
MinOrder=10;%起点
Internal=2;%间距
MaxOrder=50;%终点
Len=floor((MaxOrder-MinOrder)/Internal)+1;%扫描次数
Norder=MinOrder:Internal:MaxOrder;%扫描节点
%设置误差存储矩阵
Err=zeros(Ntest,Len); %均方误差
%% 阶数扫描
for i=1:Ntest
    for ii=1:Len
        Q=Norder(ii);%临时阶数
        %设置最佳初始极点
        startpoles=InitPoles(Q,offset,Nsmp);
        %向量拟合与迭代计算
        P=1;%P为迭代次数
        %[a2 cde2 rmserr fit]=[极点A 留数C和DE 均方根误差 拟合函数]
        [t1,t2,t3,t4,rmserr,fit]=vectfit3(Ftest(i,:),S2,startpoles,weight,opts);
        for iii=1:P-1
            [t1,t2,t3,t4,rmserr,fit]=vectfit3(Ftest(i,:),S2,t1,weight,opts);
        end
        Err(i,ii)=rmserr;%!!
    end
end
%% 作阶数-误差关系图
if(0)
figure(1);
for i=1:Ntest
    S3=MinOrder:Internal:MaxOrder;
    plot(S3,log(Err(i,:)));
    hold on;
end
end
%记录最小误差下的阶数
OrderList=zeros(Ntest,1);%阶数列表
[minErr,index]=min(Err');
for i=1:Ntest
    OrderList(i)=Norder(index(i));
end

%% 按最小阶数重新生成极点留数
%设置极点、留数、拟合函数的存储矩阵
MaxOrderList=max(OrderList);
An=zeros(Ntest,MaxOrderList);%An和Cn分别为极点和留数矩阵,
Cn=zeros(Ntest,MaxOrderList);%
D=zeros(Ntest,1);
E=zeros(Ntest,1);
Err2=zeros(Ntest,1);
Ffit=zeros(Ntest,Nsmp);
for i=1:Ntest
    Len2=OrderList(i);
    startpoles=InitPoles(Len2,offset,Nsmp);
    [An(i,1:Len2),Cn(i,1:Len2),D(i),E(i),Err2(i),Ffit(i,:)]=vectfit3(Ftest(i,:),S2,startpoles,weight,opts);
end
%% 复数按虚实和奇偶拆分
if(1)
%{
A=Re_An±j*Im_An
C=Re_Cn±j*Im_Cn
D=D
E=E
%}
Re_An=real(An(:,1:2:end));Im_An=imag(An(:,1:2:end));
Re_Cn=real(Cn(:,1:2:end));Im_Cn=imag(Cn(:,1:2:end));
OUT=[OrderList Re_An Im_An Re_Cn Im_Cn D E*1e6];    % !!!E*1e6,注意恢复
end
%% 输出
if(0)
filepath='SavedDataAsFormOf_order_An_Cn_D_E.txt';

fileID = fopen(filepath,'w');
fprintf(fileID,[repmat('%f,', 1, size(OUT,2)), '\n'], OUT');
fprintf('/n');
fclose(fileID);
elseif(1) %% 按.csv格式保存数据
    writematrix(OUT,'Labels_order_An_Cn_D_E.csv');
    writematrix(m,'Features_m.csv');
    writematrix(Ffit,'Ffit.csv');
end

%% 极点和留数复原
if(0)
%对k需要前置处理
[order,An,Cn,D,E,frcv]=readvect(filepath);
err=SqrtError(Ffit,frcv)
Err2
elseif(1)
[order2,An2,Cn2,D2,E2,frcv]=readvect('Labels_order_An_Cn_D_E.csv');
err=SqrtError(Ffit,frcv)
Err2;
end

%% 自定义函数
% 初始化极点矩阵
function startpoles=InitPoles(Q,offset,Nsmp)
    startpoles=zeros(1,Q);  %startpoles为初始极点矩阵
    for iiii=1:2:Q
        beta=offset+1+Nsmp*(iiii-1)/Q;
        alpha=(beta-1)/100;
        startpoles(iiii)=-alpha+1i*beta;
        startpoles(iiii+1)=-alpha-1i*beta;
    end
end
% 求均方误差函数
function err=SqrtError(f1,f2)
err=sqrt(sum(abs((f1-f2).^2),2)/size(f1,2));
end


