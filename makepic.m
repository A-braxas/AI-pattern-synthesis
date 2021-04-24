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
Nunit=10;    %阵列单元的个数
Nscan=1;
Ndata=Nscan^Nunit;  %训练数据和测试数据一共Ndata组
Nsmp=3600*1;  %采样点数
Ntest=floor(Ndata*1);    % Ntest < Ndata

smp=round(logspace(1,3,20))*50;
smp_err=zeros(length(smp),1);
for jj=1:length(smp)
    Nsmp = smp(jj);
%% 采样方向图
S1=1:1:Nsmp;   %采样序列
thta=-pi+2*pi/Nsmp:2*pi/Nsmp:pi; % -pi to pi

%% 产生阵列方向图
if(1)
     m=[0.8147, 0.9058, 0.1270, 0.9134, 0.6324, ...
         0.0975, 0.2785, 0.5469, 0.9575, 0.9649];
elseif(0)   %(1_5) %产生Ndata*Nunit的均匀分布数作为阵列的输入
    m=zeros(Ndata,Nunit);
    cols=zeros(Nscan,Nunit);
    rng(1);
    col=1/2/Nscan:1/Nscan:1;
   
    for i=1:Nunit
        cols(:,i)=normrnd(col,0.01);
    end
    % 只得手动控制输入变量
    [j i h g f e d c b a] = ndgrid(cols(:,1),cols(:,2),cols(:,3)...
        ,cols(:,4),cols(:,5),cols(:,6)...
        ,cols(:,7),cols(:,8),cols(:,9)...
        ,cols(:,10));    %cols(:,4),cols(:,5)
    m=[a(:) b(:) c(:) d(:) e(:) f(:) g(:) h(:) i(:) j(:)];     %d(:),e(:)
    m(m>=1)=0.99;
    m(m<=0)=0.01;    %确保每组幅值在0～1之间
end

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
Ftest=Ftotal(1:Ntest,:);

% 释放内存
clear -regexp 'F(total|array|unit)'

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
def.spy2=1;       % Create magnitude plot for fitting of f(s)
def.logx=0;       % Use logarithmic abscissa axis
def.logy=0;       % Use logarithmic ordinate axis
def.errplot=1;    % Include deviation in magnitude plot
def.phaseplot=0;  % exclude plot of phase angle (in addition to magnitiude)
def.legend=1;     % Do include legends in plots
opts=def;
%}
weight=ones(1,Nsmp);    %设置权重

%设置初始阶数，因极点和留数成共轭对出现，需为偶数
MinOrder=2;%起点
Internal=2;%间距
MaxOrder=100;%终点
Len=floor((MaxOrder-MinOrder)/Internal)+1;%扫描次数
Norder=MinOrder:Internal:MaxOrder;%扫描节点
%设置误差存储矩阵
Errthreshold=1e-3; 
Err=zeros(Ntest,Len); %均方误差
OrderList=zeros(Ntest,1);%阶数列表
%% 阶数扫描
for i=1:Ntest
    finishflag = 0;
    Minerr=100000;
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
        Minerr= min(Minerr,rmserr);
        if (rmserr < Errthreshold & ~finishflag)
            OrderList(i,:)=Q;
            break;
            % finishflag = 1;
        elseif ((rmserr <= Minerr) & ~finishflag)
            OrderList(i,:)=Q;
        end
    end
end
smp_err(jj) = Minerr;
end
%% 作采样数-误差关系图
if(1)
     figure(3);
    for i=1:Ntest
        S4=smp';
        plot(S4,10*log10(smp_err),'b+-');
        xlabel('Num of sample'); ylabel('rmserr(dB)'); 
        grid on;
        hold on;
    end
end
%% 作阶数-误差关系图
if(0)
    figure(2);
    for i=1:Ntest
        S3=MinOrder:Internal:MaxOrder;
        S3(Err==Minerr)
        plot(S3,10*log10(Err(i,:)),'b+-');
        xlabel('order'); ylabel('rmserr(dB)'); 
        grid on;
        hold on;
    end
end
%% 由误差确定阶数
if(1)    % (1_5)最小有效阶数，按从小到大排列
    [OrderList,index]=sort(OrderList);
    m=m(index,:);
    %Forigin=Ftest;
    Ftest=Ftest(index,:);
    %sum(sum(Forigin-Ftest))
    
    % 取最大分支
    % MaxBranch = mode(OrderList)
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
