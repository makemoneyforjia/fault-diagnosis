%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%      WPD-ASSTFT        %%%%%%%%%%%%%
%%%%%%%%%%%%%         CSQ          %%%%%%%%%%%%%
%%%%%%%%%%%%%         成都           %%%%%%%%%%%%%
%%%%%%%%%%%%%       202303      %%%%%%%%%%%%%
%%%%%%%%%%%%%        version 1     %%%%%%%%%%%%%%%
%%小波包分解：参开文献：基于轴箱振动加速度的钢轨波磨评价方法及应用_徐晓迪
clear 
close all

%% 载入数据
disp('载入数据');
Data=load('D:\桌面\WPD-ASSTFT\Data9(10m)\Data9(28555.187-28545.187).txt');%载入信号（列1时间，列2加速度，列3公里标）
% Velocity=load('D:\桌面\WPD-ASSTFT\Data\Data9_speed.txt');%%速度

data.time=Data(:,1);%时间
data.acc=Data(:,2);%轴箱垂向加速度
data.distance=Data(:,3);%路程即公里标
data.line=size(Data,1);%采样点个数
data.V=Data(:,4);%速度
data.Vm=abs((data.distance(end)-data.distance(1))/(data.time(end)-data.time(1)));%%求平均速度
offset=2000/data.Vm;%%每1m对应的点的个数
fs=2000;

%% 画原始图像
disp('画原始图像');
figure;
set(gcf,'Position',[100,100,600,400]);
subplot(211)
hold on


xlabel('\fontsize{10}\fontname{宋体}里程 \fontsize{10}\fontname{Times New Roma}/m');	%设置x轴字体
ylabel('\fontsize{10}\fontname{宋体}加速度')
plot(data.distance,data.acc);

yyaxis right
plot(data.distance,data.V,'color','r');
set(gca,'ycolor','k');
ylabel('\fontsize{10}\fontname{宋体}速度 \fontname{Times}km/h');
if data.distance(end)<data.distance(1)
    set(gca,'xdir','reverse');
end

plotAmplitude(data.acc,'原始信号',0.5,1000,fs)
% figure;[Pxx,f] = pwelch(Data(:,2),[],[],[],fs);
% hold on
% plot(f(10:end),Pxx(10:end));
% [Max,pos]=max(Pxx(10:end));
% plot(f(pos,1),Max,'k-o','MarkerFacecolor','k');
% text(f(pos,1)+0.02*(max(xlim)-min(xlim)),Max,[num2str(round(f(pos,1),1)),'Hz'],'color','k','fontsize',6);
% xlabel('frequency(hz)');
% ylabel('power spectral density(db/hz)');
% title('period psd estimate');
hold off
%% 数据前处理,以路程为parameter，将数据分段到两站之间
%Data_preprocessing1 (data);

%进一步分段

%% 小波分解  修改小波类型看效果
wpt=wpdec(data.acc,3,'dmey','shannon') %wpdec(信号源，分解层数，小波类型，采用熵类型。）

plot(wpt);
% cfs_1=wpcoef(wpt,[3,0]);    %%获取小波树上某个节点的小波包系数
% cfs_2=wpcoef(wpt,[3,1]);
% cfs_3=wpcoef(wpt,[3,2]);
% cfs_4=wpcoef(wpt,[3,3]);
% cfs_5=wpcoef(wpt,[3,4]);
% cfs_6=wpcoef(wpt,[3,5]);
% cfs_7=wpcoef(wpt,[3,6]);
% cfs_8=wpcoef(wpt,[3,7]);
%%Reconstruct wavelet packet coefficients.
nodes=[7;8;9;10;11;12;13;14]; %第3层的节点号
ord=wpfrqord(nodes);  %小波包系数重排，ord是重排后小波包系数索引构成的矩阵　如3层分解的[1;2;4;3;7;8;6;5]
nodes_ord=nodes(ord); %重排后的小波系数
%% 绘制小波分解后的子信号和FFT变换幅频图
disp('小波分解后子信号&FFT变换幅频图');
for i=1:8
    rex(:,i)=wprcoef(wpt,nodes_ord(i));
    figure;
    hold on
    set(gcf,'Position',[100,100,800,600]);
    subplot(2,1,1);
    plot(data.distance,rex(:,i));
    xlabel({'\fontsize{10}\fontname{宋体}里程 \fontsize{10}\fontname{Times New Roma}/m',['rex',num2str(i)]});	%设置x轴字体
    ylabel('\fontsize{10}\fontname{宋体}加速度')
  

    if data.distance(end)<data.distance(1)
        set(gca,'xdir','reverse');
    end
    plotAmplitude(rex(:,i),['rex',num2str(i)],2,1000,fs)
    
    hold off
end
%% 选择最优窗长 ASSTFT
disp('最优窗长，瞬时频率');
database=cell(1,10);
databaseS=cell(8,2);
TFS=0;
TS=0;
for i =1:8
    Optimal_windowlength(rex(:,i),2000,data.Vm,offset);
    [b,Lmax,Lmin,Ampmax,fmax,s,f,t,p,window_lenth]=Optimal_windowlength(rex(:,i),2000,data.Vm,offset);
    database(i,:)={b,Lmax,Lmin,Ampmax,fmax,s,f,t,p,window_lenth};%%储存数据
    figure;
    instfreq(rex(:,i),fs)
    [itf,itt]=instfreq(rex(:,i),fs);
    ifq=instfreq(abs(s),f,t);
%%   同步压缩  

  [tfr,Ts] =SST2(rex(:,i),fs,window_lenth);%每30m一段
%   databaseS(i,:)={tfr,Ts};
  TFS=TFS+tfr;
  TS=TS+Ts;
  title(['rex',num2str(i)]);
%% 最终图像
figure;[y x]=size(TFS)
subplot(211)
imagesc(0:1/fs:1/fs*x,fs/y:fs/y:fs/2,abs(TFS));
axis xy
ylabel('Freq / Hz');
xlabel('Time / Sec');
title('PSTFT');
colormap('jet');
ch1=colorbar,axis tight;
subplot(212)
imagesc(0:1/fs:1/fs*x,fs/y:fs/y:fs/2,abs(TS));
axis xy
ylabel('Freq / Hz');
xlabel('Time / Sec');
title('SST');
colormap('jet');
ch2=colorbar,axis tight;
%     figure;
%     [s,w,n] = fsst(ifq);
% 
%     mesh(n,w/pi,abs(s))
%     
%     axis tight
%     view(2)
%     colorbar
%   figure;
%   fsst(rex(:,i),fs,'yaxis');
%   [fsf(:,i),fst(:,i)] = instfreq(abs(s).^2,f,t);
%   fsst(intf(:,i));
end
%% 瞬时频率
%   for i =1:8
%  figure;
%  set(gca,'Position',[100,100,600,400]);
%  instfreq(rex(:,i),2000);
%  ifq = instfreq(rex,fs); 
% end
%% 同步压缩小波变换





