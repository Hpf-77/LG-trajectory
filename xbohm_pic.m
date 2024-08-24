%xbohm理论下画图文件。



%画某一 z 平面的光强图
% zilab = [0,500,2000];
% for i =1
% x = -2.5001:0.01:2.5;
% [x,y] = meshgrid (x,x);
% z =zilab(i);
% ct = real(eval(LGlight_c));
% %cct = conj(ct).*ct;
% cct = ct./max(max(ct));
% 
% cmax = max(max(max(cct)));
% cmin = min(min(min(cct)));
% figure 
% hold on
% A1 = pcolor(x,y,cct);
% caxis([cmin,cmax]);
% set(A1,'edgecolor','none','facecolor','interp')
%  end

%初始数据处理，将直角坐标系转换为直角坐标系
clear pi;
zr = (pi * w0^2) / lamd ;     %瑞利距离

zlab = z0:h:zmax; 
wz = w0 * sqrt(1 + (zlab/zr).^2); 
R = zlab + zr^2 ./ zlab;

 [Vtheta1,Vrlab1] = cart2pol(Vxlab,Vylab);
 [theta1,rlab1] = cart2pol(real(xlab),real(ylab));

%r方向速度
Vr_lab = Vxlab.*cos(theta1) + Vylab.*sin(theta1);
%theta方向  角速度
 Vtheta_lab = -1./rlab1 .* (Vxlab.*sin(theta1) - Vylab.*cos(theta1));

 %r 方向动量
p_theta = (2*pi)/lamd * Vtheta_lab;
% theta 方向动量
p_r = (2*pi)/lamd * Vr_lab;

collab = [1,0,0;0,0,0;0,0,1];
 
%  半径r，角度theta，r方向速度Vr，theta方向 角速度w 随z 传播的变化图
% figure
%  hold on
%  ii=1;
%  for tempi =  [30,50,70]
%     plot(zlab(1:10)./zr, theta1(96,1:10,tempi),'color',collab(ii,:)); %'-*','markerindices',200:1000:10000,
%   %legend_str{1}=['l=3,r_0=0.2',num2str(1)];
%   ii=ii+1;
%  end
% xlabel('z/zr')
%   ylabel('角度theta')
%   title('z-theta图')
% 
%  figure
%  %for fi = 1:16
%      %figure(fi)
%  hold on
%   for tempi2 = 1:6%[3,17,30]
%     plot(zlab(1:10000), rlab1(1,1:10000,tempi2));  %'color',[0.5,0.5,0.5],'-*','markerindices',800:1000:10000,'color','g'
%   end
%   xlabel('z/zr')
%   ylabel('半径r')
%   title('z-r图')

   omg = 1*Vtheta_lab;
figure 
 hold on
 ii=1;
  for tempi2 =  [4,49,92]
    plot(zlab(:), omg(1,:,tempi2),'color',collab(ii,:));%,200:1000:10000,'color'
   
    ii=ii+1;
  end
     %ylim([-4.7e-3 0.3e-3])
     xlim([400  4900])
    xlabel('z/mm')
  ylabel('角速度w/c')
  title('z-w图')
box

  figure 
 hold on
 ii = 1;
  for tempi2 =  [4,49,92]
    plot(zlab(:), Vr_lab(1,:,tempi2),'color',collab(ii,:));   
    ii=ii+1;
  end
  ylim([-0.3e-4 10.3e-4])
  xlim([400  4900])
  xlabel('z/mm')
  ylabel('Vr/c')        
  title('z-Vr图');
  box
  


%  以下为画玻姆轨迹三维空间图代码

modelmax = 9504;                                         %每个平面上的数据取样点数量                   半径数量乘角度数据量   
th =1;                                                              % 角度间隔   6是每隔60°取一个点。    
zh = 2000;                                                       %最远 z
zhmaxi = zh-z0;                                              %(用于索引的最大值)

 
picmodel(picmodel>100 ) = 10*(-15);    %奇点附近值归 0 

%每个点的光强归一化处理
for i = 1:b-1
    lg_c(i) = sum(picmodel(:,i));
end

[Lg_c,~] = meshgrid (lg_c,picmodel(:,1));
Picmodel = picmodel ./ Lg_c;                        %归一化
Picmodel = Picmodel * 3;

%三维散点图                                                     
%某个z平面上，根据玻姆理论求得的坐标模拟位置空间光强分布
figure
zzlab =([500,2000]-z0) ;                                  %待画z平面

for phi = 1:2                                                       % z 平面 索引号
    pz = zzlab(phi)+1;
hold on
%scatter函数将每个位置点对应的强度用颜色和圆圈大小表示，以模拟光强分布。
 scatter3(Xmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Xmodel(1:th:modelmax,pz))), Ymodel(1:th:modelmax,pz),7000*real(Picmodel(1:th:modelmax,pz)),100*real(Picmodel(1:th:modelmax,pz)),'filled')
end

% 三维轨迹图
%xbohm理论下位置空间轨迹
hold on
collab = [1,0,0;0,0,0;0,0,1];
 ii = 1;
for t = [4,49,92]                                                     %  半径 索引号
               for tp =1                                              %  角度 索引号
                    xyzi = tp;
                    hold on
                    plot3(xlab(tp,1:1:zhmaxi ,t),zlab( 1:1:zhmaxi ),ylab(tp,1:1:zhmaxi  ,t),'color',collab(ii,:),'linewidth',0.5);

               end
ii = ii+1;
end


%xbohm理论下 画动量空间轨迹
figure 
hold on
grid on
Pxmodel = zeros(size (Xmodel));
Pymodel  =Pxmodel;
for i = 1:1:b-1                                                     %二维动量轨迹转成一维数组
   for t = 1:rb    
        for tp = 1:1:jmax
             xyzi = (t-1)*jmax+tp;
            Pxmodel(xyzi,i) = pxlab(tp,i,t) ;          
            Pymodel(xyzi,i) =pylab(tp,i,t) ;
        end
   end
end
zzlab =([500,2000]-z0) ; 

%三维散点图                                                     
%模拟某个 z 平面上动量空间内的光强分布。
for phi = 1:2
    pz = zzlab(phi)+1;                                          % z 平面 索引号
hold on
 scatter3(Pxmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Pxmodel(1:th:modelmax,pz))), Pymodel(1:th:modelmax,pz) ,7000*real(Picmodel(1:th:modelmax,pz)),100*real(Picmodel(1:th:modelmax,pz)),'filled')
end

%三维动量空间轨迹图
ii = 1;
for t = [4,49,92]                                                 %  半径 索引号 
               for tp =[96,24,48,72]                          %  角度 索引号
                    xyzi = tp;
                    hold on
                    plot3( real(pxlab(tp,1:1:zhmaxi,t)) , zlab(1:1:zhmaxi) , real(pylab(tp,1:1:zhmaxi,t)) ,'color',collab(ii,:),'linewidth',0.5);
                    scatter3(pxlab(tp,1 ,t),zlab( 1 ),pylab(tp,1,t), 50,collab(ii,:),'*');
               end
ii=ii+1;
end
xlabel('x/mm')
ylabel('z/mm')
zlabel('y/mm')
view(45,45);
box


%画某 z 平面上角速度、速度分布

zp = 500;                                                            %所画 z 平面
zpi = zp-z0+1;
X = permute(xlab(:,zpi,:),[1,3,2]);
Y = permute(ylab(:,zpi,:),[1,3,2]);

VR = permute(Vr_lab(:,zpi,:),[1,3,2]);
VT = permute(Vtheta_lab(:,zpi,:),[1,3,2]);


figure
VRi = VR;                                                           %VRi 矩阵设置图的颜色，使图颜色均衡
VRi(VR>0.001) = 0.001;
VRi(VR<-0.001) = -0.001;
ttt = surf([X;X(1,:)], [Y;Y(1,:)],[VR;VR(1,:)],[VRi;VRi(1,:)]);
ttt.EdgeColor = 'none';
view(0,90)
box
figure
VTi = VT;                                                           %VTi 矩阵设置图的颜色，使图颜色均衡
VTi(VT>0.003) =0.003;     
VTi(VT<0.0001) =0.0001;     
ttt = surf([X;X(1,:)], [Y;Y(1,:)],[VT;VT(1,:)],[VTi;VTi(1,:)]);
ttt.EdgeColor = 'none';
view(0,90)
box
