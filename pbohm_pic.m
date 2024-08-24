%p玻姆画图

%画某一 z 平面的光强图
% x = -2.1:0.01:2.1;
% [x,y] = meshgrid (x,x);
% z = 10000;
% ct = eval(LGlight_c);
% cct = conj(ct).*ct;
% cmax = max(max(max(cct)));
% cmin = min(min(min(cct)));
% figure
% hold on
% A1 = pcolor(x,y,cct);
% caxis([cmin,cmax]);
% set(A1,'edgecolor','none','facecolor','interp')

%初始数据处理，将直角坐标系转换为直角坐标系
clear pi;
zr = (pi * w0^2) / lamd ;     %瑞利距离

zlab = z0:h:zmax; 
wz = w0 * sqrt(1 + (zlab/zr).^2);  
R = zlab + zr^2 ./ zlab;

 [theta1,rlab1] = cart2pol(real(xlab),real(ylab));                          %动量空间极坐标化
[theta_p,r_p] = cart2pol(real(x_p_lab),real(y_p_lab));                  %位置空间极坐标化
 
 %p玻姆理论下，轨迹的运行速度
 Vxp = lamd /(2*pi)*xlab;
 Vyp = lamd /(2*pi)*ylab;
 %位置空间 r方向速度
Vrp_lab = Vxp.*cos(theta_p) + Vyp.*sin(theta_p);
%位置空间 theta方向角速度
 Vthetap_lab = -1/r_p .* (Vxp.*sin(theta_p) - Vyp.*cos(theta_p));
 collab = [1,0,0;0,0,0;0,0,1];


%  r方向速度Vr，theta方向 角速度w 随z 传播的变化图
   omg = -1*Vthetap_lab;                                 
figure
 hold on
 ii = 1;
  for tempi2 =  [6,46,96]
    plot(zlab(1:4500)./zr, omg(24,1:4500,tempi2),'color',collab(ii,:));%200:1000:10000
    ii = ii+1;
  end
    xlabel('z/zr')
  ylabel('角速度w')
  title('z-w图(max)')
box

   
  figure

 hold on
 ii = 1;
  for tempi2 =   [6,46,96]
    plot(zlab(1:4500)./zr, Vrp_lab(24,1:4500,tempi2),'color',collab(ii,:));   %'color',[0.5,0.5,0.5]
    ii = ii+1;
  end
  xlabel('z/zr')
  ylabel('Vr')       
  title('z-Vr图(max)');
  box



%  
%  
%  以下为画玻姆轨迹三维空间图代码
modelmax = 9504;                                       %每个平面上的数据取样点数量                   半径数量乘角度数据量   
th =1;                                                           % 角度间隔   6是每隔60°取一个点。    
zh = 5000;                                                     %每隔z=zh取一个点的速度   %z平面间隔
zhmaxi = zh-z0;                                             %(用于索引的最大值)


picmodel(abs(picmodel)>100 ) = 10^(-15);         %奇点附近值归 0 
picmodel(picmodel<10^(-255) ) = 10^(-255);     %奇点附近值归 0 
 
%每个点的光强归一化处理
for i = 1:b-1
    lg_c(i) = sum(picmodel(:,i));
end

[Lg_c,~] = meshgrid (lg_c,picmodel(:,1));
Picmodel = picmodel ./ Lg_c;                        %归一化
Picmodel = Picmodel * 3;
Picmodel(Picmodel<10^(-500) ) = 10^(-500);


figure 
xlabel('x/mm')
ylabel('z/mm')
zlabel('y/mm')
%pbohm下动量空间 三维散点图
%某个z平面上，根据玻姆理论求得的坐标模拟动量空间光强分布
zzlab =([500,5000] -z0 ); 
for phi = [1,2]
    pz = zzlab(phi)+1;                                                     % z 平面 索引号
hold on
 scatter3(Xmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Xmodel(1:th:modelmax,pz))), Ymodel(1:th:modelmax,pz),5000*real(Picmodel(1:th:modelmax,pz)),1000*real(Picmodel(1:th:modelmax,pz)),'filled');
end
% pbohm下动量空间 三维轨迹图

collab = [1,0,0;0,0,0;0,0,1];
 ii = 1;
 
hold on
for t =[6,46,96]                                                  %  半径 索引号

               for tp = [96,24,48,72]                         %  角度 索引号
            hold on
                    plot3(xlab(tp,1:1:zhmaxi,t),zlab(1:1:zhmaxi),ylab(tp,1:1:zhmaxi,t),'color',collab(ii,:),'linewidth',0.5);
                    scatter3(xlab(tp,1 ,t),zlab( 1 ),ylab(tp,1,t), 50,collab(ii,:),'*' );
               end
ii = ii+1;
end

%pbohm下位置空间 三维散点图
figure
hold on
zzlab =([500,5000]-z0) ;
for phi = [1,2]
    pz = zzlab(phi)+1;                                                       % z 平面 索引号
hold on
  scatter3(Xpmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Xpmodel(1:th:modelmax,pz))), Ypmodel(1:th:modelmax,pz),5000*real(Picmodel(1:th:modelmax,pz)),1000*real(Picmodel(1:th:modelmax,pz)),'filled');
end

%pbohm下位置空间 三维轨迹图
ii=1;

for t = [6,46,96]                                                   %  半径 索引号

               for tp =[96,24,48,72]                            %  角度 索引号

            hold on
                    plot3(x_p_lab(tp,1:1:zhmaxi,t),zlab(1:1:zhmaxi),y_p_lab(tp,1:1:zhmaxi,t),'color',collab(ii,:),'linewidth',0.15);%5*Wlineu
                    scatter3(x_p_lab(tp,1 ,t),zlab( 1 ),y_p_lab(tp,1,t),50,collab(ii,:), '*');
                    scatter3(x_p_lab(tp,zhmaxi ,t),zlab( zhmaxi ),y_p_lab(tp,zhmaxi,t), 50,collab(ii,:),'o');
               end
ii = ii+1;
end
xlabel('x/mm')
ylabel('z/mm')
zlabel('y/mm')
%axis([-10.5 10.5 500 5000 -10.5 10.5 ])
view(45,45);
grid on
box


