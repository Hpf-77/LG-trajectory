%p��ķ��ͼ

%��ĳһ z ƽ��Ĺ�ǿͼ
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

%��ʼ���ݴ�����ֱ������ϵת��Ϊֱ������ϵ
clear pi;
zr = (pi * w0^2) / lamd ;     %��������

zlab = z0:h:zmax; 
wz = w0 * sqrt(1 + (zlab/zr).^2);  
R = zlab + zr^2 ./ zlab;

 [theta1,rlab1] = cart2pol(real(xlab),real(ylab));                          %�����ռ伫���껯
[theta_p,r_p] = cart2pol(real(x_p_lab),real(y_p_lab));                  %λ�ÿռ伫���껯
 
 %p��ķ�����£��켣�������ٶ�
 Vxp = lamd /(2*pi)*xlab;
 Vyp = lamd /(2*pi)*ylab;
 %λ�ÿռ� r�����ٶ�
Vrp_lab = Vxp.*cos(theta_p) + Vyp.*sin(theta_p);
%λ�ÿռ� theta������ٶ�
 Vthetap_lab = -1/r_p .* (Vxp.*sin(theta_p) - Vyp.*cos(theta_p));
 collab = [1,0,0;0,0,0;0,0,1];


%  r�����ٶ�Vr��theta���� ���ٶ�w ��z �����ı仯ͼ
   omg = -1*Vthetap_lab;                                 
figure
 hold on
 ii = 1;
  for tempi2 =  [6,46,96]
    plot(zlab(1:4500)./zr, omg(24,1:4500,tempi2),'color',collab(ii,:));%200:1000:10000
    ii = ii+1;
  end
    xlabel('z/zr')
  ylabel('���ٶ�w')
  title('z-wͼ(max)')
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
  title('z-Vrͼ(max)');
  box



%  
%  
%  ����Ϊ����ķ�켣��ά�ռ�ͼ����
modelmax = 9504;                                       %ÿ��ƽ���ϵ�����ȡ��������                   �뾶�����˽Ƕ�������   
th =1;                                                           % �Ƕȼ��   6��ÿ��60��ȡһ���㡣    
zh = 5000;                                                     %ÿ��z=zhȡһ������ٶ�   %zƽ����
zhmaxi = zh-z0;                                             %(�������������ֵ)


picmodel(abs(picmodel)>100 ) = 10^(-15);         %��㸽��ֵ�� 0 
picmodel(picmodel<10^(-255) ) = 10^(-255);     %��㸽��ֵ�� 0 
 
%ÿ����Ĺ�ǿ��һ������
for i = 1:b-1
    lg_c(i) = sum(picmodel(:,i));
end

[Lg_c,~] = meshgrid (lg_c,picmodel(:,1));
Picmodel = picmodel ./ Lg_c;                        %��һ��
Picmodel = Picmodel * 3;
Picmodel(Picmodel<10^(-500) ) = 10^(-500);


figure 
xlabel('x/mm')
ylabel('z/mm')
zlabel('y/mm')
%pbohm�¶����ռ� ��άɢ��ͼ
%ĳ��zƽ���ϣ����ݲ�ķ������õ�����ģ�⶯���ռ��ǿ�ֲ�
zzlab =([500,5000] -z0 ); 
for phi = [1,2]
    pz = zzlab(phi)+1;                                                     % z ƽ�� ������
hold on
 scatter3(Xmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Xmodel(1:th:modelmax,pz))), Ymodel(1:th:modelmax,pz),5000*real(Picmodel(1:th:modelmax,pz)),1000*real(Picmodel(1:th:modelmax,pz)),'filled');
end
% pbohm�¶����ռ� ��ά�켣ͼ

collab = [1,0,0;0,0,0;0,0,1];
 ii = 1;
 
hold on
for t =[6,46,96]                                                  %  �뾶 ������

               for tp = [96,24,48,72]                         %  �Ƕ� ������
            hold on
                    plot3(xlab(tp,1:1:zhmaxi,t),zlab(1:1:zhmaxi),ylab(tp,1:1:zhmaxi,t),'color',collab(ii,:),'linewidth',0.5);
                    scatter3(xlab(tp,1 ,t),zlab( 1 ),ylab(tp,1,t), 50,collab(ii,:),'*' );
               end
ii = ii+1;
end

%pbohm��λ�ÿռ� ��άɢ��ͼ
figure
hold on
zzlab =([500,5000]-z0) ;
for phi = [1,2]
    pz = zzlab(phi)+1;                                                       % z ƽ�� ������
hold on
  scatter3(Xpmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Xpmodel(1:th:modelmax,pz))), Ypmodel(1:th:modelmax,pz),5000*real(Picmodel(1:th:modelmax,pz)),1000*real(Picmodel(1:th:modelmax,pz)),'filled');
end

%pbohm��λ�ÿռ� ��ά�켣ͼ
ii=1;

for t = [6,46,96]                                                   %  �뾶 ������

               for tp =[96,24,48,72]                            %  �Ƕ� ������

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


