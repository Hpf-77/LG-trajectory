%xbohm�����»�ͼ�ļ���



%��ĳһ z ƽ��Ĺ�ǿͼ
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

%��ʼ���ݴ�����ֱ������ϵת��Ϊֱ������ϵ
clear pi;
zr = (pi * w0^2) / lamd ;     %��������

zlab = z0:h:zmax; 
wz = w0 * sqrt(1 + (zlab/zr).^2); 
R = zlab + zr^2 ./ zlab;

 [Vtheta1,Vrlab1] = cart2pol(Vxlab,Vylab);
 [theta1,rlab1] = cart2pol(real(xlab),real(ylab));

%r�����ٶ�
Vr_lab = Vxlab.*cos(theta1) + Vylab.*sin(theta1);
%theta����  ���ٶ�
 Vtheta_lab = -1./rlab1 .* (Vxlab.*sin(theta1) - Vylab.*cos(theta1));

 %r ������
p_theta = (2*pi)/lamd * Vtheta_lab;
% theta ������
p_r = (2*pi)/lamd * Vr_lab;

collab = [1,0,0;0,0,0;0,0,1];
 
%  �뾶r���Ƕ�theta��r�����ٶ�Vr��theta���� ���ٶ�w ��z �����ı仯ͼ
% figure
%  hold on
%  ii=1;
%  for tempi =  [30,50,70]
%     plot(zlab(1:10)./zr, theta1(96,1:10,tempi),'color',collab(ii,:)); %'-*','markerindices',200:1000:10000,
%   %legend_str{1}=['l=3,r_0=0.2',num2str(1)];
%   ii=ii+1;
%  end
% xlabel('z/zr')
%   ylabel('�Ƕ�theta')
%   title('z-thetaͼ')
% 
%  figure
%  %for fi = 1:16
%      %figure(fi)
%  hold on
%   for tempi2 = 1:6%[3,17,30]
%     plot(zlab(1:10000), rlab1(1,1:10000,tempi2));  %'color',[0.5,0.5,0.5],'-*','markerindices',800:1000:10000,'color','g'
%   end
%   xlabel('z/zr')
%   ylabel('�뾶r')
%   title('z-rͼ')

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
  ylabel('���ٶ�w/c')
  title('z-wͼ')
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
  title('z-Vrͼ');
  box
  


%  ����Ϊ����ķ�켣��ά�ռ�ͼ����

modelmax = 9504;                                         %ÿ��ƽ���ϵ�����ȡ��������                   �뾶�����˽Ƕ�������   
th =1;                                                              % �Ƕȼ��   6��ÿ��60��ȡһ���㡣    
zh = 2000;                                                       %��Զ z
zhmaxi = zh-z0;                                              %(�������������ֵ)

 
picmodel(picmodel>100 ) = 10*(-15);    %��㸽��ֵ�� 0 

%ÿ����Ĺ�ǿ��һ������
for i = 1:b-1
    lg_c(i) = sum(picmodel(:,i));
end

[Lg_c,~] = meshgrid (lg_c,picmodel(:,1));
Picmodel = picmodel ./ Lg_c;                        %��һ��
Picmodel = Picmodel * 3;

%��άɢ��ͼ                                                     
%ĳ��zƽ���ϣ����ݲ�ķ������õ�����ģ��λ�ÿռ��ǿ�ֲ�
figure
zzlab =([500,2000]-z0) ;                                  %����zƽ��

for phi = 1:2                                                       % z ƽ�� ������
    pz = zzlab(phi)+1;
hold on
%scatter������ÿ��λ�õ��Ӧ��ǿ������ɫ��ԲȦ��С��ʾ����ģ���ǿ�ֲ���
 scatter3(Xmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Xmodel(1:th:modelmax,pz))), Ymodel(1:th:modelmax,pz),7000*real(Picmodel(1:th:modelmax,pz)),100*real(Picmodel(1:th:modelmax,pz)),'filled')
end

% ��ά�켣ͼ
%xbohm������λ�ÿռ�켣
hold on
collab = [1,0,0;0,0,0;0,0,1];
 ii = 1;
for t = [4,49,92]                                                     %  �뾶 ������
               for tp =1                                              %  �Ƕ� ������
                    xyzi = tp;
                    hold on
                    plot3(xlab(tp,1:1:zhmaxi ,t),zlab( 1:1:zhmaxi ),ylab(tp,1:1:zhmaxi  ,t),'color',collab(ii,:),'linewidth',0.5);

               end
ii = ii+1;
end


%xbohm������ �������ռ�켣
figure 
hold on
grid on
Pxmodel = zeros(size (Xmodel));
Pymodel  =Pxmodel;
for i = 1:1:b-1                                                     %��ά�����켣ת��һά����
   for t = 1:rb    
        for tp = 1:1:jmax
             xyzi = (t-1)*jmax+tp;
            Pxmodel(xyzi,i) = pxlab(tp,i,t) ;          
            Pymodel(xyzi,i) =pylab(tp,i,t) ;
        end
   end
end
zzlab =([500,2000]-z0) ; 

%��άɢ��ͼ                                                     
%ģ��ĳ�� z ƽ���϶����ռ��ڵĹ�ǿ�ֲ���
for phi = 1:2
    pz = zzlab(phi)+1;                                          % z ƽ�� ������
hold on
 scatter3(Pxmodel(1: th:modelmax,pz), (pz+z0)*ones(size(Pxmodel(1:th:modelmax,pz))), Pymodel(1:th:modelmax,pz) ,7000*real(Picmodel(1:th:modelmax,pz)),100*real(Picmodel(1:th:modelmax,pz)),'filled')
end

%��ά�����ռ�켣ͼ
ii = 1;
for t = [4,49,92]                                                 %  �뾶 ������ 
               for tp =[96,24,48,72]                          %  �Ƕ� ������
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


%��ĳ z ƽ���Ͻ��ٶȡ��ٶȷֲ�

zp = 500;                                                            %���� z ƽ��
zpi = zp-z0+1;
X = permute(xlab(:,zpi,:),[1,3,2]);
Y = permute(ylab(:,zpi,:),[1,3,2]);

VR = permute(Vr_lab(:,zpi,:),[1,3,2]);
VT = permute(Vtheta_lab(:,zpi,:),[1,3,2]);


figure
VRi = VR;                                                           %VRi ��������ͼ����ɫ��ʹͼ��ɫ����
VRi(VR>0.001) = 0.001;
VRi(VR<-0.001) = -0.001;
ttt = surf([X;X(1,:)], [Y;Y(1,:)],[VR;VR(1,:)],[VRi;VRi(1,:)]);
ttt.EdgeColor = 'none';
view(0,90)
box
figure
VTi = VT;                                                           %VTi ��������ͼ����ɫ��ʹͼ��ɫ����
VTi(VT>0.003) =0.003;     
VTi(VT<0.0001) =0.0001;     
ttt = surf([X;X(1,:)], [Y;Y(1,:)],[VT;VT(1,:)],[VTi;VTi(1,:)]);
ttt.EdgeColor = 'none';
view(0,90)
box
