%多全息
clear all;
close all;
clc
lambda = 0.6328;
k0 = 2*pi./lambda;
pixel = 5.4;
M = 500;
N = 500;
lmax1 = M*pixel;
lmax2 = N*pixel;

x = linspace(-lmax2/2,lmax2/2,N);
y = linspace(-lmax1/2,lmax1/2,M);
[X,Y] = meshgrid(x,y);
g = double(imread('n3.tif'));
g = g(:,:,1);
g = imresize(g/max(g(:)),[M,N]);
Ex = g;
k1max = pi*N/lmax2;
k2max = pi*M/lmax1;
z = 20000;
z1 =26000;

g1 = double(imread('31.jpg'));
g1 = g1(:,:,1);
g1 = imresize(g1/max(g1(:)),[M,N]);
g2 = double(imread('33.jpg'));
g2 = g2(:,:,1);
g2 = imresize(g2/max(g2(:)),[M,N]);


mm=linspace(-k1max-pi/lmax1,k1max-pi/lmax1,N)/(2*pi);
nn=linspace(-k2max-pi/lmax2,k2max-pi/lmax2,M)/(2*pi);
[fx,fy] = meshgrid(mm,nn);
kz = sqrt(k0^2-(2*pi*fx).^2-(2*pi*fy).^2);

H = exp(1i*kz.*z);%正向传递函数
HH = exp(-1i*kz.*z);%逆向传递函数
H1 = exp(1i*kz.*z1);%正向传递函数
HH1 = exp(-1i*kz.*z1);
%%%...........让会产生倏逝波的逆向传递函数部分为零.............%%%
r = sqrt(fx.^2+fy.^2);
HH(find(r >= k0/2/pi)) = 0;
HH1(find(r >= k0/2/pi)) = 0;
for m = 1:300%迭代算法优化相位、振幅关系
    Ax = exp(1i*angle(ifft2(ifftshift(H.*(fftshift(fft2(Ex))))))).*sqrt(g1);
    phase_ij = (angle(ifft2(ifftshift((fftshift(fft2(Ax))).*HH))))/pi*180;
    Ex = sqrt(g).*exp(1i*phase_ij/180*pi);
    Ax = exp(1i*angle(ifft2(ifftshift(H1.*(fftshift(fft2(Ex))))))).*sqrt(g2);
    phase_ij = (angle(ifft2(ifftshift((fftshift(fft2(Ax))).*HH1))))/pi*180;
    Ex = sqrt(g).*exp(1i*phase_ij/180*pi);
end

figure;
imagesc(x,y,angle(Ex))%显示需要相位
phasef = angle(Ex)*180/pi;
title('相位分布')
figure;
imagesc(x,y,abs(Ex).^2);
colormap(gray)%显示振幅
title('振幅分布')%%%衍射计算验证效果（菲涅尔算法）%%%

Discrete_step = 2;% phase
Ex_Discrete = abs(Ex).*exp(1i*(ceil((angle(Ex)+pi+eps)/2/pi*Discrete_step)/Discrete_step*2*pi))/max(abs(Ex(:)));
Axx = ifft2(ifftshift(H.*(fftshift(fft2(Ex_Discrete)))));
Axx1 = ifft2(ifftshift(H1.*(fftshift(fft2(Ex_Discrete)))));
figure;imagesc(x,y,abs(Axx).^2)%显示浮现图像;
title('最后衍射浮现图像')
colormap(hot)
figure;imagesc(x,y,abs(Axx1).^2)%显示浮现图像;
title('最后衍射浮现图像')
colormap(hot)
b0=real(Ex_Discrete);
b=b0;
b1=[];
format long
for i=1:M
    b1=[b1 [b(i,:)]];
    b1=sort(b1);
end
j=[];
b1=b1.^2.*sign(b1);
e=12;
for i=1:e
    bmin=b1((i-1)*floor(M*N/e)+1);
    bmax=b1(i*floor(M*N/e));
    b11=linspace(bmin,bmax,100);
    b11=sqrt(abs(b11)).*sign(b11);
    b11=90*asin(b11)/pi;
    b11=round(b11);
    b11=sin(b11*pi/90);
    b11=b11.^2.*sign(b11);
    a=[];
    for k=1:length(b11)
        dd=sum(abs(b1((i-1)*floor(M*N/e)+1:i*floor(M*N/e))-b11(k)));
        a=[a dd];
    end
    h=min(a);
    f=find(a==h);
    j=[j b11(f(1))];
end
j=sqrt(abs(j)).*sign(j);
s=[j(1)];
for i=1:e
    m=0;
    for p=1:length(s)
        if j(i)==s(p)
            m=m+1;
        end
    end
    if m==0
        s=[s j(i)];
    end
end
j=s;
e=length(j);
j1=90*asin(j)/pi;
j1=round(j1);
j2=sin(j1*pi/90);
b2=[];
for m=1:M
    for n=1:N
        if b0(m,n)<=0
            j3=j2(j2<0);
            dd=abs(b0(m,n)-j3);
            h=min(dd);
            f=find(dd==h);
            b2(m,n)=j3(f(1));
        elseif b0(m,n)>0
            j3=j2(j2>0);
            dd=abs(b0(m,n)-j3);
            h=min(dd);
            f=find(dd==h);
            b2(m,n)=j3(f(1));
        end
    end
end
j4=abs(j1);
dis=j4-min(j4);
dis1=dis*30/max(dis);
j4=15+dis1;
j4=round(j4.*sign(j1));
j5=sin(j4*pi/90);
for m=1:M
    for n=1:N
        f=find(j2==b2(m,n));
        b2(m,n)=j5(f);
    end
end
b2(isnan(b2))=0;
figure;imagesc(x,y,abs(b2).^2);
colormap(gray)%显示振幅
title('振幅分布')
o=sum(sum(abs(b2).^2));
o0=sum(sum(abs(real(Ex_Discrete).^2)));
delta=abs((o-o0)/o0);
phase=zeros(M,N,e);
for m=1:M
    for n=1:N
        for p=1:e
            if b2(m,n)==j5(p)
                phase(m,n,p)=1;
            end
        end
    end
end
% %%
L = 600;
Ex_Discrete_extend = padarray((b2),[L,L],'both');
% Ex_Discrete_extend_phase = padarray(angle(Ex_Discrete),[L,L],'both');
figure;imagesc(Ex_Discrete_extend);title('object plane')
% Ex_Discrete_extend = Ex_Discrete_extend_amp.*exp(1i*Ex_Discrete_extend_phase);

M_new = length(Ex_Discrete_extend);
PB_modulation = Ex_Discrete_extend;
pixel = 5.4;
xmax=M_new*pixel;
ymax=M_new*pixel;
kmax = pi/pixel;
mm_kx = linspace(-kmax-pi/xmax,kmax-pi/xmax,M_new)/(2*pi);%
nn_ky = linspace(-kmax-pi/ymax,kmax-pi/ymax,M_new)/(2*pi);%
[MM_kx,NN_ky] = meshgrid(nn_ky,mm_kx);
% SAS diffraction
lambda = 0.6328;% um
% zz = z;
k0 = 2*pi/lambda;
kz = sqrt(k0^2-(2*pi*MM_kx).^2-(2*pi*NN_ky).^2);
H = exp(1i*kz.*z);%正向传递函数
H1 = exp(1i*kz.*z1);
Ex_imag1 = ifft2(ifftshift(H.*(fftshift(fft2(PB_modulation)))));
I_imag1 = abs(Ex_imag1).^2/max(abs(Ex_imag1(:)).^2);
figure;imagesc(I_imag1)
colormap(hot);
Ex_imag2 = ifft2(ifftshift(H1.*(fftshift(fft2(PB_modulation)))));
I_imag2 = abs(Ex_imag2).^2/max(abs(Ex_imag2(:)).^2);
figure;imagesc(I_imag2)
colormap(hot);
% I_imag1 = abs(Ex_imag1).^2;

% Ex_imag2 = ifft2(ifftshift(H.*(fftshift(fft2(exp(-1i*phase_matrix))))));
% I_imag2 = abs(Ex_imag2).^2;