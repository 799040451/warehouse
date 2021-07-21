N=128;                         % ��������
T=1;                                 % ����ʱ���յ�
t=linspace(0,T,N);                % ����N������ʱ��ti(I=1:N)
x=12*sin(2*pi*10*t+pi/4)+5*cos(2*pi*40*t);  % �������������ֵx
dt=t(2)-t(1);                   % ��������
f=1/dt;                          % ����Ƶ��(Hz)
X=fft(x);                        % ����x�Ŀ��ٸ���Ҷ�任X
F=X(1:N/2+1);                   % F(k)=X(k)(k=1:N/2+1)
f=f*(0:N/2)/N;                  % ʹƵ����f���㿪ʼ
plot(f,abs(F),'-*')            % �������-Ƶ��ͼ
xlabel('Frequency');
ylabel('|F(k)|')