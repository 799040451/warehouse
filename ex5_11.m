x=0:0.1:10;
y=10*x.*x;
subplot(2,2,1);plot(x,y);
title('plot(x,y)');grid on;
subplot(2,2,2);semilogx(x,y);
title('semilogx(x,y)');grid on;
subplot(2,2,3);semilogy(x,y);
title('semilogy(x,y)');grid on;
subplot(2,2,4);loglog(x,y);
title('loglog(x,y)');grid on;
