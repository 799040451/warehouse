x=[0:0.5:360]*pi/180;
plot(x,sin(x),x,cos(x));
hold on
plot(x,x.^1.1)
plot(x,x.^1.3)