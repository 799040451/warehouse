x=0:pi/50:2*pi;
y=sin(x);
z=cos(x);
plot(x,y,'r',x,z,'g');               %绘制两根不同曲线
Hl=get(gca,'Children');           %获取两曲线句柄向量Hl
for k=1:size(Hl)
   if get(Hl(k),'Color')==[0 1 0]    %[0 1 0]代表绿色
       Hlg=Hl(k);               %获取绿色线条句柄
   end
end
pause                              %便于观察设置前后的效果
set(Hlg, 'LineStyle',':', 'Marker','p');       %对绿色线条进行设置
