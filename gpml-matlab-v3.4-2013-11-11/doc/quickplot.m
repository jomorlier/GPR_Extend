plot(xprev(:,1),mF_vel,'r');
hold on;
plot(xprev(:,2),mF_vel,'b');
plot(xprev(:,1),mF_accel,'g');
plot(xprev(:,2),mF_accel,'m');
xlabel('input');
ylabel('output');
legend('pos vs vel','vel vs vel','pos vs accel','vel vs accel;');