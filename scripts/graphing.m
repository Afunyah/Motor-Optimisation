

% Motion Profile Graph
% figure('DefaultAxesFontSize',16);
% hold on;
% title('Trapezoidal Motion Profile');
% plot(out.tout, out.mp_data.ang_acc.Data,'b-','LineWidth',1.5, 'Color', "#0072BD");
% plot(out.tout, out.mp_data.ang_vel.Data,'LineWidth',1.5, 'Color', "#EDB120");
% plot(out.tout, out.mp_data.ang_pos.Data,'LineWidth',1.5, 'Color', "#D95319");
% legend({'Acceleration','Velocity','Position'},'Location','northwest');
% xlabel('Time/s');
% hold off;
% grid on;
% 
% %------------------
% figure('DefaultAxesFontSize',16);
% hold on;
% title('Load Trajectory');
% plot(out.tout, out.load_data.signal_0__.Data,'LineWidth',1.5, 'Color', "#0072BD");
% plot(out.tout, out.load_data.signal_0_.Data,'LineWidth',1.5, 'Color', "#EDB120");
% plot(out.tout, out.load_data.signal_0.Data,'LineWidth',1.5, 'Color', "#D95319");
% legend({'Acceleration','Velocity','Position'},'Location','northwest');
% xlabel('Time/s');
% hold off;
% grid on;
% 
% 
% 
% %------------------
% 
% figure('DefaultAxesFontSize',16);
% subplot(3,2,3);
% plot(out.tout, out.mp_data.ang_vel.Data,'LineWidth',1.5, 'Color', "#EDB120");
% title('Expected Velocity');
% xlabel('Time/s');
% ylabel('Velocity/rads^{-1}');
% grid on;
% 
% subplot(3,2,4);
% plot(out.tout, out.load_data.signal_0_.Data,'LineWidth',1.5, 'Color', "#EDB120");
% title('Load Velocity');
% xlabel('Time/s');
% ylabel('Velocity/rads^{-1}');
% grid on;
% 
% subplot(3,2,1);
% plot(out.tout, out.mp_data.ang_pos.Data,'LineWidth',1.5, 'Color', "#D95319");
% title('Expected Position');
% xlabel('Time/s');
% ylabel('Position/rad');
% grid on;
% 
% subplot(3,2,2);
% plot(out.tout, out.load_data.signal_0.Data,'LineWidth',1.5, 'Color', "#D95319");
% title('Load Position');
% xlabel('Time/s');
% ylabel('Position/rad');
% grid on;
% 
% subplot(3,2,5);
% plot(out.tout, out.mp_data.ang_acc.Data,'LineWidth',1.5, 'Color', "#0072BD");
% title('Expected Acceleration');
% xlabel('Time/s');
% ylabel('Acceleration/rads^{-2}');
% grid on;
% 
% subplot(3,2,6);
% plot(out.tout, out.load_data.signal_0__.Data,'LineWidth',1.5, 'Color', "#0072BD");
% title('Load Acceleration');
% xlabel('Time/s');
% ylabel('Acceleration/rads^{-2}');
% grid on;
% 
% %------------------
% % Error Graphs
% figure('DefaultAxesFontSize',16);
% subplot(3,1,2);
% plot(out.tout, out.mp_data.ang_vel.Data-out.load_data.signal_0_.Data,'LineWidth',1.5, 'Color', "#EDB120");
% title('Velocity Error');
% xlabel('Time/s');
% ylabel('Velocity/rads^{-1}');
% grid on;
% 
% subplot(3,1,1);
% plot(out.tout, out.mp_data.ang_pos.Data-out.load_data.signal_0.Data,'LineWidth',1.5, 'Color', "#D95319");
% title('Position Error');
% xlabel('Time/s');
% ylabel('Position/rad');
% grid on;
% 
% subplot(3,1,3);
% plot(out.tout, out.mp_data.ang_acc.Data-out.load_data.signal_0__.Data,'LineWidth',1.5, 'Color', "#0072BD");
% title('Acceleration Error');
% xlabel('Time/s');
% ylabel('Acceleration/rads^{-2}');
% grid on;
% 
% 
% %-------------------
% figure('DefaultAxesFontSize',16);
% subplot(2,1,1);
% plot(out.tout, out.motor_data.i.Data,'LineWidth',1.5, 'Color', "#0000FF");
% title('Motor Current');
% xlabel('Time/s');
% ylabel('Current/A');
% grid on;
% 
% subplot(2,1,2);
% plot(out.tout, out.motor_data.Tm.Data,'LineWidth',1.5, 'Color', "#FF0000");
% title('Motor Torque');
% xlabel('Time/s');
% ylabel('Torque/Nm');
% grid on;


%------------------
















% figure('DefaultAxesFontSize',16);
% hold on;
% plot(kt1, k1,'LineWidth',1.5);
% plot(kt2, k2,'LineWidth',1.5);
% plot(kt3, k3,'LineWidth',1.5);
% title('Energy Used against Time');
% xlabel('Time/s');
% ylabel('Energy/J');
% legend({'k = 0.2','k = 0.55','k = 0.8'},'Location','northwest');
% grid on;
% hold off;



