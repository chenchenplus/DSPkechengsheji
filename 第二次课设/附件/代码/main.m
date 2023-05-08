A0=1;
f0=4e5;
fs=1e6;
sample_time=1.5;
T=0.1;
k=1e3;
phi0=0;
tau=0.005;
iter_num=5;
c0=3e8;
sigma=[sqrt(100),sqrt(10),sqrt(1),sqrt(0.1),sqrt(0.01)];
[N,x1]=get_FMCW(A0,f0,fs,sample_time,T,k,phi0,0);
[~,x2]=get_FMCW(A0,f0,fs,sample_time,T,k,phi0,tau);
signal=x1.*x2;
for i=1:5
    noise=sigma(i)*randn(1,N);
    [f_estimate_sigma(i),delta_f_estimate,phase_estimate_sigma(i),X_CZT] = CZT_Frequency_Phase(signal+noise,fs,N,iter_num);
    [r_estimate_sigma(i),delta_r_estimate_sigma(i)]=get_r(f_estimate_sigma(i),delta_f_estimate,k,T,N,fs,f0,c0);
end
iter_num=[1,3,5,10,20];
for i=1:5
    noise=0;
    [f_estimate_iter(i),delta_f_estimate,phase_estimate_iter(i),X_CZT] = CZT_Frequency_Phase(signal+noise,fs,N,iter_num(i));
    [r_estimate_iter(i),delta_r_estimate_iter(i)]=get_r(f_estimate_iter(i),delta_f_estimate,k,T,N,fs,f0,c0);
end
fs=[1e5,3e5,5e5,8e5,1e6];
for i=1:5
    [N,x1]=get_FMCW(A0,f0,fs(i),sample_time,T,k,phi0,0);
    [~,x2]=get_FMCW(A0,f0,fs(i),sample_time,T,k,phi0,tau);
    signal=x1.*x2;
    noise=0;
    [f_estimate_fs(i),delta_f_estimate,phase_estimate_fs(i),X_CZT] = CZT_Frequency_Phase(signal+noise,fs(i),N,5);
    [r_estimate_fs(i),delta_r_estimate_fs(i)]=get_r(f_estimate_fs(i),delta_f_estimate,k,T,N,fs(i),f0,c0);
end
figure(1)
subplot(2,1,1)
plot(1:5,r_estimate_sigma);
grid on;
xlabel('试验次数'); ylabel('估算的距离R');
legend('估算的距离R');
subplot(2,1,2)
plot(1:5,log(delta_r_estimate_sigma));
xlabel('试验次数'); ylabel('估算的误差Re'); title('方差对CZT算法估计的影响');
legend('误差Re(dB)');
figure(2)
subplot(2,1,1)
plot(1:5,r_estimate_iter);
grid on;
xlabel('试验次数'); ylabel('估算的距离R');
legend('估算的距离R');
subplot(2,1,2)
plot(1:5,log(delta_r_estimate_iter));
xlabel('试验次数'); ylabel('估算的误差Re'); title('迭代次数对CZT算法估计的影响');
legend('误差Re(dB)');
figure(3)
subplot(2,1,1)
plot(1:5,r_estimate_fs);
grid on;
xlabel('试验次数'); ylabel('估算的距离R');
legend('估算的距离R');
subplot(2,1,2)
plot(1:5,log(delta_r_estimate_fs));
xlabel('试验次数'); ylabel('估算的误差Re'); title('采样频率fs对CZT算法估计的影响');
legend('误差Re(dB)');