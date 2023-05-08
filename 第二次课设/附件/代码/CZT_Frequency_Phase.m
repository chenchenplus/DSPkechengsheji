function [f_estimate,delta_f_estimate,phase_estimate,X_CZT] = CZT_Frequency_Phase(xn,fs,N,iter_num)
X_CZT=fft(xn);
delta_f=fs;
[~,index]=max(abs(X_CZT));
f_start=0;
for i=1:iter_num
    delta_f=delta_f/N;
    f_start=f_start+(index-1.5)*delta_f;
    W=exp(-1j*2*pi*(delta_f/N)/fs);
    A=exp(1j*2*pi*f_start/fs);
    X_CZT=czt(xn,N,W,A);
    [~,index]=max(abs(X_CZT));
    f_estimate=f_start+(index-1)*delta_f/N;
end
delta_f_estimate=delta_f/N;
phase_estimate=angle(X_CZT(index));