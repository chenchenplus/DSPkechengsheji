function [r_estimate,delta_r_estimate] = get_r(f_estimate,delta_f_estimate,k,T,N,fs,f0,c0)
r_estimate = mod(f_estimate/(2*k),T)*c0/2;
delta_r_estimate=c0*pi*delta_f_estimate*(N-1)/(4*f0*fs);
end