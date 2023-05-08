function [N,y] = get_FMCW(A0,f0,fs,sample_time,T,k,phi0,tau)
N=round(fs*sample_time);
i=0:N-1;
y=A0*cos(2*pi*(f0+k*mod(i/fs-tau,T)).*mod(i/fs-tau,T)+phi0);
end

