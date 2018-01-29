function [X,Y]=simuin(Mx,Px,Nnoise,noiseX,noiseY,A)
%+++ SIMUIN simulation based on UVE paper.

if nargin<6;A=5;end
if nargin<5;noiseY=0.005;end
if nargin<4;noiseX=0.005;end
if nargin<3;Nnoise=0;end
if nargin<2;Px=100;end
if nargin<1;Mx=25;end
    
method='center';
x=rand(Mx,Px);
A=min([Mx Px A]);
B=repmat([1;-1],100,1);
x=pretreat(x,'center');

[u,s,v]=svd(x);
d=diag(s);
d=d(1:A)/sum(d);
T=x*v(:,1:A);
X0=T*v(:,1:A)';
noise=randn(Mx,Px+Nnoise);
noise=noiseX*noise/max(max(abs(noise)));
X=[X0 randn(Mx,Nnoise)*0.1]+noiseX;  % SIMUIN
Y=u(:,1:A)*B(1:A)+randn(Mx,1)*noiseY;
