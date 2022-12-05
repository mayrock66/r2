%This integral is used by charge.m sub routine for calculating 
%charge distribution in the transport model 3---semiclassical 
%version ballistic transport model

%Note:
%In order to improve the over simulation accuracy, you have to 
%decrease dx to proper values.

function [y]=integral(zeta_s,zeta_d,upper_lim,fermi_flag)
E_tail=12;
tem_lim=upper_lim^0.5;
if tem_lim>0
  dx=1e-2;
  %dx=1e-2 means dE=2xdx, so dE varies from 0eV to ~2meV@E=20KT/q 
  %the error criterion then can be around 1e-4 eV.
  nx=round(tem_lim/dx)+2;
  x=linspace(0,tem_lim,nx);
  dx=x(2);
  dummy1=fermi(zeta_s-x.^2,fermi_flag,-1/2);
  y1=(sum(dummy1)-dummy1(1)/2-dummy1(nx)/2)*dx;
else
  y1=0;
end

tem_lim_l=upper_lim^0.5;
tem_lim_h=(max(upper_lim+E_tail,zeta_d+E_tail))^0.5;
dx=1e-2;
%dx=1e-3 means dE=2xdx, so dE varies from ~0.2meV@E=20KT/q to ~0.3meV@E=40KT/q
%the error criterion then can be around 1e-4 eV.
nx=round((tem_lim_h-tem_lim_l)/dx)+2;
x=linspace(tem_lim_l,tem_lim_h,nx);
dx=x(2)-x(1);
dummy2=fermi(zeta_d-x.^2,fermi_flag,-1/2);
y2=(sum(dummy2)-dummy2(1)/2-dummy2(nx)/2)*dx;

tem_lim_l=0;
tem_lim_h=(max(upper_lim+E_tail,zeta_s+E_tail))^0.5;
dx=1e-2;
nx=round((tem_lim_h-tem_lim_l)/dx)+2;
x=linspace(tem_lim_l,tem_lim_h,nx);
dx=x(2)-x(1);
dummy3=fermi(zeta_s-x.^2,fermi_flag,-1/2);
y3=(sum(dummy3)-dummy3(1)/2-dummy3(nx)/2)*dx;

y=y1+y2+y3;


