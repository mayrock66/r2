%Dummy function used to convert quasi Fermi energy to Ne
%(Zhibin Ren 6-5-00)

function [y]=dummy_prime(x,dummy_flag,fermi_flag)

if dummy_flag==0
  if fermi_flag==0
    y=exp(x);
  elseif fermi_flag==1
    y=1./(1+exp(-x));
  end
elseif dummy_flag==1/2
  y=fermi(x,fermi_flag,-1/2);
end

