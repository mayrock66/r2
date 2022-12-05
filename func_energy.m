function N_den=func_energy(E,tt,U_bias,A,spB_s,spB_d)

global eta Nx Vs Vd k_B Temp q fermi_flag;

ee=E;ep=ee+eta;
ck=1-((ep-U_bias(1))/(2*tt));con_s=-tt*exp(i*acos(ck));
ck=1-((ep-U_bias(Nx))/(2*tt));con_d=-tt*exp(i*acos(ck));
U_eff=U_bias;
U_eff(1)=U_bias(1)+con_s;
U_eff(Nx)=U_bias(Nx)+con_d;
G_inv=(ep*eye(Nx))-A-diag(U_eff);
G_s=G_inv\spB_s;
G_d=G_inv\spB_d;
f_1=fermi(((-Vs-ee)/(k_B*Temp/q)),fermi_flag,-1/2);
f_2=fermi(((-Vd-ee)/(k_B*Temp/q)),fermi_flag,-1/2);
N_den=-abs(G_s).^2*imag(con_s)*2*f_1...
           -abs(G_d).^2*imag(con_d)*2*f_2;
