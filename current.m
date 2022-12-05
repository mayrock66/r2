%********************************************************************************
%CALCULATE THE CURRENT DENSITY [A/m] GIVEN A 1D CHARGE SHEET POTENTIAL PROFILE***
%(Zhibin Ren 7-18-00)************************************************************
%********************************************************************************

function [Ie,Ie_sub,Te_sub,Fn_sub]=current(Ne,Ec,Ne_sub,E_sub)


%FUNDAMENTAL physical constants

global eps_si eps_o psi_si q k_B h_bar m_e m_t m_l Nc Ncc mx my mz Temp;
global transport_model Vs Vd
global dummy_flag fermi_flag ox_pnt_flag t_vall max_subband;
global delta_x Lsd Lg_top Lg_bot t_top t_bot t_si dx dy;
global phi_top phi_bot eps_top eps_bot bar_top bar_bot;
global eps_si m_t m_l;
global Nx Ny Ntotal;
global N_sd N_body;    
global mu_low beta Vel_sat criterion_outer
global ELE_TAUW ELE_CQ delta_T_1 delta_T_2 dim_c;
global nu_scatter Info_scatter_old;
global Trans N_dos E;

%%%%MODEL 2 related parameters

eta=1e-6*i;
Ef_tail_low=0.01;
Ef_tail_up=0.3;
E_step=criterion_outer/2;%0.5 times criterion_outer

%%%%MODEL 4 related parameters
e_field_limit=1.0;%in V/m
%%%%MODEL 5 related parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	INPUT AND OUTPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Info_scatter_old contains all information of scatters
%2-dim matrix, dim_2:current,index, and Fermi energy, low field mobility
%Ne_old is the old electron density
%1 column of Ntotal elements
%Ec_old is the old potential energy profile in eV
%1 column of Ntotal elements
%Te_old is the old electron temperature along the channel
%1 column of Nx elements, for charge-sheet-model
%Ne_sub contains 2D electron density on each subband
%E_sub contains the 1D subband energy for all subband 
%considered
%Te_sub_old contains the 1D electron temperature for all subband
%considered, for bulk model

%Ne is the 3D electron density profile in the entire device
%1 column matrix of Ntotal elements
%Ec is the conduction band edge profile in the entire device
%1 column matrix of Ntotal elements
%Te_old is 1 column matrix of Nx elements for electron temperature
%, for charge-sheet-model

%Ie is 1 column matrix of Nx elements for current density
%A/m
%Ie_sub contains currents carried in all subbands.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	TEMPORARY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ec_channel is the given potential energy profile 
%along the channel
%Fn_channel is the given quasi Fermi level along the channel
%Ne_channel is the given electron 2D density along the channel
%max_subband is the number of subbands considered
%t_vall is the number of valleyes considered

%************************INITIALIZATION***************************************
Ie=zeros(Nx,1);
Ie_sub=zeros(Nx,max_subband,t_vall);
Te_sub=zeros(Nx,max_subband,t_vall);
Fn_sub=zeros(Nx,max_subband,t_vall);
%E_sub=zeros(Nx,max_subband,t_vall);
%Ne_sub=zeros(Nx,max_subband,t_vall);

%Info_scatter_new=zeros(nu_scatter,4);
%Info_scatter_new=Info_scatter_old;
Ec_channel=zeros(Nx,1);
Fn_channel=zeros(Nx,1);
Ne_channel=zeros(Nx,1);
%Fn_new=zeros(Ntotal,1);
%Ne_new=zeros(Ntotal,1);

%temporarily used variables
if ox_pnt_flag==0
   Np_v=round(t_si/dy)+1;
elseif ox_pnt_flag==1
   Np_v=Ny;
end

Ns=zeros(Nx,1);
N_body=zeros(Np_v,Nx);
N_body_sum=zeros(Np_v,Nx);
U_sub=zeros(Nx,max_subband,t_vall);
W_sub=zeros(Np_v,Nx,max_subband,t_vall);
U_bias=zeros(Nx,1);

%*****************************************************************************
%THE START OF TRANSPORT EQUATION PART*****************************************
%*****************************************************************************

%*****************************************************************************
 if transport_model==1%TRANSPORT MODEL 1**************************************
%*****************************************************************************

%*****************************************************************************
 elseif transport_model==2%TRANSPORT MODEL 2 DD*******************************
%*****************************************************************************

U_bias=zeros(Nx,1);
Ns=zeros(Nx,1);
Ie_tem=zeros(Nx,1);
mobility=zeros(nu_scatter-1,1);
mobility=(diff(Info_scatter_old(:,4))+2*Info_scatter_old(1:Nx-1,4))/2;

for i_val=1:t_vall
   Nccc=2*m_e*sqrt(mx(i_val)*my(i_val))*k_B*Temp/(pi*h_bar^2);
   
   for i_sub=1:max_subband
      U_bias=E_sub(:,i_sub,i_val);
      Ns=Ne_sub(:,i_sub,i_val);
      
      if fermi_flag==0
        fac_deg=ones(Nx-1,1);
      elseif fermi_flag==1
        Fn_channel=U_bias+(k_B*Temp/q)*anti_dummy(Ns./Nccc,0,fermi_flag);
        zeta_channel=(Fn_channel-U_bias)./(k_B*Temp/q);
        fac_deg_tem=log(1+exp(zeta_channel)).*(1+exp(-zeta_channel));
        fac_deg=(diff(fac_deg_tem)+2*fac_deg_tem(1:Nx-1))/2;
      end
      V_channel=-U_bias;
      E_channel=-diff(V_channel)/dx;
      mu_channel=mobility./(1+(mobility./Vel_sat.*abs(E_channel)).^beta).^(1/beta);

      for j=1:Nx-1
        if abs(E_channel(j))<=abs(e_field_limit*fac_deg(j))
           Coe_1(j)=-mu_channel(j)*(k_B*Temp/q)...
                             *fac_deg(j)/dx;
           Coe_2(j)=-Coe_1(j);
        else
           Coe_1(j)=mu_channel(j)*E_channel(j)...
                             *1/(1-exp(E_channel(j)*dx/...
                             ((k_B*Temp/q)*fac_deg(j))));
           Coe_2(j)=mu_channel(j)*E_channel(j)...
                             *1/(1-exp(-E_channel(j)*dx/...
                             ((k_B*Temp/q)*fac_deg(j))));
        end
        Ie_tem(j)=-q*(Ns(j)*Coe_1(j)+Ns(j+1)*Coe_2(j));
      end
      Ie_tem(Nx)=Ie_tem(Nx-1);
      Ie_sub(:,i_sub,i_val)=Ie_tem;
      Fn_sub(:,i_sub,i_val)=Fn_channel;
      Ie=Ie+Ie_tem;
   end
end

%************************************************************************** 
elseif transport_model==3%BALLISTIC TRANSPORT USING SEMICLASSICAL APPROACH
%**************************************************************************
U_bias=zeros(Nx,1);

for i_val=1:t_vall
   Ie_2d=2*q/h_bar^2*sqrt(mx(i_val)*m_e/2)*((k_B*Temp/q)*q/pi)^(3/2);

   for i_sub=1:max_subband
      U_bias=E_sub(:,i_sub,i_val);
      Ec_peak=max(U_bias);
      
      Ie_tem=0;
      Ie_tem=Ie_2d*(fermi(((-Vs-Ec_peak)/(k_B*Temp/q)),fermi_flag,1/2)...
         -fermi(((-Vd-Ec_peak)/(k_B*Temp/q)),fermi_flag,1/2));
      Ie_sub(:,i_sub,i_val)=Ie_tem*ones(Nx,1);
      Ie=Ie+Ie_tem*ones(Nx,1);       
   end
end

%*********************************************************************************
elseif transport_model==4%BALLISTIC TRANSPORT MODEL USING GREEN FUNCTION APPROACH*
%*********************************************************************************

U_bias=zeros(Nx,1);

Ec_peak=max(-Vs,max(E_sub(:,max_subband,t_vall)));
E_number=round((Ec_peak+Ef_tail_up-E_sub(Nx,1,1)+Ef_tail_low)/E_step)+2;
E=linspace((E_sub(Nx,1,1)-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number);
delta_E=((Ec_peak+Ef_tail_up)-(E_sub(Nx,1,1)-Ef_tail_low))/(E_number-1);

N_dos_one=zeros(Nx,E_number);
N_dos=zeros(Nx,E_number);

%used to plot the transmission coefficient with several valleys
Trans_sub=zeros(E_number,1);
Trans=zeros(E_number,1);

for i_val=1:t_vall
   Ne_2d=2*sqrt(mx(i_val)*m_e*(k_B*Temp/q)*q/(2*pi^3))/(h_bar*dx);
   Ie_2d=2*q^2/(pi^2*h_bar^2)*sqrt(mx(i_val)*m_e*(k_B*Temp/q)*q*pi/2);
   tt=(h_bar^2)/(2*my(i_val)*m_e*(dx^2)*q);
   A=tt*((2*eye(Nx))-(diag(ones(Nx-1,1),1))-(diag(ones(Nx-1,1),-1)));

   for i_sub=1:max_subband
     U_bias=E_sub(:,i_sub,i_val);
     %Ec_peak=max(-Vs,max(U_bias));
     %E_number=round((Ec_peak+Ef_tail_up-U_bias(Nx)+Ef_tail_low)/E_step)+2;
     %E=linspace((U_bias(Nx)-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number);
     %delta_E=((Ec_peak+Ef_tail_up)-(U_bias(Nx)-Ef_tail_low))/(E_number-1);
     
     B_d=zeros(Nx,1);
     B_d(Nx)=1;
     spB_d=sparse(B_d);
     B_s=zeros(Nx,1);
     B_s(1)=1;
     spB_s=sparse(B_s);
     

     Ie_tem=0;
     for k=1:E_number,
       ee=E(k);ep=ee+eta;
       ck=1-((ep-U_bias(1))/(2*tt));con_s=-tt*exp(i*acos(ck));
       ck=1-((ep-U_bias(Nx))/(2*tt));con_d=-tt*exp(i*acos(ck));
       U_eff=U_bias;
       U_eff(1)=U_bias(1)+con_s;
       U_eff(Nx)=U_bias(Nx)+con_d;
       G_inv=sparse((ep*eye(Nx))-A-diag(U_eff));
       G_s=G_inv\spB_s;
       G_d=G_inv\spB_d;
       f_1=fermi(((-Vs-ee)/(k_B*Temp/q)),fermi_flag,-1/2);
       f_2=fermi(((-Vd-ee)/(k_B*Temp/q)),fermi_flag,-1/2);
       Ie_tem=Ie_tem+abs(G_d(1))^2*imag(con_s)*imag(con_d)*4*(f_1-f_2);
       Trans(k)=abs(G_d(1))^2*imag(con_s)*imag(con_d)*4;
       N_dos_one(:,k)=-abs(G_s).^2*imag(con_s)-abs(G_d).^2*imag(con_d);
       
     end
     
     N_dos=N_dos+Ne_2d*N_dos_one;
     
     Ie_sub(:,i_sub,i_val)=Ie_2d*delta_E*Ie_tem*ones(Nx,1);
     Ie=Ie+Ie_2d*delta_E*Ie_tem*ones(Nx,1);
     Trans_sub=Trans_sub+Trans;
     
   end
end

Trans=Trans_sub;

%*********************************************************************************
elseif transport_model==5%SCATTERING APPROACH USING NEGF METHOD*******************
%*********************************************************************************

%*********************************************************************************
elseif transport_model==6%ENERGY TRANSPORT MODEL *********************************
%*********************************************************************************
U_bias=zeros(Nx,1);
Ns=zeros(Nx,1);
Ie_tem=zeros(Nx,1);

for i_val=1:t_vall
   Nccc=2*m_e*sqrt(mx(i_val)*my(i_val))*k_B*Temp/(pi*h_bar^2);
   
   for i_sub=1:max_subband
      U_bias=E_sub(:,i_sub,i_val);
      Ns=Ne_sub(:,i_sub,i_val);
      
    %===============================
    %  Nondegenerate always for ET
    %===============================
      Ns_con=Nccc*exp((-Vs-U_bias(1))/(k_B*Temp/q));
      Nd_con=Nccc*exp((-Vd-U_bias(Nx))/(k_B*Temp/q));

      V_channel=-U_bias;

      [Ne_tem,Ie_tem,Te_tem]=ET(V_channel,mu_low,Ns_con,Nd_con,...
				Temp,Temp,ELE_TAUW,ELE_CQ,Vel_sat,beta);

      Ie_tem(Nx)=Ie_tem(Nx-1);
      Ie_sub(:,i_sub,i_val)=Ie_tem;
      Ie=Ie+Ie_tem;
      Te_sub(:,i_sub,i_val)=Te_tem;
   end
end


end

%***************************************************************************
%****************************END OF CURRENT SOLVER**************************
%***************************************************************************
