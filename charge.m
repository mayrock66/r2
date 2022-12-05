%*****************************************************************************
%**********1D transport equation slover(Zhibin Ren 8-02-00)*******************
%*****************************************************************************

function [Fn_new,Ne_new,Ne_sub,E_sub]=charge(Ne_old,Ec_old,Ne_sub_old,E_sub_old)

%*************************GLOBAL VARIABLES************************************

global q k_B h_bar m_e m_t m_l Nc Ncc mx my mz Temp;
global transport_model
global dummy_flag fermi_flag ox_pnt_flag t_vall max_subband;
global delta_x Lsd Lg_top Lg_bot t_top t_bot t_si dx dy;
global Vs Vd;
global junction_l junction_r;
global Nx Ny Ntotal;
global div_avd mu_low beta Vel_sat criterion_outer;
global eta;
global N_dos E;
%Energy transport parameters

global ELE_TAUW ELE_CQ delta_T_1 delta_T_2 dim_c;

%NEGF scattering method parameters

global nu_scatter Info_scatter_new Info_scatter_old;
global Is Id;

Lsda=round(Lsd/dx);
Lg_topa=round(Lg_top/dx);
Lg_bota=round(Lg_bot/dx);
t_topa=round(t_top/dy);
t_bota=round(t_bot/dy);
t_sia=round(t_si/dy);
%**********************MODEL 4 related parameters*****************************
eta=1e-6*i;
Ef_tail_low=0.01;
Ef_tail_up=0.3;
E_step=criterion_outer/2;%0.5 times criterion_outer 
zeta_self=-i*(25*120/(mu_low*1e4))*1.0e-3;%eV
decay_fac=1;%dimensionless
%%%%********************MODEL 3 related parameters****************************
%%%%********************MODEL 2 related parameters****************************
e_field_limit=1.0;%in V/m
%*****************************************************************************
%**************************INPUT AND OUTPUT VARIABLES*************************
%*****************************************************************************
%Ne_old is the old electron density 
%1 column of Ntotal elements
%Ec_old is the old potential energy profile in eV
%1 column of Ntotal elements
%Te_old is the old electron temperature along the channel
%1 column of Nx elements, for charge-sheet-model

%Ne_sub_old contains 2D electron density on each subband
%E_sub_old contains the 1D subband energy for all subband 
%considered
%Te_sub_old contains temperture data of electron on each subband
%Nx by m, for bulk model

%Fn_new 1 column of Ntotal elements, dummy varible 
%for Poisson Equation solving
%Ne_new is the new 3D density in m^-3
%1 column of Ntotal elements
%Ne_sub contains 2D electron density on each subband
%E_sub contains the 1D subband energy for all subband considered

%Fn is dummy Fermi energy level used by models 1,4 and 5 for 
%current calculation
%max_subband is the number of subbands considered
%t_vall is the number of ladders considered
%****************************************************************************

%************************INITIALIZATION**************************************
Ec_old=real(Ec_old);%1 column of Ntotal elements
Ne_old=real(Ne_old);%1 column of Ntotal elements
%Ne_sub_old=zeros(Nx,max_subband,t_vall);
%E_sub_old=zeros(Nx,max_subband,t_vall);

Fn_new=zeros(Ntotal,1);%1 column of Ntotal elements
Ne_new=zeros(Ntotal,1);%1 column of Ntotal elements
Ne_sub=zeros(Nx,max_subband,t_vall);
E_sub=zeros(Nx,max_subband,t_vall);

% MODIFIED BY RAMESH, JAN 13th 03. Comment out as
% already initialized in main.m
%Info_scatter_new=zeros(nu_scatter,4);
%Info_scatter_new=Info_scatter_old;

%temporarily used variables
if ox_pnt_flag==0
   Np_v=round(t_si/dy)+1;
elseif ox_pnt_flag==1
   Np_v=Ny;
end
Fn=zeros(Nx,1);
Ns=zeros(Nx,1);
N_body=zeros(Np_v,Nx);
N_body_sum=zeros(Np_v,Nx);
U_sub=zeros(Nx,max_subband,t_vall);
W_sub=zeros(Np_v,Nx,max_subband,t_vall);
U_bias=zeros(Nx,1);

%***************************************************************************
%THE START OF TRANSPORT EQUATION PART***************************************
%***************************************************************************

%***************************************************************************
 if transport_model==1%TRANSPORT MODEL 1************************************
%***************************************************************************
%INITIALIZE THE REAL FERMI ENERGY LEVEL AND COMPUTE
N_row_mid=round((Ntotal-Nx*(t_topa+t_bota+1))/Nx/2)-1;
Nc_start=(Nx*(t_topa+1))+N_row_mid*Nx+1;
Nc_end=(Nx*(t_topa+1))+(N_row_mid+1)*Nx;
U_bias=Ec_old(Nc_start:Nc_end);
Ez=(h_bar*pi/t_si)^2/(2*mz(1)*m_e)/q;

channel_s=junction_l;
channel_e=junction_r;
Fn_slope=(-Vd-(-Vs))/(channel_e-channel_s);

for i_node=1:Nx
  if i_node<=channel_s
    Fn(i_node)=-Vs-Ez;
  elseif i_node>=channel_e
    Fn(i_node)=-Vd-Ez;
  else
    Fn(i_node)=-Vs-Ez+(i_node-channel_s)*Fn_slope;
  end
  if fermi_flag==1
    Ns(i_node)=Ncc...
              *log(1+exp((Fn(i_node)-U_bias(i_node))/(k_B*Temp/q)));
  elseif fermi_flag==0
    Ns(i_node)=Ncc...
              *exp((Fn(i_node)-U_bias(i_node))/(k_B*Temp/q));
  end
end

for iii_row=(Nx*(t_topa+1))/Nx:(Ntotal-Nx*t_bota)/Nx-2
   for iii_col=1:Nx
    i_node=iii_row*Nx+iii_col;
    Ne_new(i_node)=(sin((iii_row+1-(Nx*(t_topa+1))/Nx)...
               *dy/t_si*pi)^2)...
               *Ns(iii_col)*2/t_si;
   end
end

%****************************************************************************
 elseif transport_model==2 % DRIFT DIFFUSION ********************************
%****************************************************************************
mobility=zeros(nu_scatter-1,1);
mobility=(diff(Info_scatter_old(:,4))+2*Info_scatter_old(1:Nx-1,4))/2;

[U_sub,W_sub]=schred(Ec_old);

E_sub=U_sub;

for i_val=1:t_vall
   Nccc=2*m_e*sqrt(mx(i_val)*my(i_val))*k_B*Temp/(pi*h_bar^2);

   for i_sub=1:max_subband
      U_bias=U_sub(:,i_sub,i_val);
      Ns=Ne_sub_old(:,i_sub,i_val);
      if fermi_flag==1
        Ns_con=Nccc*log(1+exp((-Vs-U_bias(1))/(k_B*Temp/q)));
        Nd_con=Nccc*log(1+exp((-Vd-U_bias(Nx))/(k_B*Temp/q)));
      elseif fermi_flag==0
        Ns_con=Nccc*exp((-Vs-U_bias(1))/(k_B*Temp/q));
        Nd_con=Nccc*exp((-Vd-U_bias(Nx))/(k_B*Temp/q));
      end
      
      if fermi_flag==0
        fac_deg=ones(Nx-1,1);
      elseif fermi_flag==1
        arg=anti_dummy(Ns./Nccc,0,fermi_flag);
        Fn=U_bias+(k_B*Temp/q)*arg;
        zeta_channel=(Fn-U_bias)./(k_B*Temp/q);
        fac_deg_tem=log(1+exp(zeta_channel)).*(1+exp(-zeta_channel));
        fac_deg_tem(find(arg<-5))=1;% To ensure could scalability when subbands are empty
        fac_deg=(diff(fac_deg_tem)+2*fac_deg_tem(1:Nx-1))/2;
      end
      V_channel=-U_bias;
      E_channel=-diff(V_channel)/dx;
      
      %If you want constant mobility uncomment the following line
      %mu_channel=mu_low./(1+(mu_low/Vel_sat*abs(E_channel)).^beta).^(1/beta);
      
      %If you want doping dependent mobility uncomment the following line
      mu_channel=mobility./(1+(mobility./Vel_sat.*abs(E_channel)).^beta).^(1/beta);
      
      for j=1:Nx-1
        if abs(E_channel(j))<=abs(e_field_limit*fac_deg(j))
           Coe_1(j)=-mu_channel(j)*(k_B*Temp/q)*fac_deg(j)/dx;
           Coe_2(j)=-Coe_1(j);
        else
           Coe_1(j)=mu_channel(j)*E_channel(j)*1/(1-exp(E_channel(j)*dx/...
                  ((k_B*Temp/q)*fac_deg(j))));
           Coe_2(j)=mu_channel(j)*E_channel(j)*1/(1-exp(-E_channel(j)*dx/...
                  ((k_B*Temp/q)*fac_deg(j))));
        end
      end

      AAA=diag(-Coe_1(2:Nx-1)+Coe_2(1:Nx-2))...
         -diag(Coe_2(2:Nx-2),1)+diag(Coe_1(2:Nx-2),-1);
      CCC=zeros(Nx-2,1);CCC(1)=-Coe_1(1)*Ns_con/Ncc;
      CCC(Nx-2)=Coe_2(Nx-1)*Nd_con/Ncc;

      BBB=sparse(AAA)\sparse(CCC);
      Ne_sub(:,i_sub,i_val)=[Ns_con;full(BBB)*Ncc;Nd_con];

      for i_node=1:Nx
      	N_body(:,i_node)=Ne_sub(i_node,i_sub,i_val)...
            *W_sub(:,i_node,i_sub,i_val)/dy;
      end
   
      N_body_sum=N_body_sum+N_body;
    end
end

Info_scatter_new(:,3)=Fn;

if ox_pnt_flag==0   
   Ne_new=[Ne_old(1:(Nx*(t_topa+1)));...
      reshape((N_body_sum(2:Np_v-1,:))',...
      Nx*(Np_v-2),1);...
      Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)];
elseif ox_pnt_flag==1
   Ne_new=reshape(N_body_sum',Ntotal,1);
end

%*******************************************************************************
elseif transport_model==3 %BALLISTIC TRANSPORT USING SEMICLASSICAL APPROACH*****
%*******************************************************************************

[U_sub,W_sub]=schred(Ec_old);

E_sub=U_sub;

for i_val=1:t_vall
   Ne_2d_1=2*(sqrt(mx(i_val)*my(i_val))*m_e*k_B*Temp)/(2*pi*h_bar^2);
   Ne_2d_2=Ne_2d_1*2/pi^0.5;

   for i_sub=1:max_subband
    U_bias=U_sub(:,i_sub,i_val);
    [Ec_peak,i_peak]=max(U_bias);
    for i_node=1:Nx
      MEc_peak=(Ec_peak-U_bias(i_node))/(k_B*Temp/q);
      zeta_s=(-Vs-U_bias(i_node))/(k_B*Temp/q);
      zeta_d=(-Vd-U_bias(i_node))/(k_B*Temp/q);
      if zeta_s==zeta_d+10000
        Ne_sub(i_node,i_sub,i_val)=2*Ne_2d_1*fermi(zeta_s,fermi_flag,0);
      else   
        if i_node<=i_peak
            Ne_sub(i_node,i_sub,i_val)=...
            Ne_2d_2*integral(zeta_s,zeta_d,MEc_peak,fermi_flag);
   	    elseif i_node>i_peak
            Ne_sub(i_node,i_sub,i_val)=...
            Ne_2d_2*integral(zeta_d,zeta_s,MEc_peak,fermi_flag);
        end
      end     
      N_body(:,i_node)=Ne_sub(i_node,i_sub,i_val)...
                       *W_sub(:,i_node,i_sub,i_val)/dy;
    end
    N_body_sum=N_body_sum+N_body;
   end
end

if ox_pnt_flag==0   
   Ne_new=[Ne_old(1:(Nx*(t_topa+1)));...
      reshape((N_body_sum(2:Np_v-1,:))',...
      Nx*(Np_v-2),1);...
      Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)];
elseif ox_pnt_flag==1
   Ne_new=reshape(N_body_sum',Ntotal,1);
end
 
%*************************************************************************
elseif transport_model==4%BALLISTIC TRANSPORT MODEL USING GREEN FUNCTION APPROACH
%*************************************************************************

[U_sub,W_sub]=schred(Ec_old);

E_sub=U_sub;

%E , E_number are defined outside the loop to remain constant 
% When simulating several valleys and subband so the size of the DOS matrix
% Remain constant...This maybe more PCU intensive...

[Ec_peak,i_peak]=max(U_sub(:,max_subband,t_vall));
Ec_peak=max(-Vs,Ec_peak);
E_number=round((Ec_peak+Ef_tail_up-U_sub(Nx,1,1)+Ef_tail_low)/E_step)+2
E=linspace((U_sub(Nx,1,1)-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number);
delta_E=((Ec_peak+Ef_tail_up)-(U_sub(Nx,1,1)-Ef_tail_low))/(E_number-1);

for i_val=1:t_vall
   Ne_2d=2*sqrt(mx(i_val)*m_e*(k_B*Temp/q)*q/(2*pi^3))/(h_bar*dx);
   tt=(h_bar^2)/(2*my(i_val)*m_e*(dx^2)*q);
   A=tt*((2*eye(Nx))-(diag(ones(Nx-1,1),1))-(diag(ones(Nx-1,1),-1)));

   for i_sub=1:max_subband
      U_bias=U_sub(:,i_sub,i_val);
      %[Ec_peak,i_peak]=max(U_bias);
      %Ec_peak=max(-Vs,max(U_bias));
      %E_number=round((Ec_peak+Ef_tail_up-U_bias(Nx)+Ef_tail_low)/E_step)+2
      %E=linspace((U_bias(Nx)-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number);
      %delta_E=((Ec_peak+Ef_tail_up)-(U_bias(Nx)-Ef_tail_low))/(E_number-1);
      N_den=zeros(Nx,1);
      B_s=zeros(Nx,1);
      B_d=zeros(Nx,1);
      B_s(1)=1;
      B_d(Nx)=1;
      spB_s=sparse(B_s);
      spB_d=sparse(B_d);
    
      %for k=1:E_number,
        %ee=E(k);ep=ee+eta;
        %ck=1-((ep-U_bias(1))/(2*tt));con_s=-tt*exp(i*acos(ck));
        %ck=1-((ep-U_bias(Nx))/(2*tt));con_d=-tt*exp(i*acos(ck));
        %U_eff=U_bias;
        %U_eff(1)=U_bias(1)+con_s;
        %U_eff(Nx)=U_bias(Nx)+con_d;
        %G_inv=sparse((ep*eye(Nx))-A-diag(U_eff));
        %G_s=G_inv\spB_s;
        %G_d=G_inv\spB_d;
        %f_1=fermi(((-Vs-ee)/(k_B*Temp/q)),fermi_flag,-1/2);
        %f_2=fermi(((-Vd-ee)/(k_B*Temp/q)),fermi_flag,-1/2);
        %N_den=N_den-abs(G_s).^2*imag(con_s)*2*f_1...
        %      -abs(G_d).^2*imag(con_d)*2*f_2;
        %N_dos_one(:,k)=-abs(G_s).^2*imag(con_s)-abs(G_d).^2*imag(con_d);     
      %end
      
      [N_den1 Nquad]=myquad(@func_energy,E(1),E(E_number),1e-6,[],tt,U_bias,A,spB_s,spB_d);
      N_den=N_den1/(E(2)-E(1));
      
      Ne_sub(:,i_sub,i_val)=full(N_den)*Ne_2d*delta_E;
      for i_node=1:Nx
      	N_body(:,i_node)=Ne_sub(i_node,i_sub,i_val)...
            *W_sub(:,i_node,i_sub,i_val)/dy;
      end
   
      N_body_sum=N_body_sum+N_body;
	end
end

if ox_pnt_flag==0   
   Ne_new=[Ne_old(1:(Nx*(t_topa+1)));...
      reshape((N_body_sum(2:Np_v-1,:))',...
      Nx*(Np_v-2),1);...
      Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)];
elseif ox_pnt_flag==1
   Ne_new=reshape(N_body_sum',Ntotal,1);
end

%*******************************************************************************
elseif transport_model==5 %SCATTERING MODEL USING GREEN"S FUNCTION METHOD*******
%*******************************************************************************

[U_sub,W_sub]=schred(Ec_old);
E_sub=U_sub;

U_scatter=zeros(Nx,1);
mu_scatter_old=zeros(nu_scatter+2,1);
mu_scatter_new=zeros(nu_scatter+2,1);
delta_mu=zeros(nu_scatter+2,1);
i_scatter=zeros(nu_scatter+2,1);
BB_dummy=eye(Nx);
Iin=zeros(nu_scatter,1);

%comment the two following line for constant mobility
%zeta_self=zeros(nu_scatter,1);
%decay_fac=1;%dimensionless

for i_s=1:nu_scatter
 i_scatter(i_s)=Info_scatter_old(i_s,2);

 % CHANGED BY RAMESH TO TAKE OLD GUESS
 mu_scatter_old(i_s)=Info_scatter_new(i_s,3);
end

mu_scatter_old(nu_scatter+1)=-Vs;
mu_scatter_old(nu_scatter+2)=-Vd;
i_scatter(nu_scatter+1)=1;
i_scatter(nu_scatter+2)=Nx;
BB=sparse(Nx,nu_scatter+2);

BB=[BB_dummy(:,i_scatter(1:nu_scatter)),...
    BB_dummy(:,1),BB_dummy(:,Nx)];

Ec_peak=max(-Vs,max(U_sub(:,1,1)));
Ec_bottom=U_sub(Nx,1,1);
E_number=round((Ec_peak+Ef_tail_up-Ec_bottom+Ef_tail_low)/E_step)+2;
E=linspace((Ec_bottom-Ef_tail_low),(Ec_peak+Ef_tail_up),E_number);
delta_E=((Ec_peak+Ef_tail_up)-(Ec_bottom-Ef_tail_low))/(E_number-1);

GG=sparse(Nx,nu_scatter+2);
T_E=zeros(E_number,nu_scatter+2,nu_scatter+2);

%comment the next line and uncomment the previous line for constant mobility

E_self=ones(E_number,nu_scatter+2);%correction on March 14th to include mobility
U_tem=zeros(nu_scatter,1);

%end of initialization

% ADDED BY RAMESH JAN 13th 03.
% We need to compute a mean free path to set
% zeta_self based on the Fermi-level of the
% probe, subband energy and grid spacing.
sum1 = zeros(nu_scatter,1);
sum2 = zeros(nu_scatter,1);
lambda = zeros(nu_scatter,1);
for i_val=1:t_vall
    for i_sub=1:max_subband
	for i_s=1:nu_scatter
	    factor = exp((Info_scatter_new(i_s,3)-U_sub(i_scatter(i_s),i_sub,i_val))...
		 /(k_B*Temp/q));
	    sum1(i_s) = sum1(i_s)...
		 +my(i_val)*m_e*log(1+factor);
            % Prevent divide by zero
	    factor1 = log(1+factor);
	    if(factor1<=1e-20)
	       factor1 =  factor;
            end 
		 factor2 = factor/(factor1*(1+factor))*fermi(((Info_scatter_new(i_s,3)-...
                 U_sub(i_scatter(i_s),i_sub,i_val))/(k_B*Temp/q)),fermi_flag,1/2);
	    sum2(i_s) = sum2(i_s)...
                 +my(i_val)*m_e*sqrt(2*k_B*Temp/(pi*mx(i_val)*m_e))*factor2;		 
        end
    end
end
for i_s=1:nu_scatter
      lambda(i_s)=2*k_B*Temp/q*Info_scatter_old(i_s,4)*...
      sum1(i_s)/sum2(i_s);
end

for i_val=1:t_vall
   Ie_2d=2*q^2/(pi^2*h_bar^2)*sqrt(my(i_val)*m_e*(k_B*Temp/q)*q*pi/2);
   tt=(h_bar^2)/(2*mx(i_val)*m_e*(dx^2)*q);
   A=tt*((2*eye(Nx))-(diag(ones(Nx-1,1),1))-(diag(ones(Nx-1,1),-1)));

   for i_sub=1:max_subband
      U_bias=U_sub(:,i_sub,i_val);
      U_tem=U_bias(i_scatter(1:nu_scatter));

      for k=1:E_number,
        ee=E(k);ep=ee+eta;
        ck=1-((ep-U_bias(1))/(2*tt));con_s=-tt*exp(i*acos(ck));
        ck=1-((ep-U_bias(Nx))/(2*tt));con_d=-tt*exp(i*acos(ck));
        
        % ADDED BY RAMESH JAN 13 03.
	for i_s=1:nu_scatter
	    elutot = ee-(A(i_scatter(i_s),i_scatter(i_s))+U_bias(i_scatter(i_s)));
	    g1 = (elutot+sqrt(elutot^2-4*tt^2))/(2*tt^2);
            g2 = (elutot-sqrt(elutot^2-4*tt^2))/(2*tt^2);
	    g3 = (imag(g2)<=0)*g2+(imag(g1<=0))*g1;
	    E_self(k,i_s) = i*imag(g3*2*dx*tt^2/lambda(i_s));
        end
        
        E_self(k,nu_scatter+1)=con_s;
        E_self(k,nu_scatter+2)=con_d;
        U_scatter(i_scatter)=transpose(E_self(k,:));
               
        U_eff=U_bias+U_scatter;
        G_inv=sparse((ep*eye(Nx))-A-diag(U_eff));
        GG=G_inv\BB;
        T_E(k,:,:)=4*Ie_2d*delta_E*...
        spdiags(imag(E_self(k,:))',0,nu_scatter+2,nu_scatter+2)*...
        abs(GG(i_scatter,:)).^2*...
        spdiags(imag(E_self(k,:))',0,nu_scatter+2,nu_scatter+2)+squeeze(T_E(k,:,:));
      end

    end%subband
end%valley

%calculate current, and search for mu_scatter
%based on current continuity

[Is,Id,Iin,mu_scatter_new]=current_mat(mu_scatter_old,T_E,E);

%end of mu_scatter search

%Calculate charge density
for i_val=1:t_vall
   Ne_2d=2*sqrt(my(i_val)*m_e*(k_B*Temp/q)*q/(2*pi^3))/(h_bar*dx);
   tt=(h_bar^2)/(2*mx(i_val)*m_e*(dx^2)*q);
   A=tt*((2*eye(Nx))-(diag(ones(Nx-1,1),1))-(diag(ones(Nx-1,1),-1)));

   for i_sub=1:max_subband
      U_bias=U_sub(:,i_sub,i_val);
      U_tem=U_bias(i_scatter(1:nu_scatter));
     
      for k=1:E_number,
        
        % All this needs to be recomputed for multiple bands
        %--------------------------------------------------
        ee=E(k);ep=ee+eta;
        ck=1-((ep-U_bias(1))/(2*tt));con_s=-tt*exp(i*acos(ck));
        ck=1-((ep-U_bias(Nx))/(2*tt));con_d=-tt*exp(i*acos(ck));
             
         % ADDED BY RAMESH JAN 13 03.
	for i_s=1:nu_scatter
	    elutot = ee-(A(i_scatter(i_s),i_scatter(i_s))+U_bias(i_scatter(i_s)));
	    g1 = (elutot+sqrt(elutot^2-4*tt^2))/(2*tt^2);
            g2 = (elutot-sqrt(elutot^2-4*tt^2))/(2*tt^2);
	    g3 = (imag(g2)<=0)*g2+(imag(g1<=0))*g1;
	    E_self(k,i_s) = i*imag(g3*2*dx*tt^2/lambda(i_s));
        end
        
        E_self(k,nu_scatter+1)=con_s;
        E_self(k,nu_scatter+2)=con_d;
        U_scatter(i_scatter)=transpose(E_self(k,:));
                
        U_eff=U_bias+U_scatter;
        G_inv=sparse((ep*eye(Nx))-A-diag(U_eff));
        GG=G_inv\BB;
        Ne_sub(:,i_sub,i_val)=Ne_sub(:,i_sub,i_val)...
        -2*Ne_2d*delta_E*abs(GG).^2*(imag(E_self(k,:))'.*...
        fermi(((mu_scatter_new-ee)/(k_B*Temp/q)),fermi_flag,-1/2));
        
        GGG = zeros(Nx,Nx);
        GGG(:,1)=GG(:,Nx-1);
        GGG(:,Nx)=GG(:,Nx);
        GGG(1:Nx,2:Nx-1)=GG(1:Nx,1:Nx-2);
        N_dos(:,k)=-2*imag(diag(GGG));
        
      end
      
      for i_node=1:Nx
         N_body(:,i_node)=Ne_sub(i_node,i_sub,i_val)...
            *W_sub(:,i_node,i_sub,i_val)/dy;
      end
      N_body_sum=N_body_sum+N_body;

    end%subband
    
end%valley

Info_scatter_new(:,1)=Iin;
Info_scatter_new(1,1)=Is;
Info_scatter_new(:,3)=mu_scatter_new(1:nu_scatter);

if ox_pnt_flag==0   
   Ne_new=[Ne_old(1:(Nx*(t_topa+1)));...
      reshape((N_body_sum(2:Np_v-1,:))',...
      Nx*(Np_v-2),1);...
      Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)];
elseif ox_pnt_flag==1
   Ne_new=reshape(N_body_sum',Ntotal,1);
end

%*******************************************************************************
elseif transport_model==6 %ENERGY TRANSPORT MODEL*******************************
%*******************************************************************************

[U_sub,W_sub]=schred(Ec_old);
E_sub=U_sub;

for i_val=1:t_vall
   Nccc=2*m_e*sqrt(mx(i_val)*my(i_val))*k_B*Temp/(pi*h_bar^2);

   for i_sub=1:max_subband
      U_bias=U_sub(:,i_sub,i_val);
      Ns=Ne_sub_old(:,i_sub,i_val);

    %===============================
    %  Nondegenerate always for ET
    %===============================
      Ns_con=Nccc*exp((-Vs-U_bias(1))/(k_B*Temp/q));
      Nd_con=Nccc*exp((-Vd-U_bias(Nx))/(k_B*Temp/q));
      
      V_channel=-U_bias;

      Ne_sub(:,i_sub,i_val)=ET(V_channel,mu_low,Ns_con,Nd_con,...
				Temp,Temp,ELE_TAUW,ELE_CQ,Vel_sat,beta);

      for i_node=1:Nx
      	N_body(:,i_node)=Ne_sub(i_node,i_sub,i_val)...
            *W_sub(:,i_node,i_sub,i_val)/dy;
      end
   
      N_body_sum=N_body_sum+N_body;
      
    end
end

if ox_pnt_flag==0   
   Ne_new=[Ne_old(1:(Nx*(t_topa+1)));...
      reshape((N_body_sum(2:Np_v-1,:))',...
      Nx*(Np_v-2),1);...
      Ne_old((Ntotal-Nx*(t_bota+1)+1):Ntotal)];
elseif ox_pnt_flag==1
   Ne_new=reshape(N_body_sum',Ntotal,1);
end

%*******************************************************************************
%***************THE END OF TRANSPORT EQUATION PART******************************
%*******************************************************************************
end
%*******************************************************************************
%*********************START OF VARIABLE CHANGE PART*****************************
%*******************************************************************************

 if ox_pnt_flag==0%No electron penetration into oxide
    
   for iii_row=(Nx*(t_topa+1))/Nx:(Ntotal-Nx*t_bota)/Nx-2
    for iii_col=1:Nx
      i_node=iii_row*Nx+iii_col;
      Fn_new(i_node)=anti_dummy(Ne_new(i_node)/Nc,dummy_flag,fermi_flag)...
                     *(k_B*Temp/q)+Ec_old(i_node);              
    end
   end
   
elseif ox_pnt_flag==1%assuming electron penetration into oxide

   Fn_new=anti_dummy((Ne_new+div_avd)./Nc,dummy_flag,fermi_flag)...
          *(k_B*Temp/q)+Ec_old;              
       
 end

%****************************************************************************
%**********************END OF VARIABLE CHANGE PART***************************
%****************************************************************************

%****************************************************************************
%***********************END OF FUNCTION CHARGE.M*****************************
%****************************************************************************
