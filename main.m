%***************************************************************************************
%******************2D_Poisson+1D_Schroedinger+1D_transport******************************
%***************************************************************************************

function [Ie,Ie_sub_body,Te_sub_body,Ne_sub_body,E_sub_body,Ne_3d,Ec_3d,conv]=main


%***************************INPUT VARIABLES**********************************************
global transport_model;
global DG_flag fermi_flag ox_pnt_flag dummy_flag;
global t_vall max_subband Temp;
global dopslope_s dopslope_d overlap_s overlap_d;
global Lsd Lg_top Lg_bot t_top t_bot t_si dx dy refine Te;
global Vg1 Vg2 Vs Vd Vd_initial Vg_step Ng_step Vd_step Nd_step;
global phi_top phi_bot eps_top eps_bot bar_top bar_bot eps_si m_t m_l;
global mu_low beta Vel_sat charge_fac;
global N_sd N_body criterion_outer criterion_inner;
global junction_l junction_r;
global Eg1 Eg2 Es Ed Vd_temp;
global Nx Ny Ntotal;

%Energy model parameters

global ELE_TAUW ELE_CQ delta_T_1 delta_T_2 dim_c;

%NEGF scattering model parameters

global nu_scatter Info_scatter_new Info_scatter_old;
global Is Id;
global N_dos Trans E;

%************************FUNDAMENTAL physical constants**********************************
global eps_o psi_si q k_B h_bar m_e ...
       Nc Ncc mx my mz Temp;
global F_prime Nd;
global div_avd

Lsd=round(Lsd/dx)*dx;
Lg_top=round(Lg_top/dx)*dx;
Lg_bot=round(Lg_bot/dx)*dx;
t_top=round(t_top/dy)*dy;
t_bot=round(t_bot/dy)*dy;
t_si=round(t_si/dy)*dy;

Lsda=round(Lsd/dx);
Lg_topa=round(Lg_top/dx);
Lg_bota=round(Lg_bot/dx);
t_topa=round(t_top/dy);
t_bota=round(t_bot/dy);
t_sia=round(t_si/dy);

%Parameters for ET model
delta_T_1=5/2;%energy flux parameter one
delta_T_2=5/2;%energy flux parameter two
dim_c=2;%degree of carrier freedom 

%*****************************************************************************************
%***************************Step FUNCTION profile for Nsd*********************************
%*****************************************************************************************
junction_l=round((Lsd+overlap_s)/dx)+1;
junction_r=round((Lsd+Lg_top-overlap_d)/dx)+1;
%*****************************************************************************************
mx=zeros(3,1);my=zeros(3,1);mz=zeros(3,1);
Temp=Te;
eps_o=8.85e-12;
psi_si=4.05;
q=1.6e-19;
k_B=1.38e-23;
h_bar=1.05e-34;

m_e=0.91e-30;
Nc=2.8e25;
Ncc=2*m_e*m_t*k_B*Temp/(pi*h_bar^2);

mx(1)=m_t;mx(2)=m_l;mx(3)=m_t;
my(1)=m_t;my(2)=m_t;my(3)=m_l;
mz(1)=m_l;mz(2)=m_t;mz(3)=m_t;
%****************************************************************************************
%SPECIFY THE NEUTRAL BOUNDARY ***********************************************************
%Calculate boundary Ec based neutral charge and Fermi-Dirac statistics*******************
%****************************************************************************************
if ox_pnt_flag==0
  N_sd=((t_si/dy)/(t_si/dy-1))*N_sd;
  N_body=((t_si/dy)/(t_si/dy-1))*N_body;
end
  Eg1=-Vg1+phi_top-psi_si;
  Eg2=-Vg2+phi_bot-psi_si;

    if fermi_flag==0
       Es=-Vs-k_B*Temp/q*log((N_sd-N_body)/Ncc);
       Ed=-Vd-k_B*Temp/q*log((N_sd-N_body)/Ncc);
   elseif fermi_flag==1
       Es=-Vs-k_B*Temp/q*log(exp((N_sd-N_body)/Ncc)-1);
       Ed=-Vd-k_B*Temp/q*log(exp((N_sd-N_body)/Ncc)-1);
   end

%****************************************************************************************
%*************************END OF SPECIFY THE NEUTRAL BOUNDARY**************************** 
%****************************************************************************************

%*********************************ASSIGNING VARIABLES************************************

%NORMALIZATION PARAMETERS
charge_fac=dx*dy*q/(eps_si*eps_o)*Nc;
div_avd=1e-10*Nc;%a parameter used to avoid 
%divergence in converting electron density to dummy quantity

Nx=round((2*Lsd+Lg_top)/dx)+1;
Ny=round((t_top+t_bot+t_si)/dy)+1;
Ntotal=Nx*Ny;

%******************************************************************************
%*********************END OF ASSIGNING VARIABLES*******************************
%******************************************************************************

%******************************************************************************
%****************************START OF INITIALIZATION***************************
%******************************************************************************
Nd=zeros(Ntotal,1);%unchanged through the entire calculation
F_prime=zeros(Ntotal,Ntotal);%unchanged through the entire 
                               %calculation
                               
Ne_old=zeros(Ntotal,1);
Ne_new=zeros(Ntotal,1);
Ec_old=zeros(Ntotal,1);
Ec_new=zeros(Ntotal,1);

Fn_new=zeros(Ntotal,1);

Ne_sub=zeros(Nx,max_subband,t_vall);
E_sub=zeros(Nx,max_subband,t_vall);
Ne_sub_old=zeros(Nx,max_subband,t_vall);
E_sub_old=zeros(Nx,max_subband,t_vall);

Ie_tem=zeros(Nx,1);
Ie_sub=zeros(Nx,max_subband,t_vall);
Mu_sub=zeros(Nx,max_subband,t_vall);
Te_sub=zeros(Nx,max_subband,t_vall);

%***************************START OF SPECIFYING Nd****************************
doping;
%**************************END OF SPECIFING Nd*******************************

%******************Preparing F_prime(one time evaluation)********************
fprime;
%**************************END OF SPECIFIING F_prime*************************

%****************************************************************************
%****************************END OF INITIALIZATION***************************
%****************************************************************************

%****************************************************************************
%************START OF SELF CONSISTENT CALCULATION OF POISSON AND ************
%****************************TRANSPORT EQUATIONS*****************************
%****************************************************************************
if transport_model==5
nu_scatter=Nx-2;
elseif transport_model==2
nu_scatter=Nx;
else
nu_scatter=1;
end

Info_scatter_old=zeros(nu_scatter,4);
Info_scatter_new=zeros(nu_scatter,4);

%see reference, MEDICI manual, p2-15
mu_min=55*1e-4;
mu_max=300*1e-4;
Nref=1e22;
alpha=0.73;

%============Modified. Mar 18, 2002==================
Nd2D=reshape(Nd,Nx,Ny);
%============Modified. Mar 18, 2002==================

for i=1:nu_scatter
 Info_scatter_old(i,2)=i+1;
 %Info_scatter_old(i,4)=1/(1/mu_low+...
 %1/(mu_min+(mu_max-mu_min)./(1+(abs(Nd(Nx*round(t_top/dy)+1+Nx+i))/Nref).^alpha)));
 %Info_scatter_old(i,4)=1/(1/mu_low+1/(mu_min+(mu_max-mu_min)./(1+(abs(Nd2D(i,round(Ny/2)))/Nref).^alpha)));
 %============No Methiessen's rule========================================================
 Info_scatter_old(i,4)=mu_min+(mu_low-mu_min)./(1+(abs(Nd2D(i,round(Ny/2)))/Nref).^alpha);
 %========================================================================================
end

%keyboard
%****************************compress matrix*********************************
  spEc=sparse(Ec_old);
  spNe=sparse(Ne_old);
  spNd=sparse(Nd);
  F_prime=sparse(F_prime);
  
%***********************START OF INITIAL GUESS ******************************
  trans_temp=transport_model;
  fermi_temp=fermi_flag;
  transport_model=1;
  fermi_flag=1;
  Vd_temp=Vd;
  Ed_temp=Ed;
  Ed=Ed_temp+Vd-Vd_initial;
  Vd=Vd_initial;

  [Fn_new,Ne_new,Ne_sub,E_sub]=charge(spNe,spEc,Ne_sub_old,E_sub_old);
  %Info_scatter_old=Info_scatter_new;
  spFn=sparse(Fn_new);
  spNe=sparse(Ne_new);
  Ne_sub_old=Ne_sub;
  E_sub_old=E_sub;
  
  [Ec_new]=poisson(spNd,spFn,spEc);

  spEc=sparse(Ec_new);                   
                      
  transport_model=trans_temp
  fermi_flag=fermi_temp
  Ntotal
  
  
   if ((transport_model~=3) & fermi_flag==1)
        
     transport_model=3;

    [Fn_new,Ne_new,Ne_sub,E_sub]=charge(spNe,spEc,Ne_sub_old,E_sub_old);
   %Info_scatter_old=Info_scatter_new;
   
     spFn=sparse(Fn_new);
     spNe=sparse(Ne_new);
     Ne_sub_old=Ne_sub;
     E_sub_old=E_sub;

     [Ec_new]=poisson(spNd,spFn,spEc);
 
     spEc=sparse(Ec_new);
     
     transport_model=trans_temp;
     iter_outer=0;
     error_outer=1;
     while(error_outer>=criterion_outer)

      [Fn_new,Ne_new,Ne_sub,E_sub]=charge(spNe,spEc,Ne_sub_old,E_sub_old);
     %Info_scatter_old=Info_scatter_new;
     
      spEc_old=spEc;
      spFn=sparse(Fn_new);
      spNe=sparse(Ne_new);
      Ne_sub_old=Ne_sub;
      E_sub_old=E_sub;

      [Ec_new]=poisson(spNd,spFn,spEc);

      spEc=sparse(Ec_new);
      iter_outer=iter_outer+1;
      error_outer=max(abs(full(real(spEc-spEc_old))));
      fprintf ('%s %e \n','error_outer = ',error_outer);
     end
   end
   
    SpNein=spNe;
    SpEcin=spEc;
    Ne_sub_oldin=Ne_sub_old;
    E_sub_oldin=E_sub_old;
    SpNdin=spNd;
    SpFnin=spFn; 
   
%***************************END OF INITIAL GUESS OF Ec******************************

%*************************START OF CURRENT CALCULATION LOOP*************************

%******************************GATE BIAS LOOP***************************************
%transport_model=trans_temp;
Eg1_temp=Eg1;
Eg2_temp=Eg2;
for ii_vg=1:Ng_step+1
  Vg_bias(ii_vg)=Vg1+Vg_step*(ii_vg-1);
  Eg1=Eg1_temp-Vg_step*(ii_vg-1);
  if DG_flag==1    
     Eg2=Eg2_temp-Vg_step*(ii_vg-1);
  end

    spNe=SpNein;
    spEc=SpEcin;
    Ne_sub_old=Ne_sub_oldin;
    E_sub_old=E_sub_oldin;
    spNd=SpNdin;
    spFn=SpFnin;    
%**********************************DRAIN BIAS LOOP**********************************

     for ii_vd=1:Nd_step+1
        Vd_bias(ii_vd)=Vd_temp+Vd_step*(ii_vd-1);
        Ed=Ed_temp-Vd_step*(ii_vd-1);
        Vd=Vd_bias(ii_vd);

%***************************START OF SELF CONSISTENT LOOP***************************
        iter_outer=0;
        error_outer=1;
        converge=[error_outer];
        
            while(error_outer>=criterion_outer)

            [Fn_new,Ne_new,Ne_sub,E_sub]=charge(spNe,spEc,Ne_sub_old,E_sub_old);
            %Info_scatter_old=Info_scatter_new;
            
            spEc_old=spEc;
            spFn=sparse(Fn_new);
            spNe=sparse(Ne_new);
            Ne_sub_old=Ne_sub;
            E_sub_old=E_sub;
            
            [Ec_new]=poisson(spNd,spFn,spEc);

            spEc=sparse(Ec_new);
            iter_outer=iter_outer+1
            error_outer=max(abs(full(real(spEc-spEc_old))))
            
            converge=[converge;[error_outer]];
            
                if iter_outer>50
                    ver='******Converge Problem!!! Please step down DVMAX******';
                    disp(ver)
                    break;
                end
            end
 
        if (transport_model==5)
            Ie_tem=Is;
        else
            [Ie_tem,Ie_sub,Te_sub,Mu_sub]=current(spNe,spEc,Ne_sub,E_sub);
        end
%*************************END OF SELF CONSISTENT LOOP******************************

  Vggg=Vg_bias(ii_vg)
  Vddd=Vd_bias(ii_vd)

  Ie(ii_vg,ii_vd)=mean(real(Ie_tem));
  Mu_sub_body(:,:,:,ii_vg,ii_vd)=full(Mu_sub);
  Ie_sub_body(:,:,:,ii_vg,ii_vd)=full(Ie_sub);
  Ne_sub_body(:,:,:,ii_vg,ii_vd)=full(Ne_sub);
  Te_sub_body(:,:,:,ii_vg,ii_vd)=full(Te_sub);
  E_sub_body(:,:,:,ii_vg,ii_vd)=full(E_sub);
  Ne_3d(:,ii_vg,ii_vd)=full(Ne_new);
  Ec_3d(:,ii_vg,ii_vd)=full(Ec_new);
  conv(1:length(converge),ii_vg,ii_vd)=converge;
   
end
%***************************END OF DRAIN BIAS LOOP*********************************
end
%***************************END OF GATE BIAS LOOP**********************************

%**********************************************************************************
%******************END OF SELF CONSISTENT CALCULATION OF POISSON AND***************
%*******************************TRANSPORT EQUATIONS********************************
%**********************************************************************************



