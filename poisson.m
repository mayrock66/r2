%****************************************************************************
%********2D Poisson equation solver using Newton approach(Zhibin Ren 7-18-00)
%****************************************************************************
function [Ec]=poisson(spNd,spFn,Ec_old)

global F_prime transport_model;
global dummy_flag fermi_flag ox_pnt_flag;
global Lsd Lg_top Lg_bot t_top t_bot t_si dx dy;
global eps_top eps_bot eps_si;
global div_avd charge_fac criterion_inner;
global eps_si q k_B Nc Temp;
global Eg1 Eg2 Es Ed;  
global Nx Ny Ntotal;

Lsda=round(Lsd/dx);
Lg_topa=round(Lg_top/dx);
Lg_bota=round(Lg_bot/dx);
t_topa=round(t_top/dy);
t_bota=round(t_bot/dy);
t_sia=round(t_si/dy);
%**************************INPUT AND OUTPUT VARIABLES************************
%spNd is the 3d net dopant density
%1 column of Ntotal element
%spFn is the 3d electron density converted quasi Fermi energy
%from the solutions of transport equation, 
%1 column of Ntotal elements
%Ec_old is the input old potential energy profile
%1 column of Ntotal elements
%F_prime is the derivative matrix of Ne=0 F (Ntotal X Ntotal)
%criterion_inner is the poisson solver convergence criterion in eV
%Ec is the new potential profile computed
%1 column of Ntotal elements
%All equations involved assume ISS UNITS

%******************************INITIALIZATION************************************
delta_Ec=zeros(Ntotal,1);
F=zeros(Ntotal,1);
MF_prime=zeros(Ntotal,Ntotal);
Ec=zeros(Ntotal,1);
Fn=zeros(Ntotal,1);
Nd=zeros(Ntotal,1);
dummy_fun=zeros(Ntotal,1);
dummy_fun_prime=zeros(Ntotal,1);

iter_inner=0;
error_inner=1;
Ec=real(Ec_old);
Nd=real(spNd);
Fn=real(spFn);

if ox_pnt_flag==0
   div_avdt=0;
elseif ox_pnt_flag==1
   div_avdt=div_avd;
end

%*************************START OF INNER LOOP**************************************

while(error_inner>=criterion_inner)
iter_inner=iter_inner+1;
fprintf ('%s %i \n','iter_inner = ',iter_inner);

%*******************THE START OF DUMMY VARIABLE DEFINITION**************************
dummy_fun=charge_fac*((Nd+div_avdt)/Nc-dummy((Fn-Ec)/(k_B*Temp/q)...
            ,dummy_flag,fermi_flag));
dummy_fun_prime=charge_fac*dummy_prime((Fn-Ec)/(k_B*Temp/q)...
            ,dummy_flag,fermi_flag)/(k_B*Temp/q);

%***********************THE END OF DUMMY VARIABLE DEFINITION************************

if ox_pnt_flag==0%(NO ELECTRON PENETRATION INTO OXIDE REGIONS)

%********************************EVALUATE F*****************************************

%****************************Top gate insulator region******************************
for i=1:Nx*(t_topa+1)
  if(i>=1 & i<=Lsda+1-1) 
    F(i)=Ec(i)-Ec(i+Nx);
  elseif(i>=Lsda+1 & i<=(Lsda+Lg_topa)+1)
    F(i)=Ec(i)-Eg1;
  elseif(i>=(Lsda+Lg_topa)+1+1 & i<=Nx)
    F(i)=Ec(i)-Ec(i+Nx);
  elseif(i>=(Nx*t_topa+1) & i<=Nx*(t_topa+1))
    F(i)=-1/8*(dy/dx)*eps_top/eps_si*Ec(i-Nx-1)...
              -(dx/dy-1/4/(dx/dy))*eps_top/eps_si*Ec(i-Nx)...
              -1/8*(dy/dx)*eps_top/eps_si*Ec(i-Nx+1)...
              -3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec(i-1)...
              +(dx/dy+3/4/(dx/dy))*(eps_top+eps_si)/eps_si*Ec(i)...
              -3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec(i+1)...
              -1/8*(dy/dx)*Ec(i+Nx-1)...
              -(dx/dy-1/4/(dx/dy))*Ec(i+Nx)...
              -1/8*(dy/dx)*Ec(i+Nx+1);
  else
    F(i)=-(dx/dy)*eps_top/eps_si*Ec(i-Nx)...
              -(dy/dx)*eps_top/eps_si*Ec(i-1)...
              +2*(dx/dy+dy/dx)*eps_top/eps_si*Ec(i)...
              -(dy/dx)*eps_top/eps_si*Ec(i+1)...
              -(dx/dy)*eps_top/eps_si*Ec(i+Nx);
  end
end

%************************Bottom gate insulator region******************************
for i=(Ntotal-Nx*(t_bota+1)+1):Ntotal
  if(i>=(Ntotal-Nx*(t_bota+1)+1) & i<=(Ntotal-Nx*t_bota))
    F(i)=-1/8*(dy/dx)*Ec(i-Nx-1)...
              -(dx/dy-1/4/(dx/dy))*Ec(i-Nx)...
              -1/8*(dy/dx)*Ec(i-Nx+1)...
              -3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec(i-1)...
              +(dx/dy+3/4/(dx/dy))*(eps_bot+eps_si)/eps_si*Ec(i)...
              -3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec(i+1)...
              -1/8*(dy/dx)*eps_bot/eps_si*Ec(i+Nx-1)...
              -(dx/dy-1/4/(dx/dy))*eps_bot/eps_si*Ec(i+Nx)...
              -1/8*(dy/dx)*eps_bot/eps_si*Ec(i+Nx+1);
  elseif(i>=(Ntotal-Nx+1) & i<=(Ntotal-Nx+1+Lsda)-1)
    F(i)=Ec(i)-Ec(i-Nx);
  elseif(i>=(Ntotal-Nx+1+Lsda) & i<=(Ntotal-Nx+1+Lsda+Lg_bota))
    F(i)=Ec(i)-Eg2;
  elseif(i>=(Ntotal-Nx+1+Lsda+Lg_bota)+1 & i<=Ntotal)
    F(i)=Ec(i)-Ec(i-Nx);
  else
    F(i)=-(dx/dy)*eps_bot/eps_si*Ec(i-Nx)...
              -(dy/dx)*eps_bot/eps_si*Ec(i-1)...
              +2*(dx/dy+dy/dx)*eps_bot/eps_si*Ec(i)...
              -(dy/dx)*eps_bot/eps_si*Ec(i+1)...
              -(dx/dy)*eps_bot/eps_si*Ec(i+Nx);
  end
end

%********************Specify the F matrix in the silicon film region****************
for i=Nx*(t_topa+1)+1:(Ntotal-Nx*(t_bota+1)+1)-1
  F(i)=-(dx/dy)*Ec(i-Nx)-(dy/dx)*Ec(i-1)+2*(dx/dy+dy/dx)*Ec(i)...
            +dummy_fun(i)...
            -(dy/dx)*Ec(i+1)-(dx/dy)*Ec(i+Nx);
end

%****************Modify the F matrix at the right and left boundaries***************
i_l=1;
i_r=Nx;
for j=1:Ny
  if(j==1)
      
    F(i_l)=2*Ec(i_l)-Ec(i_l+1)-Ec(i_l+Nx);
    F(i_r)=2*Ec(i_r)-Ec(i_r-1)-Ec(i_r+Nx);
    
elseif(j>1 & j<round(Nx*(t_topa+1)/Nx))
    
    F(i_l)=Ec(i_l)-Ec(i_l+1);
    F(i_r)=Ec(i_r)-Ec(i_r-1);  
    
  elseif(j>=round(Nx*(t_topa+1)/Nx) & j<=round(Ntotal-Nx*t_bota/Nx))
        
  F(i_l)=Ec(i_l)-Ec(i_l+1);
  F(i_r)=Ec(i_r)-Ec(i_r-1);
        
  elseif(j>round(Ntotal-Nx*t_bota/Nx) & j<Ny)
      
    F(i_l)=Ec(i_l)-Ec(i_l+1);
    F(i_r)=Ec(i_r)-Ec(i_r-1);

  elseif(j==Ny & ((Ntotal-Nx+1)<(Ntotal-Nx+1+Lsda)) &...
          ((Ntotal-Nx+1+Lsda+Lg_bota)<Ntotal))
      
    F(i_l)=2*Ec(i_l)-Ec(i_l+1)-Ec(i_l-Nx);
    F(i_r)=2*Ec(i_r)-Ec(i_r-1)-Ec(i_r-Nx);
 
  end
  i_l=1+j*Nx;
  i_r=(j+1)*Nx;
end

%*****************************END OF EVALUATING F**********************************

%******************************EVALUATE MF_prime***********************************
%MF_prime matrix in the silicon film region
for j_row=Nx*(t_topa+1)/Nx:(Ntotal-Nx*t_bota)/Nx-2
  for j_col=2:Nx-1
    ii=j_row*Nx+j_col;
    MF_prime(ii,ii)=dummy_fun_prime(ii);
  end
end
%**************END OF EVALUATION FOR NO PENETRATION INTO THE OXIDE*****************

elseif ox_pnt_flag==1 %(ACCOUNTING FOR ELECTRON PENETRATION INTO OXIDE REGIONS)

%********************************EVALUATE F****************************************

%*************************Top gate insulator region********************************
for i=1:Nx*(t_topa+1)
  if(i>=1 & i<=Lsda+1-1) 
    F(i)=Ec(i)-Ec(i+Nx);
  elseif(i>=Lsda+1 & i<=(Lsda+Lg_topa)+1)
    F(i)=Ec(i)-Eg1;
  elseif(i>=(Lsda+Lg_topa)+1+1 & i<=Nx)
    F(i)=Ec(i)-Ec(i+Nx);
  elseif(i>=(Nx*t_topa+1) & i<=Nx*(t_topa+1))
    F(i)=-1/8*(dy/dx)*eps_top/eps_si*Ec(i-Nx-1)...
              -(dx/dy-1/4/(dx/dy))*eps_top/eps_si*Ec(i-Nx)...
              -1/8*(dy/dx)*eps_top/eps_si*Ec(i-Nx+1)...
              -3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec(i-1)...
              +(dx/dy+3/4/(dx/dy))*(eps_top+eps_si)/eps_si*Ec(i)...
              +dummy_fun(i)...
              -3/8*(dy/dx)*(eps_top+eps_si)/eps_si*Ec(i+1)...
              -1/8*(dy/dx)*Ec(i+Nx-1)...
              -(dx/dy-1/4/(dx/dy))*Ec(i+Nx)...
              -1/8*(dy/dx)*Ec(i+Nx+1);
  else
    F(i)=-(dx/dy)*eps_top/eps_si*Ec(i-Nx)...
              -(dy/dx)*eps_top/eps_si*Ec(i-1)...
              +2*(dx/dy+dy/dx)*eps_top/eps_si*Ec(i)+dummy_fun(i)...
              -(dy/dx)*eps_top/eps_si*Ec(i+1)...
              -(dx/dy)*eps_top/eps_si*Ec(i+Nx);
  end
end

%***************************Bottom gate insulator region***************************
for i=(Ntotal-Nx*(t_bota+1)+1):Ntotal
  if(i>=(Ntotal-Nx*(t_bota+1)+1) & i<=Ntotal-Nx*t_bota)
    F(i)=-1/8*(dy/dx)*Ec(i-Nx-1)...
              -(dx/dy-1/4/(dx/dy))*Ec(i-Nx)...
              -1/8*(dy/dx)*Ec(i-Nx+1)...
              -3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec(i-1)...
              +(dx/dy+3/4/(dx/dy))*(eps_bot+eps_si)/eps_si*Ec(i)...
              +dummy_fun(i)...
              -3/8*(dy/dx)*(eps_bot+eps_si)/eps_si*Ec(i+1)...
              -1/8*(dy/dx)*eps_bot/eps_si*Ec(i+Nx-1)...
              -(dx/dy-1/4/(dx/dy))*eps_bot/eps_si*Ec(i+Nx)...
              -1/8*(dy/dx)*eps_bot/eps_si*Ec(i+Nx+1);
  elseif(i>=(Ntotal-Nx+1) & i<=(Ntotal-Nx+1+Lsda)-1)
    F(i)=Ec(i)-Ec(i-Nx);
  elseif(i>=(Ntotal-Nx+1+Lsda) & i<=(Ntotal-Nx+1+Lsda+Lg_bota))
    F(i)=Ec(i)-Eg2;
  elseif(i>=(Ntotal-Nx+1+Lsda+Lg_bota)+1 & i<=Ntotal)
    F(i)=Ec(i)-Ec(i-Nx);
  else
    F(i)=-(dx/dy)*eps_bot/eps_si*Ec(i-Nx)...
              -(dy/dx)*eps_bot/eps_si*Ec(i-1)...
              +2*(dx/dy+dy/dx)*eps_bot/eps_si*Ec(i)+dummy_fun(i)...
              -(dy/dx)*eps_bot/eps_si*Ec(i+1)...
              -(dx/dy)*eps_bot/eps_si*Ec(i+Nx);
  end
end

%*****************Specify the F matrix in the silicon film region*****************
for i=Nx*(t_topa+1)+1:(Ntotal-Nx*(t_bota+1)+1)-1
  F(i)=-(dx/dy)*Ec(i-Nx)-(dy/dx)*Ec(i-1)+2*(dx/dy+dy/dx)*Ec(i)...
            +dummy_fun(i)...
            -(dy/dx)*Ec(i+1)-(dx/dy)*Ec(i+Nx);
end

%***************Modify the F matrix at the right and left boundaries**************
i_l=1;
i_r=Nx;
for j=1:Ny
  if(j==1)
      
    F(i_l)=2*Ec(i_l)-Ec(i_l+1)-Ec(i_l+Nx);
    F(i_r)=2*Ec(i_r)-Ec(i_r-1)-Ec(i_r+Nx);
    
  elseif(j>1 & j<round(Nx*(t_topa+1)/Nx))
      
    F(i_l)=Ec(i_l)-Ec(i_l+1);
    F(i_r)=Ec(i_r)-Ec(i_r-1);
    
  elseif(j>=round(Nx*(t_topa+1)/Nx) & j<=round(Ntotal-Nx*t_bota/Nx))
        
    F(i_l)=Ec(i_l)-Ec(i_l+1);
    F(i_r)=Ec(i_r)-Ec(i_r-1);
       
  elseif(j>round(Ntotal-Nx*t_bota/Nx) & j<Ny)
      
    F(i_l)=Ec(i_l)-Ec(i_l+1);
    F(i_r)=Ec(i_r)-Ec(i_r-1);
    
  elseif(j==Ny & ((Ntotal-Nx+1)<(Ntotal-Nx+1+Lsda)) &...
          ((Ntotal-Nx+1+Lsda+Lg_bota)<Ntotal))
      
    F(i_l)=2*Ec(i_l)-Ec(i_l+1)-Ec(i_l-Nx);
    F(i_r)=2*Ec(i_r)-Ec(i_r-1)-Ec(i_r-Nx);
    
  end
  i_l=1+j*Nx;
  i_r=(j+1)*Nx;
end

%*************************END OF EVALUATING F***********************************

%****************************EVALUATE MF_prime**********************************
%MF_prime matrix in the silicon film region
for i=Nx+1:(Ntotal-Nx+1)-1
   MF_prime(i,i)=dummy_fun_prime(i);
end

for j=2:(Ny-1) 
  MF_prime((j-1)*Nx+1,(j-1)*Nx+1)=0;
  MF_prime(j*Nx,j*Nx)=0;
end
%**************END OF EVALUATION FOR PENETRATION INTO OXIDE *******************
end
%*************END OF ELECTRON PENETRATION INTO OXIDE OPTION SWITCH*************

MF_prime=F_prime+sparse(MF_prime);

%**********************END OF EVALUATING MF_prime******************************

%***************************SOLVING FOR delta_Ec*******************************
delta_Ec=-sparse(MF_prime)\sparse(F);

for i=1:Ntotal
  if(abs(delta_Ec(i))<=1)
     delta_Ec(i)=delta_Ec(i);
  elseif(1<abs(delta_Ec(i)) & abs(delta_Ec(i)) <3.7)
%    delta_Ec(i)=sign(delta_Ec(i))*power(abs(delta_Ec(i)),1/5);
% octave does not have power function
     delta_Ec(i)=sign(delta_Ec(i))*(abs(delta_Ec(i)).^1/5);
  elseif(abs(delta_Ec(i))>=3.7)
     delta_Ec(i)=sign(delta_Ec(i))*log(abs(delta_Ec(i)));
  end
end
%*************************END OF SOLVING FOR delta_Ec*************************

Ec=Ec+delta_Ec;
error_inner=max(abs(real(F)));
fprintf ('%s %e \n','error_inner = ',error_inner);
max_delta_Ec=max(abs(real(full(delta_Ec))));
MF_prime=zeros(Ntotal,Ntotal);
F=zeros(Ntotal,1);
end 
%*************************END OF INNER LOOP (WHILE) **************************

%*****************************************************************************
%**************************END OF FUNCTION POISSON****************************
%*****************************************************************************
