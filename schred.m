%****************************************************************************
%******Solving the 1D Schroedinger's equation within vertical slices.********
%*************************(Zhibin Ren 7-28-00)*******************************
%****************************************************************************

function [E_v, W_v]=schred(Ec_old)

%**********************FUNDAMENTAL physical constants************************
global eps_si eps_o psi_si q k_B h_bar m_e m_t m_l ...
       Nc Ncc mx my mz Temp;
%***************************BOUNDARY AND BIAS********************************
global ox_pnt_flag t_vall max_subband;
global delta_x Lsd Lg_top Lg_bot t_top t_bot t_si dx dy refine;
global phi_top phi_bot eps_top eps_bot bar_top bar_bot;
global eps_si m_t m_l Temp;
global Nx Ny Ntotal;
global N_sd N_body; 

%	INPUT AND OUTPUT VARIABLES
%Ec_old is the old potential energy profile in eV
%1 column of Ntotal elements
%E_v: bands formed by subband energies in vertical direction 
%W_v: distribution function in vertical direction

%***************************TEMPORARY VARIABLES******************************

if ox_pnt_flag==0%(NO ELECTRON PENETRATION INTO OXIDE REGIONS)
   t_sch=t_si;
   Ec_start=Nx*t_top/(dy/refine)+1;
   Ec_end=Ntotal-Nx*t_bot/(dy/refine);
elseif ox_pnt_flag==1%(ACCOUNTING FOR ELECTRON PENETRATION INTO OXIDE REGIONS)
   t_sch=t_top+t_si+t_bot;
   Ec_start=1;
   Ec_end=Ntotal;
end

Np_old=round(t_sch/dy)+1;
x_dummy_old=(linspace(0,t_sch,Np_old))';
Np_new=round(t_sch/(dy/refine))+1;
x_dummy_new=(linspace(0,t_sch,Np_new))';

%******************************INITIALIZATION********************************

Ec_old=real(Ec_old);
E_v=zeros(Nx,max_subband,t_vall);
W_v=zeros(Np_old,Nx,max_subband,t_vall);
W_v_tem_1=zeros(Np_new,1);
W_v_tem_2=zeros(Np_old,1);
MEc=zeros(Np_old,Nx);%Potential in the silicon region
Ec_start=round(Ec_start);
Ec_end=round(Ec_end);
MEc=(reshape(Ec_old(Ec_start:Ec_end),Nx,Np_old))';

if ox_pnt_flag==0
   Ec_mod=zeros(Np_new,1);
elseif ox_pnt_flag==1
   Np_top=round(t_top/(dy/refine));
   Np_bot=round(t_bot/(dy/refine));
   Np_si=round(t_si/(dy/refine))+1;
   Ec_top=bar_top*ones(Np_top+1,1);
   Ec_bot=bar_bot*ones(Np_bot+1,1);
   Ec_si=0*ones(Np_si-2,1);
   Ec_mod=[Ec_top;Ec_si;Ec_bot];
end
%*****************************************************************************
%*******************************MAIN COMPUTATION******************************
%*****************************************************************************

for iii_vall=1:t_vall
 m_ee=mz(iii_vall)*m_e; 
 if iii_vall==3
   E_v(:,:,3)=E_v(:,:,2);
   W_v(:,:,:,3)=W_v(:,:,:,2);
   break;
 end
  
 tt=(h_bar^2)/(2*m_ee*((dy/refine)^2)*q);

for iii_col=1:Nx
    
    if refine==1
     U_vertical=MEc(:,iii_col);
    else
     U_vertical=interp1(x_dummy_old,MEc(:,iii_col),x_dummy_new,'spline');
    end

    U_vertical=U_vertical+Ec_mod;
    H=tt*((2*eye(Np_new-2))-(diag(ones(Np_new-1-2,1),1))...
       -(diag(ones(Np_new-1-2,1),-1)))+diag(U_vertical(2:Np_new-1));
    [evac,eval]=eig(H);
    [meval,i_order]=sort(diag(eval));
    E_v(iii_col,:,iii_vall)=(meval(1:max_subband))';
    for i_counter=1:max_subband
       W_v_tem_1(2:Np_new-1)=conj(evac(:,i_order(i_counter)))...
          .*evac(:,i_order(i_counter));
       if refine==1
        W_v_tem_2=W_v_tem_1;
       else
        W_v_tem_2=interp1(x_dummy_new,W_v_tem_1,x_dummy_old,'spline');
       end
       W_v(:,iii_col,i_counter,iii_vall)=W_v_tem_2./sum(W_v_tem_2);
    end
 end
end

%**************************************************************************
%************************END OF OF SCHRED**********************************
%**************************************************************************

