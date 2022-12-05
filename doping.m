% A M file to specify the function that generates the doping profile
% With different overlaps and different doping gradient on the drain and source side
% By sebastien Goasguen October 2001/ Purdue

function doping

global dopslope_s dopslope_d overlap_s overlap_d;
global junction_l junction_r ox_pnt_flag;
global Nx Ny Ntotal;
global Lsd Lg_top Lg_bot t_top t_bot t_si dx dy refine Te;
global Nd N_sd N_body ;

Lsda=round(Lsd/dx);
Lg_topa=round(Lg_top/dx);
Lg_bota=round(Lg_bot/dx);
t_topa=round(t_top/dy);
t_bota=round(t_bot/dy);
t_sia=round(t_si/dy);


if (dopslope_s~=0 | dopslope_d~=0)
    decay_lens=(2/sqrt(2.3)*dopslope_s);
    decay_lend=(2/sqrt(2.3)*dopslope_d);
%***********************ASSUMING GAUSSIAN DISTRIBUTION************************
 for iii_row=((Nx*(t_topa+1))/Nx):((Ntotal-Nx*t_bota)/Nx-2)
  for iii_col=1:junction_l
    i_node=iii_row*Nx+iii_col;
    Nd(i_node)=N_sd-N_body;
  end
  
  for iii_col=junction_l+1:junction_r-1
    i_node=iii_row*Nx+iii_col;
    if (dopslope_s~=0 & dopslope_d~=0)
        Nd(i_node)=N_sd*exp(-((iii_col-junction_l)*dx/decay_lens)^2)+...
               N_sd*exp(-((iii_col-junction_r)*dx/decay_lend)^2)-...
               N_body;
    elseif (dopslope_s~=0 & dopslope_d==0)
         Nd(i_node)=N_sd*exp(-((iii_col-junction_l)*dx/decay_lens)^2)-...
               N_body;  
    elseif (dopslope_s==0 & dopslope_d~=0)
         Nd(i_node)=N_sd*exp(-((iii_col-junction_r)*dx/decay_lend)^2)-...
               N_body;         
    end
    
  end

  for iii_col=junction_r:Nx
    i_node=iii_row*Nx+iii_col;
    Nd(i_node)=N_sd-N_body;
  end
 end

 if ox_pnt_flag==1
     
  Nd((Nx*t_topa+1):(Nx*t_topa+Nx))=...
      Nd((Nx*t_topa+1+Nx):(Nx*t_topa+2*Nx));  
  Nd((Ntotal-Nx*(t_bota+1)+1):(Ntotal-Nx*(t_bota+1)+Nx))=...
      Nd((Ntotal-Nx*(t_bota+1)+1-Nx):(Ntotal-Nx*(t_bota+1)));  
 end
 
%*********************ABRUPT PROFILE on BOTH SIDE*************************
elseif (dopslope_s==0 & dopslope_d==0)

 for iii_row=((Nx*(t_topa+1))/Nx):((Ntotal-Nx*t_bota)/Nx-2)
  for iii_col=1:junction_l
    i_node=iii_row*Nx+iii_col;
    Nd(i_node)=N_sd-N_body;
  end
  
  for iii_col=junction_l+1:junction_r-1
    i_node=iii_row*Nx+iii_col;
    Nd(i_node)=-N_body;
  end

  for iii_col=junction_r:Nx
    i_node=iii_row*Nx+iii_col;
    Nd(i_node)=N_sd-N_body;
  end
 end

 if ox_pnt_flag==1
  Nd((Nx*t_topa+1):(Nx*t_topa+1)+junction_l-1)=(N_sd-N_body)/2;
  Nd((Nx*t_topa+1)+junction_l:(Nx*t_topa+1)+junction_r-2)=-N_body/2;
  Nd((Nx*t_topa+1)+junction_r-1:(Nx*(t_topa+1)))=(N_sd-N_body)/2;
  Nd((Ntotal-Nx*(t_bota+1)+1):(Ntotal-Nx*(t_bota+1)+1)+junction_l-1)=(N_sd-N_body)/2;
  Nd((Ntotal-Nx*(t_bota+1)+1)+junction_l:(Ntotal-Nx*(t_bota+1)+1)+junction_r-2)=-N_body/2;
  Nd((Ntotal-Nx*(t_bota+1)+1)+junction_r-1:(Ntotal-Nx*t_bota))=(N_sd-N_body)/2;
 end
end


%************************* THE END OF FUNCTION DOPING ******************************************************