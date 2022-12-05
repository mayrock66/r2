%**************************************************************************
%*A M file function to read the input file using parser.m *****************
%**************************************************************************

function readinput

global transport_model;
global DG_flag fermi_flag ox_pnt_flag dummy_flag;
global t_vall max_subband;
global dopslope_s dopslope_d overlap_s overlap_d;
global Lsd Lg_top Lg_bot t_top t_bot t_si dx dy refine;
global Vg1 Vg2 Vs Vd Vd_initial Vg_step Ng_step Vd_step Nd_step;
global phi_top phi_bot eps_top eps_bot bar_top bar_bot;
global eps_si m_t m_l Te;
global mu_low beta Vel_sat;
global N_sd N_body criterion_outer criterion_inner;
global filename dirname;

%Energy transport parameters

global ELE_TAUW ELE_CQ;

%options for plotting

global plot_IV plot_Ec3d plot_Ne3d plot_Ecsub plot_Te plot_Nesub;
global plot_Ec_IV plot_Ne_IV;

filename=input('Enter filename to be run: ','s');
dirname=input('Enter name of output directory: ','s');

fin=fopen(filename,'rt');

p.err=1;
while p.err ~= -1
    p=parser(fin);

    if p.err == 1 | p.err == 999

       if strcmpi(p.ncard,'device')==1
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'nsd')==1 
             if strcmpi(p.var(i).type,'number')==1
               N_sd=1e6*p.var(i).val;
             else
               dump='Invalid assignment for nsd'
               return
             end
           elseif strcmpi(p.var(i).name,'nbody')==1 
             if strcmpi(p.var(i).type,'number')==1
               N_body=1e6*p.var(i).val;
             else
               dump='Invalid assignment for nbody'
               return
             end
           elseif strcmpi(p.var(i).name,'lgtop')==1 
             if strcmpi(p.var(i).type,'number')==1
               Lg_top=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for lgtop'
               return
             end
           elseif strcmpi(p.var(i).name,'lgbot')==1 
             if strcmpi(p.var(i).type,'number')==1
               Lg_bot=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for lgbot'
               return
             end
           elseif strcmpi(p.var(i).name,'lsd')==1 
             if strcmpi(p.var(i).type,'number')==1
               Lsd=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for lsd'
               return
             end
           elseif strcmpi(p.var(i).name,'overlap_s')==1 
             if strcmpi(p.var(i).type,'number')==1
               overlap_s=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for overlap_s'
               return
             end
            elseif strcmpi(p.var(i).name,'overlap_d')==1 
             if strcmpi(p.var(i).type,'number')==1
               overlap_d=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for overlap_d'
               return
             end 
           elseif strcmpi(p.var(i).name,'dopslope_s')==1 
             if strcmpi(p.var(i).type,'number')==1
               dopslope_s=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for dopslope_s'
               return
             end 
            elseif strcmpi(p.var(i).name,'dopslope_d')==1 
             if strcmpi(p.var(i).type,'number')==1
               dopslope_d=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for dopslope_d'
               return
             end  
           elseif strcmpi(p.var(i).name,'tsi')==1 
             if strcmpi(p.var(i).type,'number')==1
               t_si=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for tsi'
               return
             end
           elseif strcmpi(p.var(i).name,'tox_top')==1 
             if strcmpi(p.var(i).type,'number')==1
               t_top=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for tox_top'
               return
             end
           elseif strcmpi(p.var(i).name,'tox_bot')==1 
             if strcmpi(p.var(i).type,'number')==1
               t_bot=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for tox_bot'
               return
             end
           elseif strcmpi(p.var(i).name,'temp')==1 
             if strcmpi(p.var(i).type,'number')==1
               Te=p.var(i).val;
             else
               dump='Invalid assignment for temp'
               return
             end
           end
         end

       elseif strcmpi(p.ncard,'grid')==1 
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'dx')==1 
             if strcmpi(p.var(i).type,'number')==1
               dx=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for dx'
               return
             end
           elseif strcmpi(p.var(i).name,'dy')==1 
             if strcmpi(p.var(i).type,'number')==1
               dy=1e-9*p.var(i).val;
             else
               dump='Invalid assignment for dy'
               return
             end
           elseif strcmpi(p.var(i).name,'refine')==1 
             if strcmpi(p.var(i).type,'number')==1
               refine=p.var(i).val;
             else
               dump='Invalid assignment for refine'
               return
             end
           end
         end

       elseif strcmpi(p.ncard,'transport')==1 
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'model')==1 
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'dd')==1
                 transport_model=2;
               elseif strcmpi(p.var(i).val{1},'clbte')
                 transport_model=3;
               elseif strcmpi(p.var(i).val{1},'qbte')
                 transport_model=4;
               elseif strcmpi(p.var(i).val{1},'qdte') 
                 transport_model=5;  
               elseif strcmpi(p.var(i).val{1},'et') 
                 transport_model=6;  
               else
                 disp('****** ERROR !!! MODEL CAN ONLY BE DD/CLBTE/QBTE/ET/QDTE *******')
                 disp(p.var(i).name)
                 disp(p.var(i).type)
                 disp(p.var(i).val{1})
                 quit
               end
             else
               dump='Invalid assignment for model'
               return
             end
           elseif strcmpi(p.var(i).name,'mu_low')==1 
             if strcmpi(p.var(i).type,'number')==1
               mu_low=1e-4*p.var(i).val;
             else
               dump='Invalid assignment for mu_low'
               return
             end
           elseif strcmpi(p.var(i).name,'beta')==1 
             if strcmpi(p.var(i).type,'number')==1
               beta=p.var(i).val;
             else
               dump='Invalid assignment for beta'
               return
             end
           elseif strcmpi(p.var(i).name,'vsat')==1 
             if strcmpi(p.var(i).type,'number')==1
               Vel_sat=1e-2*p.var(i).val;
             else
               dump='Invalid assignment for vsat'
               return
             end
            elseif strcmpi(p.var(i).name,'ELE_TAUW')==1 
             if strcmpi(p.var(i).type,'number')==1
               ELE_TAUW=p.var(i).val;
             else
               dump='Invalid assignment for ELE_TAUW'
               return
             end 
             elseif strcmpi(p.var(i).name,'ELE_CQ')==1 
             if strcmpi(p.var(i).type,'number')==1
               ELE_CQ=p.var(i).val;
             else
               dump='Invalid assignment for ELE_CQ'
               return
             end
             
           end
         end

       elseif strcmpi(p.ncard,'bias')==1 
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'vgtop')==1
             if strcmpi(p.var(i).type,'number')==1
               Vg1=p.var(i).val;
             else
               dump='Invalid assignment for vgtop'
               return
             end
           elseif strcmpi(p.var(i).name,'vgbot')==1
             if strcmpi(p.var(i).type,'number')==1
               Vg2=p.var(i).val;
             else
               dump='Invalid assignment for vgbot'
               return
             end
           elseif strcmpi(p.var(i).name,'vs')==1
             if strcmpi(p.var(i).type,'number')==1
               Vs=p.var(i).val;
             else
               dump='Invalid assignment for vs'
               return
             end
           elseif strcmpi(p.var(i).name,'vd')==1
             if strcmpi(p.var(i).type,'number')==1
               Vd=p.var(i).val;
             else
               dump='Invalid assignment for vd'
               return
             end
           elseif strcmpi(p.var(i).name,'vgstep')==1
             if strcmpi(p.var(i).type,'number')==1
               Vg_step=p.var(i).val;
             else
               dump='Invalid assignment for vgstep'
               return
             end
           elseif strcmpi(p.var(i).name,'vdstep')==1
             if strcmpi(p.var(i).type,'number')==1
               Vd_step=p.var(i).val;
             else
               dump='Invalid assignment for vdstep'
               return
             end
           elseif strcmpi(p.var(i).name,'ngstep')==1
             if strcmpi(p.var(i).type,'number')==1
               Ng_step=p.var(i).val;
             else
               dump='Invalid assignment for ngstep'
               return
             end
           elseif strcmpi(p.var(i).name,'ndstep')==1
             if strcmpi(p.var(i).type,'number')==1
               Nd_step=p.var(i).val;
             else
               dump='Invalid assignment for ndstep'
               return
             end
           elseif strcmpi(p.var(i).name,'vd_initial')==1
             if strcmpi(p.var(i).type,'number')==1
               Vd_initial=p.var(i).val;
             else
               dump='Invalid assignment for vd_initial'
               return
             end
           end
         end

       elseif strcmpi(p.ncard,'material')==1 
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'wfunc_top')==1
             if strcmpi(p.var(i).type,'number')==1
              phi_top =p.var(i).val;
             else
               dump='Invalid assignment for wfunc_top'
               return
             end
           elseif strcmpi(p.var(i).name,'wfunc_bot')==1
             if strcmpi(p.var(i).type,'number')==1
               phi_bot=p.var(i).val;
             else
               dump='Invalid assignment for wfunc_bot'
               return
             end
           elseif strcmpi(p.var(i).name,'mlong')==1
             if strcmpi(p.var(i).type,'number')==1
               m_l=p.var(i).val;
             else
               dump='Invalid assignment for mlong'
               return
             end
           elseif strcmpi(p.var(i).name,'mtran')==1
             if strcmpi(p.var(i).type,'number')==1
               m_t=p.var(i).val;
             else
               dump='Invalid assignment for mtran'
               return
             end
           elseif strcmpi(p.var(i).name,'kox_top')==1
             if strcmpi(p.var(i).type,'number')==1
               eps_top=p.var(i).val;
             else
               dump='Invalid assignment for kox_top'
               return
             end
           elseif strcmpi(p.var(i).name,'kox_bot')==1
             if strcmpi(p.var(i).type,'number')==1
               eps_bot=p.var(i).val;
             else
               dump='Invalid assignment for kox_bot'
               return
             end
           elseif strcmpi(p.var(i).name,'dec_top')==1
             if strcmpi(p.var(i).type,'number')==1
               bar_top=p.var(i).val;
             else
               dump='Invalid assignment for dec_top'
               return
             end
           elseif strcmpi(p.var(i).name,'dec_bot')==1
             if strcmpi(p.var(i).type,'number')==1
               bar_bot=p.var(i).val;
             else
               dump='Invalid assignment for dec_bot'
               return
             end
           elseif strcmpi(p.var(i).name,'ksi')==1
             if strcmpi(p.var(i).type,'number')==1
               eps_si=p.var(i).val;
             else
               dump='Invalid assignment for ksi'
               return
             end
           end
         end

       elseif strcmpi(p.ncard,'solve')==1 
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'dvmax')==1
             if strcmpi(p.var(i).type,'number')==1
               criterion_outer=p.var(i).val;
             else
               dump='Invalid assignment for dvmax'
               return
             end
           end
           if strcmpi(p.var(i).name,'dvpois')==1
             if strcmpi(p.var(i).type,'number')==1
               criterion_inner=p.var(i).val;
             else
               dump='Invalid assignment for dvmax'
               return
             end
           end
         end
%*****************************************************************************         
        elseif strcmpi(p.ncard,'plots')==1 
         for i=1:p.nvar
             %Iv plot
           if strcmpi(p.var(i).name,'I_V')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_IV=1;
                else
                    plot_IV=0;
                end
             else
               dump='Invalid assignment for I_V'
               return
             end
           end
           %Ec plot
           if strcmpi(p.var(i).name,'Ec3d')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Ec3d=1;
                else
                    plot_Ec3d=0;
                end
             else
               dump='Invalid assignment for Ec'
               return
             end
           end
           %Ne plot
           if strcmpi(p.var(i).name,'Ne3d')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Ne3d=1;
                else
                    plot_Ne3d=0;
                end
             else
               dump='Invalid assignment for Ne'
               return
             end
           end
           %Ec_sub plot
           if strcmpi(p.var(i).name,'Ec_sub')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Ecsub=1;
                else
                    plot_Ecsub=0;
                end
             else
               dump='Invalid assignment for Ecsub'
               return
             end
           end
           %Ne_sub plot
           if strcmpi(p.var(i).name,'Ne_sub')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Nesub=1;
                else
                    plot_Nesub=0;
                end
             else
               dump='Invalid assignment for Ne_sub'
               return
             end
           end
           %Te plot
           if strcmpi(p.var(i).name,'Te')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Te=1;
                else
                    plot_Te=0;
                end
             else
               dump='Invalid assignment for Te'
               return
             end
           end
           %Ec_IV plot
           if strcmpi(p.var(i).name,'Ec_IV')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Ec_IV=1;
                else
                    plot_Ec_IV=0;
                end
             else
               dump='Invalid assignment for Ec_IV'
               return
             end
           end
           %Ne_IV plot
           if strcmpi(p.var(i).name,'Ne_IV')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'y')==1
                    plot_Ne_IV=1;
                else
                    plot_Ne_IV=0;
                end
             else
               dump='Invalid assignment for Ne_IV'
               return
             end
           end
         end
%*****************************************************************************         
       elseif strcmpi(p.ncard,'options')==1 
         for i=1:p.nvar
           if strcmpi(p.var(i).name,'valleys')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'unprimed')==1
                 t_vall=1;
               else
                 t_vall=3;
               end
             else
               dump='Invalid assignment for valleys'
               return
             end
           elseif strcmpi(p.var(i).name,'num_subbands')==1
             if strcmpi(p.var(i).type,'number')==1
               max_subband=p.var(i).val;
             else
               dump='Invalid assignment for num_subbands'
               return
             end
           elseif strcmpi(p.var(i).name,'dg')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'true')==1
                 DG_flag=1;
               else
                 DG_flag=0;
               end
             else
               dump='Invalid assignment for dg'
               return
             end
           elseif strcmpi(p.var(i).name,'fermi')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'true')==1
                 fermi_flag=1;
               else
                 fermi_flag=0;
               end
             else
               dump='Invalid assignment for fermi'
               return
             end
           elseif strcmpi(p.var(i).name,'ox_penetrate')==1
             if strcmpi(p.var(i).type,'string')==1
               if strcmpi(p.var(i).val{1},'true')==1
                 ox_pnt_flag=1;
               else
                 ox_pnt_flag=0;
               end
             else
               dump='Invalid assignment for ox_penetrate'
               return
             end
           end
         end
       end
    end
end

p.err=1;

if transport_model~=3
  if (Vg1>1|Vg2>1|Vs>1|Vd>1)
    dump='Too High of Voltage'
    return
  end 
end

%if (Nd_step>0&Ng_step>0)
%    dump='Only one voltage sweep can be defined'
%    return
%end

if (Ng_step>20|Nd_step>20)
  dump='More than 20 bias points specified'
  return
end

%CHECKING ALL INPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DG_flag~=1 & DG_flag~=0
  ver='******ERROR, DG_flag can only be 0 or 1!!!******';
  disp(ver)
  quit
end

if fermi_flag~=1 & fermi_flag~=0
  ver='******ERROR, fermi_flag can only be 0 or 1!!!******';
  disp(ver)
  quit
end

if ox_pnt_flag~=1 & ox_pnt_flag~=0
  ver='******ERROR, ox_pnt_flag can only be 0 or 1!!!******';
  disp(ver)
  quit
end

if Lsd>15e-9
  Lsd=15e-9;
  ver='******NOTE! LSD is reduced to 15nm !!!******';
  disp(ver)
  pause(2)
end

if Lg_top>50e-9
  Lg_top=50e-9;
  ver='******NOTE! LGTOP is reduced to 50nm !!!******';
  disp(ver)
  pause(2)
end

if Lg_bot>1.2*Lg_top;
  Lg_bot=1.2*Lg_top;
  ver='******NOTE! LGBOT is reduced to 1.2*LGTOP !!!******';
  disp(ver)
  pause(2)
end

if t_si>5e-9;
  ver='******Please specify a TSI less than 5nm !!!******';
  disp(ver)
  pause(2)
end

if t_vall~=1 & t_vall~=3
  ver='******ERROR, VALLEYS can only be UNPRIMED or ALL!!!******';
  disp(ver)
  quit
end

if max_subband>3
  max_subband=3;
  ver='******Note! NUM-SUBBAND is limited to 3 for each VALLEY!!!******';
  disp(ver)
  pause(2)
end

%if t_si<=2e-9
%  t_vall=1;
%  max_subband=1;
%  ver='******Note! For TSI<2nm, only one subband is assumed!!!******';
%  disp(ver)
%  pause(2)
%end

if criterion_outer<=1e-4 & (transport_model==2|transport_model==3);
 criterion_outer=1e-4;
 ver='******Note! DVMAX for models CLBTE and QBTE is limited to 0.1meV!!!******';
 disp(ver)
 pause(2)
end

if t_si<2e-9 & Te>250 
 dummy_flag=0;
else
 dummy_flag=1/2;
end

if transport_model==6
    fermi_flag=0;
end
%********************************************************************************
%**************** THE END OF READINPUT FUNCTION *********************************
%********************************************************************************
