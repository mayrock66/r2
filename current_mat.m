%****************************************************************************************
%***************CALCULATE THE CURRENT MATRIX Ii(muj) GIVEN TRANSMISSIONS***************** 
%***********************and Fermi energies at each scatters******************************
%*****************************(Zhibin Ren 10-12-00)**************************************
%****************************************************************************************

function [Is,Id,Iin,mu_scatter_new]=current_mat(mu_scatter_old,T_E,E);

%**************************FUNDAMENTAL physical constants********************************
global eps_si eps_o psi_si q k_B h_bar m_e m_t m_l ...
       Nc Ncc mx my mz Temp;
global nu_scatter Nx fermi_flag;

%**********************************INITIALIZATION****************************************
dmu=1e-6;
I_criterion=1e-2;%1E-8A/um

Iin=zeros(nu_scatter,1);
mu_scatter_new=zeros(nu_scatter+2,1);

delta_mu=zeros(nu_scatter+2,1);
I_tem=zeros(2);
IMU=zeros(nu_scatter,nu_scatter);
IMU_dummy=zeros(nu_scatter,nu_scatter);
T_dummy=sum(T_E,3);

%***************************************************************************************
%*************************************EVALUATE Iin**************************************
%***************************************************************************************
 for i_node=1:nu_scatter
  I_dummy1=...
  fermi(((mu_scatter_old(i_node)-E)/(k_B*Temp/q)),fermi_flag,-1/2)*...
  T_dummy(:,i_node);
  I_dummy2=0;
  for j_node=1:nu_scatter+2
    I_dummy2=I_dummy2+...
    fermi(((mu_scatter_old(j_node)-E)/(k_B*Temp/q)),fermi_flag,-1/2)*...
    T_E(:,i_node,j_node);
  end
  Iin(i_node)=I_dummy1-I_dummy2;
 end
 Isc=max(abs(Iin));
 
%***************************************************************************************
%*************************************EVALUATE IMU**************************************
%***************************************************************************************

if Isc>=I_criterion

 for i_node=1:nu_scatter
  IMU_dummy(i_node,i_node)=...
    ((fermi(((mu_scatter_old(i_node)+dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2)...
    -fermi(((mu_scatter_old(i_node)-dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2))/(2*dmu))*...
    T_dummy(:,i_node);
  for j_node=1:nu_scatter
    IMU(i_node,j_node)=...
    ((fermi(((mu_scatter_old(j_node)+dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2)...
    -fermi(((mu_scatter_old(j_node)-dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2))/(2*dmu))*...
    T_E(:,i_node,j_node);
  end
 end
 IMU=IMU_dummy-IMU;

end 

%***************************************************************************************
%**********************************END OF EVALUATE IMU**********************************
%***************************************************************************************

mu_scatter_new=mu_scatter_old;

%*********************************Newton searching loop*********************************

iiii=0;
while(Isc>=I_criterion)

 disp('Entering Jacobian loop in Current_mat');
 delta_mu(1:nu_scatter)=-sparse(IMU)\sparse(Iin);
 
 % The following IF statement performs a check to insure that the correction does not
 % bring the fermi level in the device above the source fermi level or below the drain fermi level.
 % If it does, then we averaged the fermi level of the scatterer to make it fit within that physical range of fermi levels.
 % If we don't perform that check, in some cases we can get a scatterer fermi level well below the drain fermi level
 % which forces charges to flow from the drain to fill the available state.......as a result -> NO convergence.
 
 mu_scatter_new=mu_scatter_new+delta_mu;
 
 for i=1:nu_scatter
 
 if mu_scatter_new(i)>mu_scatter_new(nu_scatter+1)
 	mu_scatter_new(i)=mu_scatter_new(nu_scatter+1);
 elseif mu_scatter_new(i)<mu_scatter_new(nu_scatter+2)
 	mu_scatter_new(i)=mu_scatter_new(nu_scatter+2);
 else
 	mu_scatter_new(i)=mu_scatter_new(i);
 end
 
 %if (max(mu_scatter_new(1:nu_scatter))>mu_scatter_new(nu_scatter+1) ...
 %| min(mu_scatter_new(1:nu_scatter))<mu_scatter_new(nu_scatter+2))
 
 %   mu_scatter_new(1:nu_scatter)=(mu_scatter_old(1:nu_scatter)+mu_scatter_new(1:nu_scatter))/2.0;
 %   fprintf(1,'%s%s%e\n','AVERAGED ','MAX CHANGE ',max(abs(mu_scatter_new(1:nu_scatter)-mu_scatter_old(1:nu_scatter)))); 
 
 %else
 %    fprintf(1,'%s%s%e\n','NOT AVERAGED ','MAX CHANGE ',max(abs(mu_scatter_new(1:nu_scatter)-mu_scatter_old(1:nu_scatter))));
 
 %end
 
 end %for i nu_scatter
 
%  for i_node=1:nu_scatter
%    if(abs(delta_mu(i_node))<=1)
%       delta_mu(i_node)=delta_mu(i_node);
%    elseif(1<abs(delta_mu(i_node)) & abs(delta_mu(i_node)) <3.7)
%       delta_mu(i_node)=sign(delta_mu(i_node))*power(abs(delta_mu(i_node)),1/5);
%    elseif(abs(delta_mu(i_node))>=3.7)
%       delta_mu(i_node)=sign(delta_mu(i_node))*log(abs(delta_mu(i_node)));
%    end
%  end

%  mu_scatter_new=mu_scatter_new+delta_mu;
  
  %****************************************************************************************
  %************************************ EVALUATE Iin***************************************
  %****************************************************************************************
  for i_node=1:nu_scatter
   I_dummy1=... 
   fermi(((mu_scatter_new(i_node)-E)/(k_B*Temp/q)),fermi_flag,-1/2)*...
   T_dummy(:,i_node);
   I_dummy2=0;
   for j_node=1:nu_scatter+2
     I_dummy2=I_dummy2+...
     fermi(((mu_scatter_new(j_node)-E)/(k_B*Temp/q)),fermi_flag,-1/2)*...
     T_E(:,i_node,j_node);
   end
   Iin(i_node)=I_dummy1-I_dummy2;
  end 
  Isc=max(abs(Iin))
  %****************************************************************************************
  %***********************************END OF EVALUATE Iin**********************************
  %****************************************************************************************

  if Isc>=I_criterion
  
  %****************************************************************************************
  %**************************************EVALUATE IMU**************************************
  %****************************************************************************************  
   for i_node=1:nu_scatter
    IMU_dummy(i_node,i_node)=...
      ((fermi(((mu_scatter_new(i_node)+dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2)...
      -fermi(((mu_scatter_new(i_node)-dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2))/(2*dmu))*...
      T_dummy(:,i_node);
    for j_node=1:nu_scatter
      IMU(i_node,j_node)=...
      ((fermi(((mu_scatter_new(j_node)+dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2)...
      -fermi(((mu_scatter_new(j_node)-dmu-E)/(k_B*Temp/q)),fermi_flag,-1/2))/(2*dmu))*...
      T_E(:,i_node,j_node);
    end
   end
   IMU=IMU_dummy-IMU;
  %****************************************************************************************
  %***********************************END OF EVALUATE IMU**********************************
  %****************************************************************************************
  end
  iiii=iiii+1
  
  % Copy old vals
  mu_scatter_old=mu_scatter_new;
  
end

for i_node=1:2
 I_dummy1=...
 fermi(((mu_scatter_new(i_node+nu_scatter)-E)/(k_B*Temp/q)),fermi_flag,-1/2)*...
 T_dummy(:,i_node+nu_scatter);
 I_dummy2=0;
 for j_node=1:nu_scatter+2
   I_dummy2=I_dummy2+...
   fermi(((mu_scatter_new(j_node)-E)/(k_B*Temp/q)),fermi_flag,-1/2)*...
   T_E(:,i_node+nu_scatter,j_node);
 end
 I_tem(i_node)=I_dummy1-I_dummy2;
end

Is=I_tem(1);
Id=I_tem(2);

%******************************************************************************************
%*************************** THE END OF FUNCTION CURRENT_MAT*******************************
%******************************************************************************************
