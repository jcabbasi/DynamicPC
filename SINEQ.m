%% Imbibition in presence of gravity forces considering non-equilibrium forces.
%% This Function uses Incompressible assumption to simulate fluid flow properties
% Incompressible problems: Pressure is not time dependent, Saturations are
% time dependent, 
% Start date of coding: 6/25/2016 - 5 Tir 1395
tic
%% This function optimizes grid dimensions for specifed matrix dimensions
clc
clear all
run_num_total=10            ;
scen_num=1                  ;
lower_side_p=1              ;
upper_side_p=1              ;
left_lateral_side_p=1       ;
right_lateral_side_p=1      ;
frontal_lateral_side_p=1    ;
backside_lateral_side_p=1   ;
fracture_water_saturation=1 ;
matrix_water_saturation=0.2 ;
system_pressure_upper=1000  ;
%%%%
matrix_length_x=.1           ;
matrix_length_y=.1           ;
%
matrix_length_z=1           ;


%%%%
fracture_length_top=5;
fracture_length_bottom=5;
%%%%


matrix_porosity=.2
matrix_permx=10
matrix_permy=10
matrix_permz=10
fracture_permx=100000
fracture_permy=100000
fracture_permz=100000
fracture_poro=1

% Fluid Properties:
density_oil=800; % Kg/m3
density_water=800; % Kg/m3

viscosityOil=10;  % cp
viscosityWater=1; % cp
%

conv_condition=0.001;



swof=[0 0 1 250  ; 1 1 0 0];
% Grid Numbers for first study:
% Matrix Grid Numbers:
matrix_grid_numbers=10;


% Adding MRST MODULES:
mrstModule add ad-core ad-blackoil ad-props mrst-gui incomp

%% Equalizing parameters:
% matrix_permy=matrix_permx;
% matrix_permz=matrix_permx;

%%
goptimizing_try=0;
error=1;

    
goptimizing_try=goptimizing_try+1;    
% Total Numbers of Runs:
run_num_total=run_num_total+1;
%  grid_dimensions_l=exp((1:matrix_grid_numbers));
% grid_dimensions_r=exp(1./(1:matrix_grid_numbers));
% grid_dimensions=grid_dimensions_r.*grid_dimensions_l
 
 grid_number=matrix_grid_numbers+2;
 grid_dimensions=matrix_length_z/matrix_grid_numbers;
 grid_dimensions=grid_dimensions*ones(matrix_grid_numbers,1);
 %
 %Grid Dimensions:
 grid_dimensions_x=matrix_length_x*ones(grid_number,1);
 grid_dimensions_y=matrix_length_y*ones(grid_number,1);
 grid_dimensions_z=[fracture_length_top; grid_dimensions; fracture_length_bottom];
 %
  
%% Preparing Input Values for MRST:
grid_num=matrix_grid_numbers+2;
% DX Values:
% Assuming this trend for cycling:
% x > y > z


perm_x_tot=[fracture_permx; matrix_permx*ones(matrix_grid_numbers,1);fracture_permx];
perm_y_tot=[fracture_permy; matrix_permy*ones(matrix_grid_numbers,1);fracture_permy];
perm_z_tot=[fracture_permz; matrix_permz*ones(matrix_grid_numbers,1);fracture_permz];


sw_init_tot=[fracture_water_saturation; matrix_water_saturation*ones(matrix_grid_numbers,1);fracture_water_saturation];
poro_tot=[fracture_poro; matrix_porosity*ones(matrix_grid_numbers,1);fracture_poro];




%% Pressure Calculator
depths=[fracture_length_top; (fracture_length_top+((1:matrix_grid_numbers)/matrix_grid_numbers)*matrix_length_z)'; (fracture_length_top+matrix_length_z+fracture_length_bottom)];

Pressure=system_pressure_upper*ones(grid_number,1)+density_oil*(1/144)*depths;



% Preparing Static Model
system_z_dims_cum=0;
  
for ii=1:grid_num
system_z_dims_cum(1,ii+1)=system_z_dims_cum(1,ii)+grid_dimensions_z(ii);
end
    

%% Preparing MRST Simulator:

% Building Model Structure
G = tensorGrid([0 matrix_length_x]*ft, [0 matrix_length_y]*ft, system_z_dims_cum*ft);
G = computeGeometry(G);
% Building Static PRoperties
rock = makeRock(G,[perm_x_tot perm_y_tot perm_z_tot]*milli*darcy, poro_tot);
% 
% 
% clf;
% plotGrid(G)
% view(30, 40)
% 
% 
% clf;
% plotCellData(G, rock.perm(:,1))
% colorbar


% Fluid Properties:
% For incompressible cases:
props = constantProperties([   viscosityWater,  viscosityOil] .* centi*poise, ...
                           [density_water, density_oil] .* kilogram/meter^3);

[kr, pc] = tabulatedSatFunc([swof(:,1),swof(:,2),swof(:,3),swof(:,4)*psia]);                      
                       
% Then make another fluid object identical to the one above except for the
% capillary pressure term 'pc'.
                 
                       
fluid = struct('properties', props                  , ...
                  'saturation', @(x, varargin)    x.s  , ...
                  'relperm'   , kr                     , ...
                  'pc'        , @(x, varargin) pc(x.s));

% Enabling Gravity:
gravity reset on


% Initialization
state = initState(G,[], Pressure*psia, [sw_init_tot (1-sw_init_tot)]);



%  for ii=1:12
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
% plotCellData(G,state.pressure(:,1) , j == ii)
% % Plot the wells
% view(30,50)
% colorbar
% pause(0.5)
% end



% Computing Transmissibilities:

% Verbose: Shows information about the convergence conditions
verbose = true;


trans  = computeTrans(G, rock, 'Verbose', verbose);

% Initial Conditions:
stateNEQ=state;
%



%% Solving Non Equilibrium Parameters:



% Scheduling:
% Preparing timesteps:

timesteps=[
    %0.00001
% 1.29953E-05 
% 1.68877E-05 
% 2.1946E-05 
% 2.85194E-05 
% 3.70617E-05
% 4.81627E-05
% 6.25887E-05
% 8.13357E-05
% 0.000105698
% 0.000137357
% 0.000178499
% 0.000231965
% 0.000301444
% 0.000391735
% 0.00050907
% 0.00066155
% 0.000859701
% 0.001117205
% 0.001451837
% 0.001886701
% 0.002451818
% 0.003186203
% 0.004140555
% 0.005380761
% 0.006992442
% 0.009086864
% 0.01180862
0.015345616
0.019942035
0.025915204
0.033677495
0.043764799
0.056873518
0.073908646
0.096046247
0.124814647
0.162199947
0.210783136
0.27391828
0.355964075
 0.462584764
 0.601141179
 0.781198918]
% 1.015188729
% 1.319264699
% 1.714419493
% 2.227933635
% 2.895258894
% 3.76246578
% 4.889424146
% 6.353936455
% 8.257109072
% 10.7303324
% 13.94435177
% 18.12105524
% 23.54879226
% 30.60228058
% 39.76847588
% 51.68018995]
% 67.15977856
% 87.27591482
% 113.4173678
% 147.3888799
% 191.5357614
% 248.9058058
% 323.4597012]
% 420.3444672
% 546.248792
% 709.864804
% 922.4881544
% 1198.797842
% 1557.869614
% 2024.492913
% 2630.882276
% 3418.901348]

% timesteps=0.01*ones(200,1)
%% Incompressibile Solver:

psolve  = @(state, stateNEQ, fluid) incompTPFANEQ(state,stateNEQ, G, trans, fluid);
tsolve  = @(state, stateNEQ, dT, fluid) implicitTransportNEQ(state, stateNEQ, G, dT, rock, ...
                                                fluid, ...
                                                'verbose', verbose);

currenttime=0;
conv_condition=0.1e-27;
iter_total=0;
IterMax=15;
% Non Equilibrium Parameters:
tau=000;
cTau=0e8;

rSol=state; % In Order to get prepared for the loops
rSol_NEQ=rSol;
for ii=1:size(timesteps)
 
 % Current Simulation Time:
currenttime(ii+1)=currenttime(ii)+timesteps(ii);
 % Sets The Known Parameters to be constant during NEQ-iterations:
rSol_tsinit=rSol;
iter_local(ii)=0;
% In the first NEQ-Iteration EQ Conditions must be defined
rSol_NEQ=rSol;
% Into Entering While loop
conv_criteria=conv_condition+eps;
while abs(conv_criteria)>conv_condition
 
      iter_total=iter_total+1;    
      iter_local(ii)=iter_local(ii)+1;    
 
      % Solving Pressure Equations:
      rSol = psolve(rSol_tsinit,rSol_NEQ, fluid);
      
      % If it is the first iteration, Saves initial Saturations:
      if iter_total==1
      Sw_mat_dyna(:,1) = rSol.s(:,1)    
      end
      % TRANSPORT SOLVE
      rSol = tsolve(rSol,rSol_NEQ, timesteps(ii)*day, fluid);

 if iter_local(ii)==1
    rSolEQ{ii}=rSol;
 end
 
  % Check for inconsistent saturations
   s = [rSol.s(:,1)];                       %
   assert(max(s) < 1+eps && min(s) > -eps); % Throws and error when condition is false  
      
      % Saving Matrix saturations for each total_iterations. In order to
      % calculate _etta.
      Sw_mat_dyna_init=rSol_tsinit.s(:,1);
      Sw_mat_dyna(:,iter_total+1)=rSol.s(:,1);
 
   
      pcTau(:,iter_total)=fluid.pc(rSol);
      
      tauDynamic=(cTau*viscosityWater./pcTau(:,1));
      % by a way calculate Ds/Dt 
      % Ds/dt t:second
      sat_change_vs_time(:,iter_total)=((Sw_mat_dyna(:,iter_total+1)-Sw_mat_dyna_init)/(timesteps(ii)*86400));

%       if iter_local(ii)>2
% sat_change_vs_time(:,iter_total)=(sat_change_vs_time(:,iter_total)+sat_change_vs_time(:,iter_total-1))/2;
%       
%       end          
          
 %% Updating Non_Equilibrium Effects in into modify effective water saturations:
 
     % Effective Saturations:
     Sw_etta_NEQ(:,iter_total)=rSol_tsinit.s(:,1)+tauDynamic.*sat_change_vs_time(:,end); % Etta!
     Sw_etta_NEQCheck(:,1)=Sw_etta_NEQ(:,iter_total);
       Sw_etta_NEQCheckSPG(iter_total,1)=Sw_etta_NEQ(2,iter_total);

       
       % Checking if the new saturations are in the correct range of 0 to 1:
     Sw_etta_NEQ(find(Sw_etta_NEQ(:,iter_total)>1),iter_total)=1;
     Sw_etta_NEQ(find(Sw_etta_NEQ(:,iter_total)<0),iter_total)=0;
     
     
     if iter_local(ii)>1
          Sw_etta_NEQ(:,iter_total)=(Sw_etta_NEQ(:,iter_total)+Sw_etta_NEQ(:,iter_total-1))/2;
      end
     
      Sw_etta_NEQCheck(:,2)=Sw_etta_NEQ(:,iter_total);
      Sw_etta_NEQCheckSPG(iter_total,2)=Sw_etta_NEQ(2,iter_total);


     % Preparing
     rSol_NEQ=rSol;
     rSol_NEQ.s=[Sw_etta_NEQ(:,iter_total) , (1-Sw_etta_NEQ(:,iter_total))];

     pc_eff(:,iter_total)=fluid.pc(rSol_NEQ);
     iter_local
     
   % Check for inconsistent saturations
   sneq = [rSol_NEQ.s(:,1)];                       %
   assert(max(sneq) < 1+eps && min(sneq) > -eps); % Throws and error when condition is false    
     
   %% Checking Convergence:
   % Considering Saturation Changes Just in the Matrix Flow !!!!!!
   if tau==0 | cTau==0
       
   conv_criteria=0;

   else
       
       if iter_local(ii)==1    
       conv_criteria=conv_condition+eps;
       else
       d2sw_dt2_mean(iter_total)=mean( (sat_change_vs_time(2:end-1,iter_total)-sat_change_vs_time(2:end-1,iter_total-1))/timesteps(ii)  );  
       SwChangeNEQIter(iter_total)=mean(Sw_mat_dyna(:,iter_total+1)-Sw_mat_dyna(:,iter_total));
       conv_criteria=SwChangeNEQIter(iter_total)   
       end   
         
       
   end   
    
   %% In order to check the convergence rate of the simulation 
   
       conv_criteriaStorage(iter_total)=conv_criteria;
       
       if iter_local(ii)==1
           
       else 
           ConvRate(iter_total)=((abs(conv_criteriaStorage(iter_total))-abs(conv_criteriaStorage(iter_total-1)))/conv_criteriaStorage(iter_total));
      
         if   abs(ConvRate(iter_total))<0.000001
            conv_criteria=0
            
         end
         
       if iter_local(ii)>2
             
             ConvRate2(iter_total)=((abs(ConvRate(iter_total))-abs(ConvRate(iter_total-1)))/ConvRate(iter_total));
        
            if   abs(ConvRate2(iter_total))<0.000001
            conv_criteria=0
            
            end
       end
             
       end
       
       if iter_local(ii)<IterMax
        
       else
     conv_criteria=0      
           
       end
           
       
end
SwNEQEffect(:,ii)=rSol.s(:,1)-rSolEQ{ii}.s(:,1);
SwNEQEffectNormed(:,ii)=SwNEQEffect(:,ii)./rSolEQ{ii}.s(:,1);
states{ii}=rSol; 

 % Update solution of pressure equation.
end  



% Averaging:
for ii=1:size(timesteps)
SwMatrixResult(ii)=mean(states{ii}.s(2:end-1,1));
SoMatrixResult(ii)=mean(states{ii}.s(2:end-1,2));
end

resultsNEQsW=[currenttime(2:end)',SwMatrixResult'];
  
hold on  
plot(currenttime(2:end),SwMatrixResult)
% plot(currenttime(2:end),SoMatrixResult)
% plot(currenttime(2:end),iter_local(:))
hold off

% for ii=1:14 ; plot(1:12,states{1,ii}.s(:,2)); pause(0.5) ; end

% for ii=1:size(timesteps)
% plot(currenttime(ii),states{ii}.s(1,1),'Marker','*')
% 
% pause(0.4)
% end   
%% Simulate base case
% Once we have defined the schedule (dynamic controls and time), model
% (mathematical description of how to advance the solution) and the initial
% solution, we can simulate the problem.
% [wellreps , states] = simulateScheduleAD(initialize, model, schedule);

% for ii=1:size(timesteps)
%     
% saturation_3d(:,:,:,ii)=reshape(states{ii}.s(:,1),[grid_num grid_num grid_num]);
% 
% 
% end

% for ii=1:size(timesteps)
% saturation_water_avg(:,:,ii)=mean(saturation_3d(2:end-1,2:end-1,2:end-1,ii),3);
% end

% for ii=1:12
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% plotGrid(G, 'FaceAlpha', 0, 'EdgeAlpha', .1)
% plotCellData(G,rock.perm(:,1) , j == ii)
% view(30,50)
% colorbar
% pause(0.5)
% end


% 
% [i j k] = ind2sub(G.cartDims, 1:G.cells.num);
% clf;
% plotGrid(G, 'FaceAlpha', 0.1, 'EdgeAlpha', .1)
% plotCellData(G, states{1}.s(:,1), j == 5)
% view(30,50)

%SUMMARY OF RUN CSISIM1                                                                                                                                                                                                                                                                                                                                                                                                                                                                   		
%TIME        	YEARS       	ROE                                                                                                    
%DAYS        	YEARS       	                                                                                                       


%% Change properties from vector to matrix:


% 
% for tt=1:1:numel(timesteps)
% [X_index,Y_index] = meshgrid(1:grid_num,1:grid_num);
% 
% surf(X_index,Y_index,results_saturations.water(:,:,grid_num/2,tt));
% 
% pause(0.5);
% 
% end



toc
