clear variables
close all

% CTRL stands for Control.
% FN stands for File Name.
FLEXIBLE_SATELLITE_SIMULINK_FN = 'flexible_satellite';
REACTION_WHEEL_SIMULINK_FN = 'reaction_wheel';
MODAL_CTRL_OPEN_LOOP_SIMULINK_FN = 'modal_control_open_loop_system';
MODAL_CTRL_CLOSED_LOOP_SIMULINK_FN = 'modal_control_closed_loop_system';
MODAL_CTRL_SISO_OPEN_LOOP_SIMULINK_FN = 'modal_control_siso_open_loop_system';
MODAL_CTRL_WITH_INTEGRATOR_OPEN_LOOP_SIMULINK_FN = 'modal_control_with_integrator_open_loop_system';
PHASE_LEAD_CTRL_SISO_OPEN_LOOP_SIMULINK_FN = 'phase_lead_control_open_loop_system';
PHASE_LEAD_CTRL_SISO_CLOSED_LOOP_SIMULINK_FN = 'phase_lead_control_closed_loop_system';
PHASE_LEAD_CTRL_SISO_CONTROLLER_OPEN_LOOP_SIMULINK_FN = 'phase_lead_control_controller_open_loop_system';
PHASE_LEAD_CTRL_WITH_PI_SISO_CONTROLLER_OPEN_LOOP_SIMULINK_FN = 'phase_lead_control_with_pi_controller_open_loop_system';
HINF_CTRL_OPEN_LOOP_SIMULINK_FN = 'hinf_control_open_loop_system';
HINF_CTRL_NO_DERIVATIVE_FEEDBACK_OPEN_LOOP_SIMULINK_FN = 'hinf_control_no_derivative_feedback_open_loop_system';
HINF_CTRL_WITH_INTEGRATOR_OPEN_LOOP_SIMULINK_FN = 'hinf_control_with_integrator_open_loop_system';
HINF_CTRL_WITH_INTEGRATOR_MULTIMODEL_OPEN_LOOP_SIMULINK_FN = 'hinf_control_with_integrator_multimodel_open_loop_system';

SATELLITE_INERTIA = 31.38;                          %Js
FLEXIBLE_MODE_FREQUENCY = 2.6268;                   %ω
FLEXIBLE_MODE_DAMPING_RATIO = 0.00495;              %ξ
COUPLING_COEFFICIENT = 0.5736;                      %α
REACTION_WHEEL_INERTIA = 0.0004;                    %Jr
COMMANDED_TORQUE_SATURATION = 0.05;                 %TcMax
REACTION_WHEEL_VELOCITY_SATURATION = 2800*2*pi/60;  %ωrMax
MEASUREMENT_DELAY = 0.8;                            %τm
MEASUREMENT_NOISE_VARIANCE = 10^(-9);               %σ²m
ESTIMATOR_TIME_CONSTANT = 0.5;                      %τe
DISTURBANCE_TORQUE = 0.001;                         %Td
INITIAL_CONDITION_ON_STATE = 0.0001;                %η0

DISTURBANCE_TORQUE_DELAY = 100;
STAR_TRACKER_NOISE_DELAY = 160;

SMALL_POINTING_ERROR = 0.5*2*pi/360;
CUSTOM_POINTING_ERROR = 2.7*2*pi/360;
LARGE_POINTING_ERROR = 10*2*pi/360;
MAXIMUM_TOLERABLE_POINTING_ERROR = 3.5*2*pi/360;

REACTION_WHEEL_DYNAMICS_NUMERATOR = REACTION_WHEEL_INERTIA*[1.214 0.763 0];	%linear dynamics of the wheel, numerator
REACTION_WHEEL_DYNAMICS_DENOMINATOR = [1 2.400 0.763];                     	%linear dynamics of the wheel, denominator

REACTION_WHEEL_FIRST_ORDER_TIME_CONSTANT = 1/0.38;
REACTION_WHEEL_FIRST_ORDER_DYNAMICS_NUMERATOR = [0 1];                                          %linear dynamics of the wheel, 1st order, numerator
REACTION_WHEEL_FIRST_ORDER_DYNAMICS_DENOMINATOR = [REACTION_WHEEL_FIRST_ORDER_TIME_CONSTANT 1];	%linear dynamics of the wheel, 1st order, denominator

ESTIMATOR_TRANSFER_FUNCTION_NUMERATOR  = [1 0];
ESTIMATOR_TRANSFER_FUNCTION_DENOMINATOR = [ESTIMATOR_TIME_CONSTANT 1];
ESTIMATOR_TRANSFER_FUNCTION = tf(...
    ESTIMATOR_TRANSFER_FUNCTION_NUMERATOR ,...
    ESTIMATOR_TRANSFER_FUNCTION_DENOMINATOR);

PURE_DERIVATOR_TRANSFER_FUNCTION  = tf([1 0],1);

UNITARY_GAIN = 1;

TARGET_SECOND_ORDER.NATURAL_FREQUENCY = 0.18;   %rad/s
TARGET_SECOND_ORDER.DAMPING = 0.7;
TARGET_SECOND_ORDER.TRANSFER_FUNCTION = tf([0 0 1],[TARGET_SECOND_ORDER.DAMPING^(2) 2*TARGET_SECOND_ORDER.DAMPING*TARGET_SECOND_ORDER.NATURAL_FREQUENCY 1]);

SATELLITE_INERTIA_FIDELITY_MARGIN = 0.3;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 1 Open-loop modeling and analysis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  1.1  Satellite model
% 
% c1 = -FLEXIBLE_MODE_FREQUENCY^(2)/(1 - COUPLING_COEFFICIENT);
% c2 = -2*FLEXIBLE_MODE_DAMPING_RATIO*FLEXIBLE_MODE_FREQUENCY/(1 - COUPLING_COEFFICIENT);
% c3 = 1/((1 - COUPLING_COEFFICIENT)*SATELLITE_INERTIA);
% c4 = COUPLING_COEFFICIENT/((1 - COUPLING_COEFFICIENT)*SATELLITE_INERTIA);
% 
% flexibleSatellite = StateSpaceRepresentation(...
%     [   0 1 0   0;...
%         0 0 c1  c2;...
%         0 0 0   1;...
%         0 0 c1  c2],...
%     [0; c3; 0; c4],...
%     [1 0 0 0],...
%     0);
% setSystemNames(flexibleSatellite,...
%     {'theta' 'theta_dot' 'eta' 'eta_dot'},...
%     {'T'},...
%     {'theta'});
% 
% damp(flexibleSatellite.system);
% 
% simulinkFlexibleSatellite = generateStateSpaceFromSimulink(FLEXIBLE_SATELLITE_SIMULINK_FN);
% 
% setSystemNames(simulinkFlexibleSatellite,...
%     {'theta' 'theta_dot' 'eta' 'eta_dot'},...
%     {'T'},...
%     {'theta'});
% 
% damp(simulinkFlexibleSatellite.system);
% 
% % figure,step(flexibleSatellite.system,10),grid on
% % figure,initial(flexibleSatellite.system,[0;0;1;0]),grid on
% 
% rigidSatellite = StateSpaceRepresentation(...
%     [0 1;0 0],...
%     [0;1/SATELLITE_INERTIA],...
%     [1 0],...
%     0);
% setSystemNames(rigidSatellite,...
%     {'theta' 'theta_dot'},...
%     {'T'},...
%     {'theta'});
% 
% % figure,bode(flexibleSatellite.system, rigidSatellite.system),grid on
% % figure,nichols(flexibleSatellite.system, rigidSatellite.system),grid on
% 
% rank(ctrb(flexibleSatellite.system))
% 
% rank(obsv(flexibleSatellite.system))
% 
% %%  1.2  Actuator model
% 
% %All in Simulink
% 
% 
% %%  1.3  Sensor and estimator model
% 
% % figure,bode(PURE_DERIVATOR_TRANSFER_FUNCTION , ESTIMATOR_TRANSFER_FUNCTION),grid on
% 
% %%  1.4  Complete open-loop model
% 
% %All in Simulink
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %2 CLOSED-LOOP CONTROL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%  2.1  Small pointing error
% %%  2.1.1  Modal control
% 
% modalControlOpenLoop = generateStateSpaceFromSimulink(MODAL_CTRL_OPEN_LOOP_SIMULINK_FN);
% 
% % damp(modalControlOpenLoop.system);
% 
% targetDampingCoefficient = 0.7;
% %%%%%%%%%
% % Attempt to determine the bandwith limit of the actuator by looking at
% % actuator's Bode. Not sure if legit.
% % figure,bode(tf(REACTION_WHEEL_DYNAMICS_NUMERATOR,REACTION_WHEEL_DYNAMICS_DENOMINATOR)),grid on;
% % Value kept for bandwith limit is 0.2 rad/s
% %%%%%%%%%
% targetNaturalFrequency = 0.2;
% 
% gainMatrixStateFeedback = generateModalControlGainMatrix(modalControlOpenLoop, targetDampingCoefficient, targetNaturalFrequency);
% damp(modalControlOpenLoop.A-modalControlOpenLoop.B*gainMatrixStateFeedback*modalControlOpenLoop.C);
% 
% % Margins calculated on the open loop Single Input Single Output system
% modalControlSisoOpenLoop = generateStateSpaceFromSimulink(MODAL_CTRL_SISO_OPEN_LOOP_SIMULINK_FN);
% 
% % figure,bode(modalControlSisoOpenLoop.system),grid on;
% % We read the value for gainMargin = 30.8dB at frequency 2.63rad/s
% % We read the value for phaseMargin = -79.3deg at frequency 2.23rad/s
% % The system is unstable because of negative phaseMargin (one would say no margin).
% 
% feedForwardGain = gainMatrixStateFeedback(1);
% % Controller implementation in simulink: MODAL_CTRL_CLOSED_LOOP_SIMULINK_FN
% % The maximum admissible pointing error is 3.5deg (roughly).
% % The system is suuuuuper slow to provide the correction.
% 
% maximumTolerableReferenceAngle = CUSTOM_POINTING_ERROR;   % As read on the scope from the closed loop simulink simulation
% 
% modalControlWithIntegratorOpenLoop = generateStateSpaceFromSimulink(MODAL_CTRL_WITH_INTEGRATOR_OPEN_LOOP_SIMULINK_FN);
% gainMatrixStateFeedBackWithIntegrator = generateModalControlGainMatrix(modalControlWithIntegratorOpenLoop, targetDampingCoefficient, targetNaturalFrequency);
% 
% modalControlWithIntegratorFeedForwardGain = 1;	%Obviously.
% damp(modalControlWithIntegratorOpenLoop.A-modalControlWithIntegratorOpenLoop.B*gainMatrixStateFeedBackWithIntegrator*modalControlWithIntegratorOpenLoop.C);
% 
% %% 2.1.2 Phase-lead control
% 
% phaseLeadControlSisoOpenLoopStateSpace = generateStateSpaceFromSimulink(PHASE_LEAD_CTRL_SISO_OPEN_LOOP_SIMULINK_FN);
% 
% % sisotool(phaseLeadControlSisoOpenLoopStateSpace.system);
% % There is no gain allowing reasonable positive phase and gain margin (seen in the Bode)
% % even though a gain of 28719 units achieves theoretical positive margins.
% % The root locus shows always some positive or negative quasi-zero roots,
% % therefore the system is unstable.
% 
% % We can read on the Bode diagram that GdB = 0 at frequency 0.173 rad/s.
% CROSSOVER_FREQUENCY = 0.173;
% 
% % The required phase lead is the distance between the current phase at
% % crossover and the desired phase margin (14.1 + 50).
% REQUIRED_PHASE_LEAD = (14.1 + 50)*(pi/180);
% 
% % We use the fact that in the Lead Controller phase is maximum at the
% % frequency 1/(tau*sqrt(a)) to get the second equation we need.
% solutions = generateRatioAndNaturalFrequencyPhaseLeadController(REQUIRED_PHASE_LEAD, CROSSOVER_FREQUENCY);
% 
% phaseLeadNaturalFrequencyRatio = double(solutions.symbolicPhaseLeadNaturalFrequencyRatio);
% phaseLeadDenominatorNaturalFrequency = double(solutions.symbolicPhaseLeadDenominatorNaturalFrequency);
% % barePhaseLeadGainAtCrossover = abs((1+phaseLeadNaturalFrequencyRatio*phaseLeadDenominatorNaturalFrequency*CROSSOVER_FREQUENCY*1i)/(1+phaseLeadDenominatorNaturalFrequency*CROSSOVER_FREQUENCY*1i));
% phaseLeadGain = UNITARY_GAIN;
% 
% phaseLeadControllerTransferFunctionNumerator = phaseLeadGain * [phaseLeadNaturalFrequencyRatio*phaseLeadDenominatorNaturalFrequency 1];
% phaseLeadControllerTransferFunctionDenominator = [phaseLeadDenominatorNaturalFrequency 1];
% 
% % Obviously equal to 1 here.
% phaseLeadFeedForwardGain = dcgain(tf(phaseLeadControllerTransferFunctionNumerator,phaseLeadControllerTransferFunctionDenominator));
% 
% phaseLeadControlSisoClosedLoopStateSpace = generateStateSpaceFromSimulink(PHASE_LEAD_CTRL_SISO_CONTROLLER_OPEN_LOOP_SIMULINK_FN);
% 
% %sisotool(phaseLeadControlSisoClosedLoopStateSpace.system);
% % We can read in the bode diagram with the controller in open loop that we
% % need a 0.23 gain to get back the crossover frequency on 0.173 rad/s.
% PHASE_LEAD_GAIN = 0.23;
% 
% phaseLeadControllerTransferFunctionNumerator = PHASE_LEAD_GAIN * [phaseLeadNaturalFrequencyRatio*phaseLeadDenominatorNaturalFrequency 1];
% phaseLeadControllerTransferFunction = tf(phaseLeadControllerTransferFunctionNumerator, phaseLeadControllerTransferFunctionDenominator);
% % Obviously equal to PHASE_LEAD_CONTROLLER_GAIN here.
% phaseLeadFeedForwardGain = dcgain(tf(phaseLeadControllerTransferFunctionNumerator,phaseLeadControllerTransferFunctionDenominator));
% phaseLeadFeedForwardGainTransferFunction = tf([phaseLeadFeedForwardGain]);
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%% ANALYSIS OF THE CLOSED LOOP SYSTEM %%%%%%%%%%%%%%%
% % To do so we open again the SISO open loop system and we add manually in
% % the Control System Designer the controller and the feedforward gain.
% % sisotool(phaseLeadControlSisoOpenLoopStateSpace.system);
% % We observe that the actual maximum of the induced phase lead is actually
% % a little behind the designed frequency at max. This is because of the
% % addition of the satellite model with the controller. Therefore, we can
% % make an adjustment to this.
% optimizedPhaseTarget = REQUIRED_PHASE_LEAD + (3.2 *(pi/180));	% Value read on the previous sisotool.
% optimizedFrequencyAtMaximumPhase = CROSSOVER_FREQUENCY + 0.13;  % Value read on the previous sisotool.
% phaseLeadControllerOptimizedGain = 0.335;                       % Value read on the sisotool based on the two valeus above.
% % For 45° margin use phase -1.3°, freq +0.13, gain 0.395.
% 
% % %%%%%% Interesting study on the best configuration %%%%%
% % % This controller has the best time response to a step command for margin
% % % 50°. However is it pratically achievable ?
% % % Natural freq is 9.169248080054350e-06, which doesn't appear reasonable.
% % optimizedPhaseTarget = 89.9 * (pi/180);
% % optimizedFrequencyAtMaximumPhase = CROSSOVER_FREQUENCY +  96;
% % phaseLeadControllerOptimizedGain = 0.433;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % We compute again the parameters of the phase lead controller. 
% solutions = generateRatioAndNaturalFrequencyPhaseLeadController(optimizedPhaseTarget, optimizedFrequencyAtMaximumPhase);
% 
% phaseLeadOptimizedNaturalFrequencyRatio = double(solutions.symbolicPhaseLeadNaturalFrequencyRatio);    
% phaseLeadOptimizedDenominatorNaturalFrequency = double(solutions.symbolicPhaseLeadDenominatorNaturalFrequency);
% 
% phaseLeadOptimizedControllerTransferFunctionNumerator = phaseLeadControllerOptimizedGain * [phaseLeadOptimizedNaturalFrequencyRatio*phaseLeadOptimizedDenominatorNaturalFrequency 1];
% phaseLeadOptimizedControllerTransferFunctionDenominator = [phaseLeadOptimizedDenominatorNaturalFrequency 1];
% phaseLeadOptimizedControllerTransferFunction = tf(phaseLeadOptimizedControllerTransferFunctionNumerator,phaseLeadOptimizedControllerTransferFunctionDenominator);
% 
% % Obviously equal to 1 here.
% phaseLeadOptimizedFeedForwardGain = dcgain(phaseLeadOptimizedControllerTransferFunction);
% phaseLeadOptimizedFeedForwardGainTransferFunction = tf([phaseLeadOptimizedFeedForwardGain]);
% 
% % We observe here for trials and validation.
% % sisotool(phaseLeadControlSisoOpenLoopStateSpace.system);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % sisotool(phaseLeadControllerTransferFunction);
% % 1. We can see that this phase lead controller is behaving as a high pass
% % filter at frequency 1/a*tau if plugged in direct link. 
% % 2. With unitary feedback the satellite acquires a non-convergent state
% % for a step command. The stationnary state is a sinusoïde of period 0.024s.
% % 3. Therefore, plugging the phase lead controller in direct link will
% % filter out the sinusoïde of the stationnary state making it effectively
% % converge onto the step command.
% 
% %  A pure "I" controller could bring the error to zero, but it would be
% %  both slow reacting at the start and brutal in the end.
% 
% % We can afford at most -5° at crossover frequency.
% TARGET_PI_PHASE_AT_CROSSOVER = - 5*(pi/180);
% syms symbolicPiTimeConstant;
% equationOnPiTimeConstant = -1/(symbolicPiTimeConstant*CROSSOVER_FREQUENCY) == tan(TARGET_PI_PHASE_AT_CROSSOVER);
% phaseLeadPiTimeConstant = double(solve(equationOnPiTimeConstant, symbolicPiTimeConstant));
% phaseLeadPiTransferFunctionNumerator = [phaseLeadPiTimeConstant 1];
% phaseLeadPiTransferFunctionDenominator = [phaseLeadPiTimeConstant 0];
% 
% phaseLeadControlWithPiSisoControllerOpenLoopStateSpace = generateStateSpaceFromSimulink(PHASE_LEAD_CTRL_WITH_PI_SISO_CONTROLLER_OPEN_LOOP_SIMULINK_FN);
% % sisotool(phaseLeadControlWithPiSisoControllerOpenLoopStateSpace.system);
% % We read a gain of 4.32 units.
% phaseLeadPiInducedAdditionalGain = 4.32;
% phaseLeadWithPiControllerTransferFunctionNumerator = phaseLeadPiInducedAdditionalGain*phaseLeadOptimizedControllerTransferFunctionNumerator;
% phaseLeadWithPiControllerTransferFunctionDenominator = phaseLeadOptimizedControllerTransferFunctionDenominator;
% phaseLeadWithPiControllerTransferFunction = tf(phaseLeadWithPiControllerTransferFunctionNumerator, phaseLeadWithPiControllerTransferFunctionDenominator);
% 
% phaseLeadWithPiFeedForwardGain = dcgain(phaseLeadWithPiControllerTransferFunction);
% phaseLeadWithPiFeedForwardGainTransferFunction = tf([phaseLeadWithPiFeedForwardGain]);

% %% 2.1.3 Structured Hinf Control
% 
% % differenceToSecondOrderWeightingStateSpace = makeweight(100,[0.032,0.9],0.1);    % Values obtained by trials and errors: 100,[0.032,0.9],0.1
% % torqueWeightingStateSpace = makeweight(0.1,[100,0.9],1);                         % Values obtained by trials and errors: 0.1,[100,0.9],1
% % 
% % hinfControlStateSpace = generateHinfStateSpace(HINF_CTRL_OPEN_LOOP_SIMULINK_FN);
% 
% % hinfControlNoDerivativeFeedbackStateSpace = generateHinfNoDerivativeFeedbackStateSpace(HINF_CTRL_NO_DERIVATIVE_FEEDBACK_OPEN_LOOP_SIMULINK_FN);
% % % Interestingly enough our previous weighting function work quite well here
% % % too. There is no need to find others (even though it's slightly less
% % % performant.
% 
% % differenceToSecondOrderWithIntegratorWeightStateSpace = makeweight(10,[0.1,0.1],0.01);        % Values obtained by trials and errors: 10,[0.1,0.1],0.01
% % torqueWithIntegratorWeightStateSpace = makeweight(0.1,[100,0.9],1);                           % Values obtained by trials and errors: 0.1,[100,0.9],1
% % positionWithIntegratorWeightStateSpace = makeweight(1,[0.1,0.9],0.01);                        % Values obtained by trials and errors: 1,[0.1,0.9],0.01
% % 
% % secondOrderGainCorrection = dcgain(TARGET_SECOND_ORDER.TRANSFER_FUNCTION);
% % 
% % hinfControlWithIntegratorStateSpace = generateHinfWithIntegratorStateSpace(HINF_CTRL_WITH_INTEGRATOR_OPEN_LOOP_SIMULINK_FN);
% 
% satelliteInertiaMaxWorstCase = SATELLITE_INERTIA*(1+SATELLITE_INERTIA_FIDELITY_MARGIN);
% satelliteInertiaMinWorstCase = SATELLITE_INERTIA*(1-SATELLITE_INERTIA_FIDELITY_MARGIN);
% 
% % Generation of CL0W for central value of the satellite inertia
% satelliteInertia = SATELLITE_INERTIA;
% % To determine the new weight functions one can use "generateHinfWithIntegratorStateSpace(HINF_CTRL_WITH_INTEGRATOR_MULTIMODEL_OPEN_LOOP_SIMULINK_FN)"
% % to observe the impact of the new inertia "satelliteInertia and the impact
% % of the new "*MultimodalWeightStateSpace" functions.
% differenceToSecondOrderWithIntegratorMultimodalWeightStateSpace = makeweight(10,[0.1,0.1],0.01);    % Same values as before of course: 10,[0.1,0.1],0.01.
% torqueWithIntegratorMultimodalWeightStateSpace = makeweight(0.1,[100,0.9],1);                       % Same values as before: 0.1,[100,0.9],1.
% positionWithIntegratorMultimodalWeightStateSpace = makeweight(1,[0.1,0.9],0.01);                    % Same values as before: 1,[0.1,0.9],0.01.
% CL0W_central = generateHinfWithIntegratorCl0wMatrix(HINF_CTRL_WITH_INTEGRATOR_MULTIMODEL_OPEN_LOOP_SIMULINK_FN);
% % Generation of CL0W for max worst case value of the satellite inertia
% satelliteInertia = satelliteInertiaMaxWorstCase;
% differenceToSecondOrderWithIntegratorMultimodalWeightStateSpace = makeweight(10,[0.1,0.057],0.01);  % Values determined by trials and validations: 10,[0.1,0.057],0.01.
% CL0W_maxWorstCase = generateHinfWithIntegratorCl0wMatrix(HINF_CTRL_WITH_INTEGRATOR_MULTIMODEL_OPEN_LOOP_SIMULINK_FN);
% % Generation of CL0W for min worst case value of the satellite inertia
% satelliteInertia = satelliteInertiaMinWorstCase;
% differenceToSecondOrderWithIntegratorMultimodalWeightStateSpace = makeweight(10,[0.1,0.225],0.01);  % Values determined by trials and validations: 10,[0.1,0.225],0.01.
% CL0W_minWorstCase = generateHinfWithIntegratorCl0wMatrix(HINF_CTRL_WITH_INTEGRATOR_MULTIMODEL_OPEN_LOOP_SIMULINK_FN);
% 
% % Generation of the multimodal HinfController, optimum coefficients are:
% % 0.4, 0.1, 0.25
% % We supposed we could make a frequency domain analysis to tune the model
% % to the frequency spectrum of inputs. Here it wouldn't be usefull.
% hinfControlWithIntegratorMultimodelStateSpace = generateHinfWithIntegratorMultimodelStateSpace(0.4,0.1,0.25, CL0W_central, CL0W_maxWorstCase, CL0W_minWorstCase);

%% Function definitions

function gainMatrix = generateModalControlGainMatrix(spaceSateRepresentation, secondOrderDampingCoefficient, secondOrderNaturalFrequency)
    n = size(spaceSateRepresentation.A,1);
    lambda=roots([1 2*secondOrderDampingCoefficient*secondOrderNaturalFrequency  secondOrderNaturalFrequency ^2]);
    VW1=null([spaceSateRepresentation.A-lambda(1)*eye(n) spaceSateRepresentation.B]);
    VW2=null([spaceSateRepresentation.A-lambda(2)*eye(n) spaceSateRepresentation.B]);
    V=[VW1(1:n) VW2(1:n)];
    W=[VW1(n+1) VW2(n+1)];
    
    gainMatrix=-real(W/(spaceSateRepresentation.C*V));
end

function stateSpaceRepresentation = generateStateSpaceFromSimulink(simulinkFileName)
    [A,B,C,D] = linmod(simulinkFileName);
    stateSpaceRepresentation = StateSpaceRepresentation(A,B,C,D);
end

function solutions = generateRatioAndNaturalFrequencyPhaseLeadController(phaseTarget, targetFrequencyForMaximumPhase)
    syms symbolicPhaseLeadNaturalFrequencyRatio symbolicPhaseLeadDenominatorNaturalFrequency;
    equationOnLeadControllerFrequency = targetFrequencyForMaximumPhase == 1 / (symbolicPhaseLeadDenominatorNaturalFrequency * sqrt(symbolicPhaseLeadNaturalFrequencyRatio));
    equationOnLeadControllerPhaseTarget = sin(phaseTarget) == (symbolicPhaseLeadNaturalFrequencyRatio - 1)/(symbolicPhaseLeadNaturalFrequencyRatio + 1);
    equations = [equationOnLeadControllerFrequency equationOnLeadControllerPhaseTarget];
    
    solutions = solve(equations, [symbolicPhaseLeadNaturalFrequencyRatio symbolicPhaseLeadDenominatorNaturalFrequency]);
end

function lowpassFilter = initializeLowPassFilter(gain, timeConstant)
    lowpassFilter.gain = gain;
    lowpassFilter.timeConstant = timeConstant;
    lowpassFilter.transferFunction = tf([0 gain], [timeConstant 1]);
end

function controllerStateSpace = generateHinfStateSpace(simulinkFileName)
    [a,b,c,d]= linmod(simulinkFileName);
    P=ss(a,b,c,d);
    ni=3;
    no=1;                                       % number  of  inputs  and  outputs  of the  controller
    K0=ltiblock.gain('K0',no,ni);               % controller  initialization: no-by-ni gain
    CL0=lft(P,K0);                              % parametric  closed -loop (genss) to be  optimized
    W1=1;                                       % weighting  functions  for z1, it's actually implemented in the simulink directly
    W2=1;                                       % weighting  functions  for z2, it's actually implemented in the simulink directly
    CL0W=blkdiag(W1,W2)*CL0;                    % weighted  standard  form
    opt=hinfstructOptions('RandomStart' ,3);	% 3 restarts (nonconvex  problem)
    [CL ,normHinf ]= hinfstruct(CL0W ,opt);     % call to  hinfstruct
    controllerStateSpace = ss(CL.Blocks.K0);              % final  controller  value
end

function controllerStateSpace = generateHinfNoDerivativeFeedbackStateSpace(simulinkFileName)
    [a,b,c,d]= linmod(simulinkFileName);
    P=ss(a,b,c,d);
    ni=2;
    no=1;                                       % controller  dimensions (now 2 inputs  instead  of 3)
    ns=1;                                       % number  of  states  of the  controller (first  order  here)
    K0=ltiblock.ss('K0',ns,no,ni);              % init: no-by-ni ss  object  with ns  states
    CL0=lft(P,K0);                              % parametric  closed -loop (genss) to be  optimized
    W1=1;                                       % weighting  functions  for z1, it's actually implemented in the simulink directly
    W2=1;                                       % weighting  functions  for z2, it's actually implemented in the simulink directly
    CL0W=blkdiag(W1,W2)*CL0;                    % weighted  standard  form
    opt=hinfstructOptions('RandomStart' ,6);	% 3 restarts (nonconvex  problem)
    [CL ,normHinf ]= hinfstruct(CL0W,opt);      % call to  hinfstruct
    controllerStateSpace = ss(CL.Blocks.K0);    % final  controller  value
end

function controllerStateSpace = generateHinfWithIntegratorStateSpace(simulinkFileName)
    [a,b,c,d]= linmod(simulinkFileName);
    P=ss(a,b,c,d);
    ni=4;
    no=1;
    K0=ltiblock.gain('K0',no,ni);
    CL0=lft(P,K0);
    W1=1; 
    W2=1; 
    W3=1;
    CL12=blkdiag(W1,W2)*CL0 (1:2 ,1);
    CL3=W3*CL0(3,2);
    CL0W=blkdiag(CL12 ,CL3);
    opt=hinfstructOptions('RandomStart',3);
    [CL ,normHinf ]= hinfstruct(CL0W ,opt);
    
    controllerStateSpace = ss(CL.Blocks.K0);
end

function CL0W = generateHinfWithIntegratorCl0wMatrix(simulinkBaseModelFileName)
    [a,b,c,d]= linmod(simulinkBaseModelFileName);
    P=ss(a,b,c,d);
    ni=4;
    no=1;
    K0=ltiblock.gain('K0',no,ni);
    CL0=lft(P,K0);
    W1=1; 
    W2=1; 
    W3=1;
    CL12=blkdiag(W1,W2)*CL0 (1:2 ,1);
    CL3=W3*CL0(3,2);
    CL0W=blkdiag(CL12 ,CL3);
end

function controllerStateSpace = generateHinfWithIntegratorMultimodelStateSpace(r0, r1, r2, CL0W_central,CL0W_wc1 ,CL0W_wc2)
    CL0W=blkdiag(r0*CL0W_central ,r1*CL0W_wc1 ,r2*CL0W_wc2 );
    opt=hinfstructOptions('RandomStart',3);
    [CL,normHinf ]= hinfstruct(CL0W ,opt);
    
    controllerStateSpace = ss(CL.Blocks.K0);
end

%%%%%%%%%
% Attempt to automate state space representation generation for augmented loop.
% It works but the newer method is more general than this particular
% example.
%
% function stateSpaceRepresentation = generateAugmentedStateSpaceFromSimulink(simulinkFileName)
%     [A,B,C,D]=linmod(simulinkFileName);
%     Aaug = [A zeros(size(A,1),1);...
%             -1 zeros(1,size(A,1))];
%     Baug = [B; 0];
%     Caug = [zeros(1,size(C,2)) -1;...
%             -C zeros(size(C,1),1)];
%     Daug = [0; -D];
%     stateSpaceRepresentation = StateSpaceRepresentation(Aaug,Baug,Caug,Daug);
% end
%%%%%%%%%