% MainModel.m Simulation of a multicellular Hes5 dynamics coupled via lateral inhibition. 
% Number of cells can be defined by rows and cols parameters,
% and will produce a hexagonal grid of cells that each interact with the
% six surrounding cells, with the option of periodic or non-periodic 
% boundary conditions.

clear;clc;
addpath('Functions')



%% Choose the colour of all graphs
GraphAppearance = 0;   % 0 = Normal, 1 = Dark grey

if GraphAppearance==1   
    BackgroundColour=38/255*[1 1 1]; TextColour=[1 1 1];
    get(0,'Factory');                           set(0,'defaultfigurecolor',BackgroundColour)
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor',BackgroundColour)
    set(0,'defaultAxesXColor',TextColour);      set(0,'defaultAxesYColor',TextColour)
    set(0,'defaultLegendTextColor',TextColour); set(0,'defaultTextColor',TextColour)  
    set(0,'defaultAxesZColor',TextColour)
elseif GraphAppearance==0
    get(0,'Factory');                           set(0,'defaultfigurecolor',[0.94 0.94 0.94])
    set(0,'DefaultAxesFontSize', 12);           set(0,'defaultAxesColor','w')
    set(0,'defaultAxesXColor','k');             set(0,'defaultAxesYColor','k')
    set(0,'defaultLegendTextColor','k');        set(0,'defaultTextColor','k')
end

%% Define grid size
rows=26;         %Number of rows of cells in hexagonal grid
cols=6;         %Number of columns of cells in hexagonal grid
cells=rows*cols; %Total number of cells

%% Define simulation time and step size
t0=0;            % Start time
tf_hours=200;    % Final time (hours)
tf=tf_hours*60;  % Final time (min)

Stochastic=1;              %1 = stochastic, 0 = deterministic  (deterministic currently has a bug) 
dt=2;            % Time step size (min)
if Stochastic==1
    dt=1;        % Stochastic uses Euler method rather than Runge-Kutta, so this requires a smaller step size      
end
Nt=(tf-t0)/dt;   % Number of time elements
T=t0:dt:tf;      % Time vector

%% Set up simulation inputs
% For each of the following assingments, 1=yes, 0=no.

%Important options!
CoupledCells=1;            %Simulate with Notch-Delta coupling = 1
Boundary=0;                %0 = 'hard' boundary, 1 = periodic boundary (set to 1)
TurnOffAutorepression=0;   %Reduce model to just lateral inhibition without autonomous Hes oscillators in each cell
VertProtrusions=0;         %Increases number of signalling neighbours in the vertical direction 
HorzProtrusions=0;
HorzProtStrength=1;

%Differentiation selection and parameters
CrudeDiff           = 0;   %Cells will be marked...
AvM                 = 0;   %Absolute vs moving mean threshold (0=Absolute thresh, 1=Moving mean thresh)
wl_mins             = 100; %Window length in mins to take the moving mean from
wl=wl_mins/dt;             %Window length converted to number of vector elements
Replace             = 0;   %Replace differentiating cells with mean population protein
DiffTime_hours=50;         %Time at which differention can start to occur (hours)
DiffTime=DiffTime_hours*60/dt;
S=0.01;                    % Rate of differentiation. Nominal value of 0.01

%Cell-movement/swapping parameters
ImplementSwapping   = 0;   %Cells will randomly swap in the horizontal direction
Pm=0.005; % Total probability of movement either left or right (max Pm value is 0.5)
SwapThresh=1-Pm/2; % Treshold to be used in probability of swapping function within diffSolver.m
Ts=5; % Time (in mins) between each swapping event (nominal Ts=5)       

if cols==1
    kymoCols=1;            % Used later for plotting kymographs        
else
    kymoCols=2;
end

%% Select simulation outputs
%Various frequency analysis options
TemporalFourier     = 1;   %Gives average Fourier period of the cell population
TemporalWavelet     = 0;   %Preliminary implementation of wavelet to examine whether Hes5 switches between periodic and noisy
SpatialFourier      = 1;   %Detect significant spatial periodic expression of Hes5

%Visualising Hes5 protein levels
AnimateGrid         = 1;   
AnimationSpeed      = 1; %Multiplier for speed of animation
AnimationSubplots   = 1;   %0=Just Hes levels plot, 1=show swapping, 2=show crude differentiation
ShowRandomCells     = 1;   %Time traces of individual cells
ShowLastFrame       = 1;   %Shows the last time point in the hexagonal lattice arrangement of cells
AnimatePhaseSpace   = 0;   %Animate Hes protein levels vs mRNA levels over time

%Phase analysis of oscillatory signals
PerformMFV          = 0;   % Calculate mean field value (degree of synchronisation)
AnimateComplexPhase = 0;   % Requires PerformMFV to be computed as well

if AnimateComplexPhase==1
    PerformMFV=1;
end


%% Produce neighbour matrix

if cells==1
    Boundary=0; % Other boundaries don't make sense for 1 cell!
end

[NM,NumNeigh]=neighbours(rows,cols,Boundary,VertProtrusions,HorzProtrusions,HorzProtStrength); 
NM=sparse(NM); % Saves a lot on computational cost for large grid sizes!

eps=(1./NumNeigh); % Make this as an output of neighbours.m!
if cells==1
    eps=0;
end

%% Model parameters
 
% Columns in accepted_parameters: 1:a_m, 2:a_p, 3:P_H0, 4:TauH, 5:n_H
load accepted_parameters %Selcted for coherence values between 0.05 and 0.015
load summary_stats

%Selecting parameter sets that give the following single cell summary
%stats:
Select(:,1)=accepted_parameters(:,4)>15;                        %TauH selection
Select(:,2)=accepted_parameters(:,5)<3;                         %nH selection
Select(:,3)=summary_stats(:,4)<0.13 & summary_stats(:,4)>0.12;  %Coherence selection
Select(:,4)=summary_stats(:,1)<5e4 & summary_stats(:,1)>3.5e4;  %Mean selection
Select(:,5)=summary_stats(:,2)>0.05 & summary_stats(:,2)<0.1;   %Coefficient of variation                         %Standard deviation selection
Select(:,6)=summary_stats(:,3)<10*60;                           %Peak period oscillation selection 

accept=sum(Select,2);
accept=find(accept==size(Select,2));
fprintf(sprintf('Number of accepted parameter sets = %.f \n',numel(accept)))

PlotSummaryStats=0;
if PlotSummaryStats==1
    
    figure(21)
    clf
    
    subplot(231); histogram(accepted_parameters(accept,4)); title('TauH')
    subplot(232); histogram(accepted_parameters(accept,5)); title('n_H')
    subplot(233); histogram(summary_stats(accept,4));       title('Coherence')
    subplot(234); histogram(summary_stats(accept,1));       title('Mean')
    subplot(235); histogram(summary_stats(accept,2));       title('Coefficient of variation')
    subplot(236); histogram(summary_stats(accept,3)/60);    title('Period')

    figure(22)
    clf
    plot(summary_stats(:,3)/60,summary_stats(:,4),'.','color','w');hold on
    plot(summary_stats(accept,3)/60,summary_stats(accept,4),'.','color',[0.2 0.9 0.4])
    xlabel('Period (hours)')
    ylabel('Coherence')
end


%% Choose single cell parameter set
%Select a parameter set from the Bayesian inferred single cell parameters 
%from (Manning et al. 2019) 
AP=accepted_parameters;
SS=summary_stats;

% s=197; %0.05 Coh
% s=663; %0.05 Coh
% s=987; %0.05 Coh

% s=300; %0.1 Coh
% s=350; %0.1 Coh
% s=371; %0.1 Coh

s=3700; %0.15 Coh
% s=616; %0.15 Coh
% s=713; %0.15 Coh

%__________________________________________________________________________
%Single-cell parameters
a_m   = AP(s,1); 
a_p=AP(s,2); 
P_H0=AP(s,3); 
TauH=AP(s,4); 
% TauH=0;
n_H=AP(s,5); 
u_m   = log(2)/30; 
u_p   = log(2)/90;

% %High coherence option
% Prmt=[58.5,  1.18,  27378, 6,    40,   1000000, 1,   169,  log(2)/30, log(2)/90]; %Coher=0.96, very weakly coupled, s of antiphase
% a_m   = Prmt(1); a_p   = Prmt(2); P_H0  = Prmt(3); n_H   = Prmt(4);     
% TauH  = Prmt(5); P_ND0 = Prmt(6); n_ND  = Prmt(7); TauND = Prmt(8);
% u_m   = Prmt(9); u_p   = Prmt(10);

%Multicellular parameters
% P_ND0 = linspace(6000,20000,cells)'; 




% u_m = linspace(log(2)/10,log(2)/50,rows)'; 
% % u_m = [linspace(0.2,50,rows/2)';linspace(50,0.2,rows/2)']; %Highest in the middle gradient
% [~,Y]=meshgrid(1:cols,u_m);
% u_m=Y(:);

% u_p = linspace(log(2)/60,log(2)/120,rows)'; 
% % u_p = [linspace(0.2,50,rows/2)';linspace(50,0.2,rows/2)']; %Highest in the middle gradient
% [~,Y]=meshgrid(1:cols,u_p);
% u_p=Y(:);


% P_ND0 = linspace(12000,6000,rows)'; 
% P_ND0 = [linspace(7000,20000,rows/2)';linspace(20000,7000,rows/2)']; %Highest in the middle gradient
% [~,Y]=meshgrid(1:cols,P_ND0);
% P_ND0=Y(:);


na=0;   % Noise amplitude - this is an absolute value, i.e. if the gradient is in TauH, then the na is in TauH units (minutes).
nst=Nt;   % Noise switching time (min)
lines=0;  % Whole lines of same noise applied = 1, if not = 0
dynamic=1;

minT=10; %Min TauH value in the gradient
maxT=110;%Max TauH value in the gradient
TauH = linspace(minT,maxT,rows)'; 

% gn=na*randn(size(TauH));
% TauH=TauH+gn;
% TauH(TauH<minT)=minT;
% TauH = [linspace(minT,maxT,rows/2)';linspace(maxT,minT,rows/2)'];
[~,Y]=meshgrid(1:cols,TauH);
TauH=Y(:);


% row1=rows/2;
% row2=row1+5;
% 
% % d=Y(row1,1)-Y(row1-1,1);
% % % Y(row1:row2,:)=Y(row1:row2,:)-2*d;
% % % Y(row1:row2,:)=flipud(Y(row1:row2,:));
% % Y(row1:row2,:)=Y(row1,1)+10*ones(row2-row1+1,cols)*d;
% Y(row1:row2,:)=10;
% % idx=[1:5:rows,2:5:rows];
% % % idx=[1:6:rows,2:6:rows,3:6:rows];
% % Y(idx,:)=10;





TauND = linspace(90,180,rows)';
% TauND = [linspace(90,180,rows/2)';linspace(180,90,rows/2)'];
[~,Y]=meshgrid(1:cols,TauND);
TauND=Y(:);

% P_H0 = linspace(2000,30000,rows)'; 
% % P_H0 = [linspace(8000,40000,rows/2)';linspace(40000,8000,rows/2)']; %Highest in the middle gradient
% [~,Y]=meshgrid(1:cols,P_H0);
% P_H0=Y(:);

% a_m = linspace(0.5,5,rows)'; 
% % a_m = [linspace(0.2,50,rows/2)';linspace(50,0.2,rows/2)']; %Highest in the middle gradient
% [~,Y]=meshgrid(1:cols,a_m);
% a_m=Y(:);

% a_p = linspace(5,50,rows)'; 
% % a_p = [linspace(5,50,rows/2)';linspace(50,5,rows/2)']; %Highest in the middle gradient
% [~,Y]=meshgrid(1:cols,a_p);
% a_p=Y(:);

% n_H = linspace(1,6,rows)'; 
% % n_H = [linspace(2,6,rows/2)';linspace(6,2,rows/2)'];
% [~,Y]=meshgrid(1:cols,n_H);
% n_H=Y(:);

% n_ND = linspace(3,6,rows)'; 
% % n_ND = [linspace(2,6,rows/2)';linspace(6,2,rows/2)'];
% [~,Y]=meshgrid(1:cols,n_ND);
% n_ND=Y(:);

% P_ND0 = linspace(8000,16000,rows)'; 
% [~,Y]=meshgrid(1:cols,P_ND0);
% P_ND0=Y(:);

%Add noise to coupling strength gradient
% gn=500*randn(size(P_ND0));
% P_ND0=P_ND0+gn;

% % % % % % gn=30*randn(size(TauH));
% % % % % % TauH=TauH+gn;
% % % % % % TauH(TauH<minT)=minT;

% gn=0.3*randn(size(TauH))+1;
% gn=1.5*(rand(size(TauH))-0.5)+1;
% gn(gn<0)=0;

% TauH=TauH.*gn;
% TauH(TauH<minT)=minT;
% TauH(TauH>maxT)=maxT;
% figure(1111)
% clf
% histogram(gn)

% TauH(TauH>
% min(TauH)
% max(TauH)
gradInput = Y(:);

P_ND0 = 12000; 
n_ND  = 4; 
% TauND = 140;


%Saving animation options
MakeGIF   = 0;          %1=yes, 0=no. Any animations that run will be made into GIFs.
filename1 = 'Travelling Waves'; %Specify the output file name
filename2 = 'GIF2.gif'; %Specify the output file name
filename3 = 'GIF3.gif'; %Specify the output file name
filename4 = 'GIF4_no_coupling.gif';

%__________________________________________________________________________

gamma=1; % Maximum intercellular Hill function value

TauH_step=round(TauH/dt);   % Conversion to simulation time steps
TauND_step=round(TauND/dt); % Conversion to simulation time steps

if CoupledCells==0
    eps=0; %Set eps=0 to decouple the cells 
end

if CrudeDiff==0
    DiffTime=Nt*2;
end

Up=1;
Low=0;
a=Low;
b=Up-a;

%% Functions/DDEs 
% %Using anonymous functions (@ symbol) rather than standard functions
% dm=@(m,p,p_delay,Psum_delay,gamma) a_m.*1./(1+(p_delay./P_H0).^n_H).*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) - u_m.*m; % Describe mRNA change in time
% 
% if TurnOffAutorepression==1
%         dm=@(m,p,p_delay,Psum_delay,gamma) a_m.*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) - u_m.*m; % No self repression version for simple Notch Delta modelling
% end
% 
% dp=@(m,p,p_delay,Psum_delay,gamma) a_p.*m - u_p.*p; % Describes protein change in time
% 
% %SDDEs_____________________________________________________________________
% 
% dm1=dm;
% dm2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_m.*1./(1+(p_delay./P_H0).^n_H).*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) + u_m.*m);
% 
% if TurnOffAutorepression==1
%    dm2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_m.*gamma./(1+(eps.*Psum_delay./P_ND0).^n_ND) + u_m.*m);
% end
% 
% dp1=dp;
% dp2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_p.*m + u_p.*p);
EPS=1;
PSUM_DELAY=linspace(1,10000);
HILL=a + b./(1+(EPS.*PSUM_DELAY./P_ND0).^n_ND);
figure(23)
clf
plot(PSUM_DELAY,HILL)
ylim([0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DDEs (using anonymous functions (@ symbol) rather than standard functions)
dm=@(m,p,p_delay,Psum_delay,gamma) a_m*1./(1+(p_delay./P_H0).^n_H).*gamma.*(  a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)  ) - u_m.*m; % Describe mRNA change in time
                                                                              

if TurnOffAutorepression==1
    dm=@(m,p,p_delay,Psum_delay,gamma) a_m.*gamma.* (a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND)) - u_m.*m; % No self repression version for simple Notch Delta modelling
end

dp=@(m,p,p_delay,Psum_delay,gamma) a_p.*m - u_p.*p; % Describes protein change in time

%% SDDEs
dm1=dm;

dm2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_m.*gamma./(1+(p_delay./P_H0).^n_H).*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);

if TurnOffAutorepression==1
   dm2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_m.*gamma.*( a + b./(1+(eps.*Psum_delay./P_ND0).^n_ND) ) + u_m.*m);
end

dp1=dp;
dp2=@(m,p,p_delay,Psum_delay,gamma) sqrt(a_p.*m + u_p.*p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialise vectors and initial values (random)____________________________

% rng(3); % Same seed for random number generator (prevents initial heterogeneities from changing on every simulation, so comment out this line if you want different initial conditions each time)
P=[rndrng(cells,1,0.2*P_H0(1),1.8*P_H0(1)),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
M=[rndrng(cells,1,0,20),zeros(cells,Nt-1),zeros(cells,1)];              %Vector that will store mRNA values from the simulation

% P=[rndrng(cells,1,0.2.*P_H0,1.8*P_H0),zeros(cells,Nt-1),zeros(cells,1)]; %Vector that will store protein values from the simulation
% P(rows*cols/2+1:(rows+1)*cols/2,Nt/2:Nt/2+100)=100000;
% P(rows*cols/2+1:(rows)*(cols/2+1),1)=100000;
% [~,Y]=meshgrid(1:rows,1:cols);
% Y(1)=60000;
% P=[Y(:)/cols*600,zeros(cells,Nt-1),zeros(cells,1)];

%% Main loop (solves differential equations)
PercentCounter=1; %Print progress of the differential solver to the command line (1=yes, 0=no)
[P, M, t_step, DiffYN, DiffYNflash, CT, DiffThresh, MM, HILL, grad]=diffSolver(P, M, Nt, TauND_step, TauH_step, NM, gamma, dm, dp, dm1, dm2, dp1, dp2, dt, Stochastic, rows, cols, DiffTime, S, ImplementSwapping, SwapThresh, Ts, PercentCounter, AvM, Replace, wl,a,b,P_ND0,n_ND,minT,maxT,na,nst,lines,dynamic, gradInput);
%     figure(103)
%     scatter(T(2:end),HILL,5,'filled');hold on
%     ylim([0 1])


% P=randn(size(P,1), size(P,2)) % Testing white noise on Fisher-G test

if dynamic==1
    grad=[grad(:,1), grad];
else
    grad=repmat(TauH,[1,Nt+1]);
end
longExposure=DiffYNflash;
DiffElems=find(DiffYNflash==1);
ExposureTime=50;
for E=1:ExposureTime
    Elems=DiffElems+cells*E;
    Elems(Elems>numel(longExposure))=[];
    longExposure(Elems)=1;
end

if eps==0
    fprintf('\nCells are uncoupled!\n')
end
if TurnOffAutorepression==1
    fprintf('\nAutorepression is turned off!\n')
end
%%
if ImplementSwapping==1
    add=meshgrid(0:cells:cells*(Nt+1)-cells, 1:cells);
    [~,sortMatrix]=sort(CT);
    recoverCellPos=sortMatrix+add;

    P_recover=P(recoverCellPos); %Recovering cell positions so that individual cell dynamics can be plotted

    cellPos=cellPosition(CT,rows,cols);
    cellPosZero=cellPos-repmat(cellPos(:,1),1,Nt+1);

    Dist=cellPos(:,end)-cellPos(:,1);
    expectedRootMeanSquare=sqrt((1-SwapThresh)*tf/Ts);
    rootMeanSquare=sqrt(mean(Dist.^2));

    figure(22123)
    clf
    subplot(221)
    plot(T/60,cellPos(1:20:cells,:))
    ylabel('Distance (cells)')
    xlabel('Time (hours)')
    title('Trajectories of cells')
    subplot(222)
    plot(T/60,cellPosZero,'w','linewidth',0.05)
    ylabel('Distance (cells)')
    xlabel('Time (hours)')
    title({'Trajectories of cells','starting from zero'})


    subplot(223)
    histogram(Dist)
    xlim([min(Dist) max(Dist)])
    xlabel('Distance (cells)')
    ylabel('Frequency')
    title({'Distance moved', 'from inital position'})
    subplot(224)
    histogram(sqrt(Dist.^2));hold on
    scatter(rootMeanSquare,0,'r','linewidth',5)
    xlim([0 max(Dist)])
    xlabel('Distance (cells)')
    ylabel('Frequency')
    title(sprintf('Expected RMS distance= %.1f cells \nSimulation RMS distance = %.1f cells',expectedRootMeanSquare, rootMeanSquare))
    drawnow

    avg_dist=zeros(1,cols);
    n_arr=zeros(1,cols);
    
    for n=1:cols
        avg_dist(n)=mean(abs(cellPosZero(cellPos(:,1)==n,end)));
        n_arr(n)=n;
    end
    figure(9786)
    plot(n_arr,avg_dist)
end

%==========================================================================
%%                        Animation of Grid
%==========================================================================

if AnimateGrid==1
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    
    colour_index1=reshape(flipud(vecTOmat(P(:,1),cols)),[1,cols*rows]);
    
    figure(1)
    clf;
    fig=gcf;
    fig.InvertHardcopy = 'on';

%     mapRB=redblue(1000);
    map=inferno(1000);
    if AnimationSubplots >0
        
        hex1=subplot(121);
        colormap(hex1,map); 
    else
        hex1=subplot(121);
        colormap(hex1,map); 
        colormap(map);
    end
    
    if AnimationSubplots == 3
        hex1=subplot(211);
    end
    
    hexagons1 = patch(X_vertices,Y_vertices,colour_index1,'edgecolor','none');
    set(gca,'xtick',[],'ytick',[]); 
    
    if AnimationSubplots~=3
        colorbar; 
    end
    axis equal; title('Hes expression')
    set(gca,'Visible','off')
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'fontsize',15)
    caxis([min(min(P(:,round(0.2*Nt):end))) max(max(P(:,round(0.2*Nt):end)))])
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    
    
    
    if AnimationSubplots==1
%         colour_index2=reshape(flipud(vecTOmat(CT(:,1),cols)),[1,cols*rows]);
        colour_index2=reshape(flipud(vecTOmat(grad(:,1),cols)),[1,cols*rows]); %Dynamic plot
%         colour_index2=gradStatic;
        hex2=subplot(122);
        hexagons2 = patch(X_vertices,Y_vertices,colour_index2,'edgecolor','none');
        hexagons2 = patch(X_vertices,Y_vertices,colour_index2);
        set(gca,'xtick',[],'ytick',[]); 
%         colormap(hex2,randColourMapSparse(cells,0.1));
        colormap(hex2,gray)
        colorbar; axis equal; title('Cell tracking')
        set(gca,'Visible','off')
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        set(gca,'fontsize',10)
        xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
        ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
        caxis([min(min(grad(:,round(0.2*Nt):end))) max(max(grad(:,round(0.2*Nt):end)))+1]) %Dynamic
%         caxis([min(gradStatic) max(gradStatic)])
        
    elseif AnimationSubplots==2
        P_exposure=P;
        P_exposure(longExposure==1)=max(max(P(:,round(0.2*Nt):end)))+1;
        exposureColorMap=[viridis(1000); ones(1,3)];
        
%         colour_index3=reshape(flipud(vecTOmat(DiffYN(:,1),cols)),[1,cols*rows]);
%         colour_index3=reshape(flipud(vecTOmat(longExposure(:,1),cols)),[1,cols*rows]);
        colour_index3=reshape(flipud(vecTOmat(P_exposure(:,1),cols)),[1,cols*rows]);
        
        hex3=subplot(122);
        hexagons3 = patch(X_vertices,Y_vertices,colour_index3,'edgecolor','none');
        set(gca,'xtick',[],'ytick',[]); 
%         colormap(RWG(0)); 
%         colormap(hex3,RWG(0));
        colormap(hex3,exposureColorMap);
        colorbar; axis equal; title('Crude diff')
        set(gca,'Visible','off')
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        set(gca,'fontsize',10)
        caxis([min(min(P(:,round(0.2*Nt):end))) max(max(P(:,round(0.2*Nt):end)))+1])
        xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
        ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
        
        elseif AnimationSubplots==3
            subplot(212)
            ymin=min(min(P(:,Nt/8:Nt)));
            ymax=max(max(P(:,Nt/8:Nt)));


            p1=plot(P(1:rows:rows*cols,t_step),'w'); hold on
            p2=plot(P(1:rows:rows*cols,t_step),'y.'); hold off
            title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
            xlabel('Cell number')
            ylabel('Hes5 abundance')
            ylim([ymin ymax])
            xlim([1 cols])

    end
      
    startT=Nt*0.5;
%     startT=Nt;
%     T_STEP=startT:round(AnimationSpeed*Nt/400):Nt;
    T_STEP=1:round(AnimationSpeed*Nt/400):Nt;
%     T_STEP=1;
    ti=length(T_STEP);
    im=struct([]);
    
    for idx=1:ti
        
        t_step=T_STEP(idx);
        
        colour_index1=reshape(flipud(vecTOmat(P(:,t_step),cols)),[1,cols*rows])';
        set(hexagons1, 'FaceVertexCData',colour_index1); 
        subplot(121)
        title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));

        
        if AnimationSubplots==1
       
%             colour_index2=reshape(flipud(vecTOmat(CT(:,t_step),cols)),[1,cols*rows])';
            colour_index2=reshape(flipud(vecTOmat(grad(:,t_step),cols)),[1,cols*rows])';
            set(hexagons2, 'FaceVertexCData',colour_index2); 
            subplot(122)
            title('Parameter gradient')

        elseif AnimationSubplots==2
%             colour_index3=reshape(flipud(vecTOmat(DiffYN(:,t_step),cols)),[1,cols*rows])';
%             colour_index3=reshape(flipud(vecTOmat(longExposure(:,t_step),cols)),[1,cols*rows])';
            colour_index3=reshape(flipud(vecTOmat(P_exposure(:,t_step),cols)),[1,cols*rows])';
            set(hexagons3, 'FaceVertexCData',colour_index3); 
%             title('Differentiation events');
            title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
        elseif AnimationSubplots==3
            set(p1, 'yData',P(1:rows:rows*cols,t_step));
            set(p2, 'yData',P(1:rows:rows*cols,t_step));
            subplot(212)
            title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
        end
        drawnow;  

        if MakeGIF==1
            F1=getframe(gcf); %gca makes frame without labels, gcf makes frame including labels
            im{idx}=frame2im(F1);
            F_all{idx}=F1;
        end

    end
    
    if MakeGIF==1
        IDX=idx;
        for idx = 1:IDX
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,append(filename1,'.gif'),'gif','LoopCount',Inf,'DelayTime',1/30);
            else
                imwrite(A,map,append(filename1,'.gif'),'gif','WriteMode','append','DelayTime',1/30);
            end
        end
        
        video = VideoWriter(filename1); %create the video object
        video.FrameRate = 30;
        open(video); %open the file for writing
        for idx = 1:IDX
        %     I = imread(im{idx}); %read the next image
            writeVideo(video,F_all{idx}); %write the image to file
        end
        close(video); %close the file
    end

end


% %% Plot of travelling wave
% ymin=min(min(P(:,Nt/8:Nt)));
% ymax=max(max(P(:,Nt/8:Nt)));
%     
% T_STEP=1:round(AnimationSpeed*Nt/400):Nt;
% ti=length(T_STEP);
%     
% figure(10)
% clf
% for idx=1:ti
%     t_step=T_STEP(idx);
%     plot(P(1:rows:rows*cols,t_step),'w'); hold on
%     plot(P(1:rows:rows*cols,t_step),'y.'); hold off
%     title(sprintf('Time: %.0f/%.0f hours',t_step*dt/60, tf/60));
%     xlabel('Cell number')
%     ylabel('Hes5 abundance')
%     ylim([ymin ymax])
%     drawnow
% end

%%
if CrudeDiff==1

    %%
     rand_ind=ceil(cells*rand(1,ceil(cells*0.1)));
    
    if cells==2
        rand_ind=[1 2];
    end
    
 
    T_plot=T(DiffTime+1:end);
    cell_num=5;
    P_plot=P(cell_num,DiffTime+1:end);
    DT=DiffThresh(cell_num,end);
    C=P_plot-DT;
    C(C>0)=0;
    C=C/min(C);
%     C=C/max(C)
    
    
    load 'cmap.mat'
    
    figure(102)
    clf
    set(gcf,'renderer','Painters')
    
    plot(T_plot./60,DiffThresh(cell_num,DiffTime+1:end),'--','color','k','linewidth',1.5);hold on
    scatter(T_plot./60,P_plot,12,C,'filled');hold on
    xlim([DiffTime_hours tf_hours])
    ylabel('Hes Protein')
    xlabel('Time (hours)')
    colormap(cmap)
    h=colorbar;
    ylim([0.95*min(P_plot) 1.05*max(P_plot)])
    set(gca,'FontSize',12)
    legend('Differentiation threshold')
    
    h=colorbar;
    set(h,'XTickLabel',{'Zero probability of diff', 'Higher probability'}); 
    set(h,'XTick',[0 1]); 
    
   
    
    
    
    NumberDiff=sum(DiffYN);

    %Estimate half life
    NumberDiffSub=NumberDiff-cells/2;
    [~,I]=min(abs(NumberDiffSub));
    t_half=T(I);

    t_half=t_half(1)-DiffTime; %Half life in minutes
    lam=log(2)/t_half;

    Td=T(DiffTime/dt+1:end);

    D=cells.*(1-exp(-lam.*(Td-DiffTime)));

%     figure(13324)
%     clf
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
% 
%     plot(T./60, linspace(cells,cells,Nt+1),'--','color',[0.6 0.6 0.6]); hold on
%     plot(Td./60, D, 'color', [0.7 0.3 0.3],'linewidth', 5)
%     plot(T./60, NumberDiff,'color','w','linewidth',1.5)
%     title(sprintf('Fitted half life = %.1f hours', t_half./60))
%     ylabel('Number of differentiated cells')
%     xlabel('Time (hours)')
%     legend('Total # cells', 'Exponential fit to data', 'Number of differentiated cells','location','southeast')
%     ylim([0 cells*1.1])

    DiffSignal=DiffYNflash(:,DiffTime/dt:end);
    WindowSize_Minutes=150;
    WindowSize=WindowSize_Minutes/dt;

    figure(13325)
    clf
    fig=gcf;
        if GraphAppearance==1
            fig.InvertHardcopy = 'off';
        else
            fig.InvertHardcopy = 'on';
        end

    h=plot(T(DiffTime/dt:end)./60,60*movmean(sum(DiffSignal),WindowSize)./dt); hold on
    xlim([DiffTime_hours tf_hours])
    ylim([0 70])
    ylabel('Differentiation rate diffs/hour')
    xlabel('Time (hours)')
    title(sprintf('Differentiation rate with sliding window size of %.f minutes', WindowSize_Minutes))
    set(h, {'color'}, num2cell(randColourMapBright(1,[0 0 0]),2));
    drawnow
    
    fprintf('Differentiation rate over whole time (%% of cell population) = %.1f %%/hour', 100*sum(DiffSignal(:))/((tf_hours-DiffTime_hours)*cells));
end


%% Final time point of hexagonal grid plot
if ShowLastFrame==1
    figure(104)
    clf;
    
    fig=gcf;
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
%     colour_index=reshape(flipud(vecTOmat(P(:,t_step),cols)),[1,cols*rows])
%     colour_index=reshape(flipud(vecTOmat(TauH(:),cols)),[1,cols*rows]);
    colour_index=flipud(grad(:,1));
%     map=viridis(1000);
    
    size(X_vertices)
    size(Y_vertices)
    size(colour_index)
%     hexagons = patch(X_vertices,Y_vertices,colour_index,'edgecolor','none');
    hexagons = patch(X_vertices,Y_vertices,colour_index);
    colormap(gray); 
    colorbar
    axis equal;
    set(gca,'Visible','off')
%     caxis([min(min(P(:,round(0.2*Nt):end))) max(max(P(:,round(0.2*Nt):end)))])
%     colour_index=reshape(flipud(vecTOmat(P(:,end),cols)),[1,cols*rows])';
    colour_index=flipud(grad(:,1));
%     reshape(flipud(vecTOmat(TauH(:),cols)),[1,cols*rows]);
    set(hexagons, 'FaceVertexCData',colour_index);
    title(sprintf('Time: %.0f hours', tf/60));
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'FontSize',13)

    drawnow
end

if ShowLastFrame==1
    figure(103)
    clf;
    
    fig=gcf;
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    n=1; %Length of hexagon side
    [X_vertices, Y_vertices]=hex(rows,cols,n); % Returns hexagonal grid vertices
    colour_index=reshape(flipud(vecTOmat(P(:,end),cols)),[1,cols*rows]);
%     colour_index=reshape(flipud(vecTOmat(TauH(:),cols)),[1,cols*rows]);
%     colour_index=flipud(grad(:,1));
%     map=viridis(1000);
    
    size(X_vertices)
    size(Y_vertices)
    size(colour_index)
    hexagons = patch(X_vertices,Y_vertices,colour_index,'edgecolor','none');
%     hexagons = patch(X_vertices,Y_vertices,colour_index);
    colormap(inferno(300)); 
    colorbar
    axis equal;
    set(gca,'Visible','off')
%     caxis([min(min(P(:,round(0.2*Nt):end))) max(max(P(:,round(0.2*Nt):end)))])
    colour_index=reshape(flipud(vecTOmat(P(:,end),cols)),[1,cols*rows])';
%     colour_index=flipud(TauH);
%     reshape(flipud(vecTOmat(TauH(:),cols)),[1,cols*rows]);
    set(hexagons, 'FaceVertexCData',colour_index);
    title(sprintf('Time: %.0f hours', tf/60));
    xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
    ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'FontSize',20)

    drawnow
end

%% Plots of Hes levels in random cells 
if ShowRandomCells==1
    
    % Figure preamble
    figure(2) 
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
    else
        fig.InvertHardcopy = 'on';
    end
    
    
    %How many cells to plot
    desired_cells=30;
    desired_fraction=desired_cells/cells;
    rand_ind=ceil(cells*rand(1,ceil(cells*desired_fraction)));
    
    if cells==2
        rand_ind=[1 2];
    end
    
    num=length(rand_ind);
    
    h=plot(T/60,P(rand_ind,:),'linewidth',1); hold on
    plot(T((DiffTime+1)/dt:end)/60,DiffThresh(1,(DiffTime+1)/dt:end),'k--','linewidth',1.5); hold on
%     set(h, {'color'}, num2cell(randColourMapBright(num,[0 0.5 0]),2));
    set(h, {'color'}, num2cell(randColourMap(num),2));
    xlabel('Time (hours)')
    ylabel(sprintf('Hes protein count'))
    set(gca,'FontSize',15)
    title('Subset of cell Hes5 time traces')
    ylim([0 5e4])

% DiffP=P(:,DiffTime/dt:end);
% SD=std(DiffP(:))
% MEAN=mean(DiffP(:))
% CoV=SD/MEAN
    
    
%Seperate trajectories on subplots
    numPlots=10;
    ymin=min(min(P(:,Nt/2:Nt)));
    ymax=max(max(P(:,Nt/2:Nt)));
    if numPlots>cells
        numPlots=cells;
    end
    rand_ind=ceil(cells*rand(1,numPlots));
%     rand_ind=randperm(size(Y_raw,1));
    t_start=180*60/dt;
    figure(8)
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    set(fig,'renderer','Painters')
    col = [0 0 0 0.2];
    
    for n=1:numPlots
        
%         subplot(4,6,n) 
        if n<numPlots
            rand_ind(n)
            h11 = plot(T(t_start:end)/60,P(rand_ind(n),t_start:end),'linewidth',1,'color', col); hold on
        else
            h11 = plot(T(t_start:end)/60,P(rand_ind(n),t_start:end),'linewidth',3,'Color', 0.8.*[103, 130, 103]./255 ); %[103  122  103]./255
        end

        %         plot(T/60,DiffThresh(1,:),'w--','linewidth',1.5)
%         set(h11, {'color'}, num2cell(randColourMapBright(1,[0.2 0.5 0.2]),2));
        xlabel('Time (hours)')
        ylabel(sprintf('HES5 abundance'))
        set(gca,'FontSize',14)
        xlim([t_start/60,tf_hours]);
%         ylim([ymin ymax]);
        ylim([0 70000]);
    end
    drawnow
    
end

%==========================================================================
%%                          Fourier transform
%==========================================================================
    
  
    %Changes fraction of signal to use (deterministic uses earlier part of
    %the signal as it is a damped oscillator)
    if Stochastic==1
        frac=0.2;
    else   
        frac=0.01;
    end
    Osc=find(max(P(:,0.7*Nt:Nt),[],2)>1000);
    Y_raw=P(Osc,:);
     
    %Detrending
    polyOrder=2;                  %Order of detrending polynomial in detrend.m function
    frameTime=40;                 %Frame length in hours for detrending window
    frameLength=frameTime*60./dt; %Conversion to window length in elements  
%     frameLength=75;
    [Ydetrend,t,Ysmooth,f,P1,coherence,f_C1,P_C1,avgFourier,ind_per,I]=tempFourier(T,Y_raw,polyOrder,frac,frameLength,Nt,dt);
    
if TemporalFourier==1   
    
    fprintf('\nExpected coherence from single cell posterior data = %.2f\n', SS(s,4))    
    fprintf('Coherence of oscillators in this model = %.2f\n',coherence)
    fprintf('\nExpected peak period from single cell summary stats = %.2f\n', summary_stats(s,3)/60)
    fprintf('Oscillators with periods below 7 hours = %.2f%%\n', 100*sum(ind_per<7)/numel(ind_per))
%__________________________________________________________________________    
%__________________________________Plots___________________________________

    figure(14)
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    fig.InvertHardcopy = 'off';
    txtSize=10;
    sampleCell=size(Y_raw,1);
    
    subplot(321)
    plot(t./60,Ysmooth(sampleCell,:),'--', 'color','k'); hold on
    plot(T(frac*Nt:end)./60,Y_raw(sampleCell,frac*Nt:end),'linewidth',2,'color', [0 0.6 0.6])
    title('Signal with detrend line')
    xlim([dt*Nt*frac/60 dt*Nt/60])
    xlabel('Time (hours)')
    set(gca,'Fontsize',txtSize)
    
    subplot(322)
    plot([t(1)/60 t(end)/60], [0 0],'w'); hold on
    plot(t./60,Ydetrend(sampleCell,:),'color',[0 0.6 0.6],'linewidth',2); 
    xlim([t(1)./60 t(end)./60])
    title('Detrended signal')
    xlabel('Time (hours)')
    set(gca,'Fontsize',txtSize)
    
    subplot(3,2,3)
    h0=plot(60*f,P1);
    xlim([0 0.5])
    xlabel('Frequency (1/h)')
    set(h0, {'color'}, num2cell(randColourMapBright(size(P1,1),[1 0 0]),2));
    ylim([0 1.2*max(max(P1))])
    title('Individual cell Fourier transforms')
    set(gca,'Fontsize',txtSize)
    
    subplot(324)
    histo=histogram(ind_per,'normalization','probability');
%     if CoupledCells == 0
%         histo=histogram(ind_per,'normalization','probability');
%     end
    set(histo,'edgecolor','w');
    set(histo,'facecolor',[0.4 0.4 0.4]);
    xlabel('Dominant period (h)')
    ylabel('Num of cells')
    title('Individual cell dominant periods')
    set(gca,'fontsize',txtSize)
    xlim([2 20])
    
    subplot(3,2,[5 6])
    h1=area(60*f,avgFourier,'LineStyle','none'); hold on;
    h1.FaceColor=[0 0.6 0.6];
    h2=area(60*f_C1,P_C1,'LineStyle','none'); hold on;
    h2.FaceColor=[0.9 0.4 0.4];
    plot(60*f,avgFourier,'.','linewidth',1,'color','w'); 
%     plot(60*f(I-ten_percent),avgFourier(I-ten_percent),'ro','linewidth',2)
%     plot(60*f(I+ten_percent),avgFourier(I+ten_percent),'ro','linewidth',2)
    plot(60*f(I),avgFourier(I),'.','markersize',20,'color', [1 1 1])
    title('All cells Fourier transform')
    xlabel('Frequency (1/h)')
    xlim([0 0.5])
    ylim([0 1.2*max(avgFourier)])
    title(sprintf('Sum of Fourier Transforms   |   Dominant period: %.2f hours   |   Coherence: %.2f', 1./(60*f(I)), coherence))
    set(gca,'Fontsize',txtSize)
    drawnow
    
    
%__________________________________________________________________________
%Light mode report version
    figure(15)
   set(gcf,'renderer','Painters')
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    txtSize=13;
    sampleCell=size(Y_raw,1);
    
%     subplot(1,2,1)
    h1=area(60*f,avgFourier,'LineStyle','none'); hold on;
    h1.FaceColor=[103  130  103]./255;
    h2=area(60*f_C1,P_C1,'LineStyle','none'); hold on;
    h2.FaceColor=[205 227 210]./255; %[205 227 210]./255
    plot(60*f,avgFourier,'.','linewidth',1,'color',0.2*[1 1 1]); 
%     plot(60*f(I-ten_percent),avgFourier(I-ten_percent),'ro','linewidth',2)
%     plot(60*f(I+ten_percent),avgFourier(I+ten_percent),'ro','linewidth',2)
    plot(60*f(I),avgFourier(I),'.','markersize',20,'color', 0.2*[1 1 1])
%     title('All cells Fourier transform')
    xlabel('Frequency (1/h)')
    xlim([0 0.5])
    ylim([0 1.2*max(avgFourier)])
    title({'Average power spectrum', sprintf('Dominant period: %.2fh   |   Coherence: %.2f', 1./(60*f(I)), coherence)})
    title({sprintf('Dominant period: %.2f hours   |   Coherence: %.2f', 1./(60*f(I)), coherence)})
    
    set(gca,'Fontsize',txtSize)
    drawnow
    
%__________________________________________________________________________
    

    numPlots=24;
    if numPlots>cells
        numPlots=cells;
    end
    
    rand_ind=randperm(size(Y_raw,1));
    rand_ind=rand_ind(1:numPlots);
    
%     t_start=80*60/dt; %fractional start time in elements
    t_start=Nt*frac; 
    t_start_h=dt*Nt*frac/60; %fractional start time in hours

    figure(9)
    clf
    fig=gcf;
    fig.InvertHardcopy = 'off';
    
    for n=1:numPlots
      
        subplot(4,6,n) 
        plot(t./60, Ysmooth(rand_ind(n),:),'w--'); hold on
        h12=plot(T(t_start:end)/60,Y_raw(rand_ind(n),t_start:end),'linewidth',1);
        set(h12, {'color'}, num2cell(randColourMapBright(1,[0.2 0.5 0.2]),2));
        xlabel('Time (hours)')
        ylabel(sprintf('Hes protein'))
        set(gca,'FontSize',8)
        xlim([t_start_h,tf_hours])
        title(sprintf('T=%.1fh',ind_per(rand_ind(n))))
        
    end
    drawnow
    
end

%==========================================================================
%%                        Synchronisation tests
%==========================================================================

if PerformMFV==1
    
    [PH,PH_unwrap,COP]=phase(t,Ydetrend); % Uses the Hilbert transform to extract phase from oscillatory signals
    Sync_val=mean(abs(COP(round(0.5*length(COP)):round(0.8*length(COP)))));
    fprintf(sprintf('\nMean field/Kuramoto order value = %.2f \n', Sync_val))
    
    mi=size(PH_unwrap,1);
    x_poly    = linspace(t0,tf,Nt+1);
    PH_smooth=zeros(mi,length(x_poly));
    
    for m=1:mi
        [p,S,mu]  = polyfit(t,PH_unwrap(m,:),8);
        PH_smooth(m,:) = polyval(p,x_poly,S,mu);
    end

    elem=1;
    PH1=PH_smooth(:,1:end-elem);
    PH2=PH_smooth(:,1+elem:end);
    inst_freq=((PH2-PH1)./(2*pi))./(elem*dt); %cycles/min

    % Filter out non-oscillating cells
    % CutOffPer=20; % Hours
    % CutOffFreq=1/(60*CutOffPer); %cycles/min
    % Osc=find(min(inst_freq,[],2)>CutOffFreq);
   
    
end

%% Scatter plot of phase
if AnimateComplexPhase==1
%     [PH,PH_unwrap,COP]=phase(T,P); % Uses Hilbert transform to extract phase from oscillatory signals
%     Sync_val=mean(abs(COP(0.5*Nt:0.8*Nt)));

    r = ones(cells,1);

    Hilbert_freq=60*(PH_smooth(:,0.8*Nt)-PH_smooth(:,0.2*Nt))./(0.5*Nt*dt*2*pi);
    Hilbert_per=1./Hilbert_freq;
    COP=(1/cells)*sum(exp(1i*PH_smooth)); % Complex order parameter
    th2=angle(COP);

    idx=0;
    for j=1:round(AnimationSpeed*Nt/400):1*Nt
        idx=idx+1;
       
        figure(8)
%         th = wrapTo2Pi(PH_unwrap(:,j));
        th =PH(:,j);
        polarplot(th,r,'o','color','w'); hold on
        polarplot(th2(j),abs(COP(j)),'o','color',[205/255 102/255 120/255],'linewidth',6);
        hold off
%         title(sprintf('Mean field value = %.2f\nTime = %.0f hours',abs(COP(j)),T(j)./60))
        title(sprintf('Mean field value = %.2f',abs(COP(j))),'color',[205/255 102/255 120/255])
        
        rlim([0 1.5])
        rticks(1)
        rticklabels({''})
        thetaticks([0 90 180 270 ])
        thetaticklabels({'0', '\pi/2', '\pi', '3/2\pi'})
        set(gcf,'color',[0.25 0.25 0.25])
        set(gca,'Color',[0.25 0.25 0.25])
        set(gca,'thetacolor', 'w')
        set(gca,'fontsize',25) 
        drawnow;
        
        if MakeGIF==1
            F1=getframe(gcf); %Make video
            im{idx}=frame2im(F1);
%             writeVideo(vidObj,F);
        end

    end
    
        if MakeGIF==1
            IDX=idx;
            for idx = 1:IDX
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename2,'gif','LoopCount',Inf,'DelayTime',1/30);
                else
                    imwrite(A,map,filename2,'gif','WriteMode','append','DelayTime',1/30);
                end
            end
        end
    
%     for j=3:10:0.8*Nt
%     PH_combinations = nchoosek(PH(:,j), 2);
%     
%     figure(7)
%     dscat=dscatter(PH_combinations(:,1),PH_combinations(:,2));
%     get(dscat)
%     title(sprintf('Time=%.0f hours', j*dt/60))
%     xlim([-pi pi])
%     ylim([-pi pi])
%     xticks([-pi -pi/2 0 pi/2 pi])
%     xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
%     yticks([-pi -pi/2 0 pi/2 pi])
%     yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
%     axis square
%     drawnow;
%     
%     if MakeVideo==1
%             F=getframe(gcf); %Make video
%             writeVideo(vidObj,F);
%     end
%         
%     end
end


%% Cluster detection
if CrudeDiff==1
    thresh=0.92;
    frac=DiffTime/Nt;

    if frac>1
        error('Ensure DiffTime is less than tf')
    end

    [Y,t,Ysmooth,f,P1,coherence,f_C1,P_C1,avgFourier,ind_per,I]=tempFourier(T,Y_raw,polyOrder,frac,frameLength,Nt,dt);
    [PH,~,~]=phase(t,Y);   % Returns wrapped phase of all cells and time points
    dtWinHr=8;            % Time window in hours in which to perform pairwise KOP
    dtWin=dtWinHr*60/dt;   % Time window in elements

    tI=floor(length(t)/dtWin);
    COLN=cols;
    % tI=1;
    % COLN=1;

    TWin=0:dtWin:tI*dtWin;
    TWinHr=TWin*dt./60;
    firstEntry=1;

    coolElHrs=16;
    coolEl=coolElHrs*60/dt;
    DiffYNflash=cooldown(DiffYNflash,coolEl);

    for ti=1:tI
        for colN=1:COLN

            winL=TWin(ti)+1;
            winR=TWin(ti+1);
            tL=winL+t(1);
            tR=winR+t(1);

            %Extract 1 column
        %     colN=2;
            ni=(colN-1)*rows+1;
            nf=colN*rows;
            PHcol=PH(ni:nf,winL:winR);
            Ycol=Y(ni:nf,winL:winR);

            % Get average KOP over time comparison matrix
            kop=KOPcomp(PHcol);
            kopFilt=kop; kopFilt(kop<thresh)=0; kopFilt(kop>thresh)=1;
            [kopNeigh]=clusterDet(kopFilt);
            % kopNeigh(eye(rows)==1)=0; %remove diagonal


            if COLN==1
                figure(132)
                clf
                subplot(131)
                imagesc(kop)
                axis square
                title('Synchronisation comparison matrix')
                ylabel('Cell i')
                xlabel('Cell j')
                colorbar
                subplot(132)
                imagesc(kopFilt)
                axis square
                colorbar
                title(sprintf('Logical comparison matrix for KOP values above %.1f',thresh))
                ylabel('Cell i')
                xlabel('Cell j')
                subplot(133)
                imagesc(kopNeigh)
                axis square
                colorbar
                title('Only showing neighbouring cells with high KOP')
                ylabel('Cell i')
                xlabel('Cell j')   
            end

            clustSize=sum(kopNeigh,2);
            clustIdx=find(clustSize>1)+ni-1; %Location of 'centre' cell in cluster
            clustSize=clustSize(clustSize>1);

            for n=1:length(clustIdx)
                cellsInClust=find(kopNeigh(clustIdx(n)-ni+1,:)==1)';
                diffEvents=sum(DiffYNflash(cellsInClust+ni-1,tL:tR),2);

                absDiff=sum(diffEvents>0);
                percentDiff=100*sum(diffEvents>0)/clustSize(n);
        %         clustDiffEvents(n)=sum(DiffYNflash(kopNeigh(clustCells(n)),tL:tR));

                if firstEntry==1
                    AD=absDiff;
                    CS=clustSize(n);
                    PD=percentDiff;
                else
                    AD=[AD; absDiff];
                    CS=[CS; clustSize(n)];
                    PD=[PD; percentDiff];
                end

                firstEntry=0;

            end

            if COLN==1
                figure(133)
                clf
                for i=1:rows
                    subplot(5,6,i)
                    yPlot=Ycol(kopNeigh(i,:)==1,:);
                    plot(t(winL:winR)./60,yPlot);hold on

                    diffExtract=DiffYNflash(find(kopNeigh(i,:)==1)+ni-1,tL:tR);
                    if sum(diffExtract(:)>0)
                        yPlotDiff=yPlot;
                        yPlotDiff(diffExtract==0)=nan;

                        tDiff=t(winL:winR);
                        tDiff=repmat(tDiff,[size(yPlotDiff,1) 1]);
                        tDiff(diffExtract==0)=nan;

                        plot(tDiff./60,yPlotDiff,'kx')
                    end

                end
            end


        end
    end


    figure(171)
    clf
    histogram(CS)
    title('Distributon of cluster sizes detected')
    ylabel('counts')
    xlabel('Cluster size')
    ylim([0 900])
    xlim([1.1 10.9])

    clustSizePlot=3;

    figure(172)
    subplot(212)
    histogram(PD(CS==clustSizePlot),20)
    title(sprintf('Percentage of cluster differentiating in cluster size = %.f.',clustSizePlot))
    xlabel('Percent differentiating in cluster')
    ylabel('Count')
    subplot(211)
    histogram(AD(CS==clustSizePlot),20)
    title(sprintf('Number of cells differentiating in cluster size = %.f.',clustSizePlot))
    xlabel('Number of cells differentiating in cluster')
    ylabel('Count')

    figure(174)
    clf
    histogram(AD(CS==clustSizePlot & AD>0),20)
    title(sprintf('Number of cells differentiating in cluster size = %.f.',clustSizePlot))
    xlabel('Number of cells differentiating in cluster')
    ylabel('Count')
end


%==========================================================================
%%                          Wavelet Transform
%==========================================================================
% 
if TemporalWavelet==1
    
% wc value is the central frequency, typical values range from 0.5 to 4: 
%                  wc  =  0.5   1   2   3   4    
% Better time resolution <-------------------> Better frequency resolution

    wc=4; %Central frequency
    
    Ntime=length(T)/100; %Reduce time resolution to speed up computation
    Nw=100;
   
    w0h=0.01;       %Lower frequency bound in 1/hours
    wfh=1;          %Upper frequency bound in 1/hours
    w0h_c=w0h/wc;   %Correcting for wc scaling
    wfh_c=wfh/wc;   %Correcting for wc scaling
    w0=w0h_c/3600;  %lower frequency bound in Hz
    wf=wfh_c/3600;  %final frequency bound in Hz
%     dw=(wf-w0)/Nw;

    ts=T*60; %Time in seconds
    
    fh=0.8;
    fs=fh/3600;
%     x=sin(2*pi*fs*ts);
%     x=P(2,:)-mean(P(2,:));
    x=P(1,:);

    % Need to subtract mean here!
    
    [WT,W,t]=wavelet(x,Nw,w0,wf,wc,Ntime,ts);

    Wh=W*3600; %Frequency in 1/hours
    Th=t/3600; %Time in hours
   
    [~,I]=max(abs(WT).^2,[],1);
    domPer=1./Wh(I);
    domPer(I==1)=NaN;
    
    if wfh==1
%         fprint('banana')
        Wh(end)=1;
    end
    
    %plot
    figure(1)
    clf
    subplot(311)
    plot(T./60,x,'k')
    ylabel('Hes protein number')
    xlabel('Time (hours)')
    subplot(312)
%     [map2]=colourMap2(0.4);
    load cmap.mat
    colormap(cmap)
%     colormap(gca,map2)
    surf(Th,Wh,abs(WT).^(2),'EdgeColor', 'None', 'facecolor', 'interp');
    xlabel('Time (hours)')
    ylabel('Frequency (1/hours)')
    xlim([t0 tf_hours])
    ylim([Wh(1) Wh(end)])
    view(2);
    set(gca,'yscale','log')
    
    subplot(313)
    plot(Th,domPer,'k')
    xlabel('Time (hours)')
    ylabel('Dominant period (hours)')
    xlim([t0 tf_hours])
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Spatial Frequency Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if SpatialFourier==1 && cols==1
    startTime=tf_hours-100;
    time_frac=startTime/tf_hours;
    
    if time_frac>=1
        error('Increase length of simulation or decrease startTime for spatial frequency analysis.')
    end
    
    t=T(time_frac*Nt:end)./60;
    ti=time_frac*Nt;

    [Pnew]=improvedKymo(P,rows,cols,kymoCols); %Improved selection of kymograph
    
    kymoWidth=1; %Number of cells to use for kymograph selection width
    K=cols/kymoWidth;
    
    t_concat=t;
    Yraw=Pnew(:,ti:end);

    map=viridis(500);
    X=t(1:10:end);
    Y=1:rows;
    plotData=Yraw(:,1:10:end);
    Y_RAW=Yraw;
    
    if CrudeDiff==1
    
        figure(301)
        clf
        sgtitle('Kymograph')
        subplot(1,4,[1 2 3 ])
        j=imagesc(X,Y,plotData);
        xlim([X(1) tf_hours])
        ylabel('Row index')
        xlabel('Time (h)')
        set(gca,'YTickLabel',[]);
        set(gca,'FontSize',10)
        colormap(gca,map)
        colorbar('Location','westoutside');
        view(0,90)
        set(gca,'fontsize',10)
        
        
        subplot(1,4,4)
        histData=sum(DiffYNflash,2);
        histData= flipud(histData);
        bar(histData,1,'facecolor',1/256.*[80 80 80],'edgecolor','w')
        set(gca,'view',[90 -90])
        xlim([0.5 rows+0.5])
        set(gca,'XTickLabel',[]);
        set(gca,'fontsize',10)
        ylabel('Diff event count')
        
        diffSignal=DiffYNflash(:,DiffTime/dt:end);
        indivDiffFreq = sum(diffSignal,2)/(tf_hours-DiffTime_hours);
        
    end
    

    %%

    %Split Y_RAW into x hour windows of expression and then average
    t_elem=(1-time_frac)*Nt; %-1 is for the fact the first column is not included in this analysis
    split_time=0.1;          %in hours (Nominal value of 2 hours)
    if split_time==0
        split_elem=1;
    else
        split_elem=round(split_time*60/dt);
    end
    
    
    M=1:split_elem:t_elem;
    M=round(M);
    II=length(M);
    
    split_t=linspace(0,tf_hours*K,II);
    
    Y_split=Y_RAW(:,M);
  
    t_split=t_concat(M);
    
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1;        %End positions marked by -1
    clipPos(clipPos>0)=0;         %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    clipEnd=find(clipPos==-1);
    
    YDETREND2=Y_split-repmat(mean(Y_split,1),[rows 1]); %Detrending in spatial direction by subrtracing poulation mean at each time point
    
    Y=fft(YDETREND2,[],1);
    winL=length(Y(:,1));
    P2 = abs(Y/winL).^2;
    P1 = P2(1:winL/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:); 

    Fs=1;
    f = Fs*(0:(winL/2))/winL;
 
    KymLong=YDETREND2;
    
    [FisherPer,Occurrence]=Fisher(P1, f, rows);
  
  
%__________________________________________________________________________

    figure(61) 
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    fontSize=10;
    
    if GraphAppearance==1
        fig.InvertHardcopy = 'off';
        
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
        xlabel('Column over time')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'w')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
        ylabel('Significant period (cells)')
        xlabel('Column over time')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
%         ylim([0 rows])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,rows],'k')
        end
        ylim([0 rows])

        drawnow

    else

        fig.InvertHardcopy = 'on';
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
        xlabel('Column')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'k')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
        ylabel('Significant period (cells)')
        xlabel('Column')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
%         ylim([0 5])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,rows],'k')
        end
        ylim([0 rows])
        drawnow
    end 
    
    %% Subplots of all columns/kymos
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1; %End positions marked by -1
    clipPos(clipPos>0)=0;  %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    
    clipStart=find(clipPos==1);
    clipEnd=find(clipPos==-1);
    
    plotCols=3;
    plotRows=cols/plotCols;
    plotCols=plotCols*2;
    
    thresh=0.6;
    kymo=KymLong;
    l=zeros(size(kymo,2),1);
    for i=1:size(kymo,2)
        vect=kymo(:,i);
        [l(i),p(i)]=microClustDetect(vect,thresh);
    end
    clusterMean=nanmean(l);
    clusterOcc=1-sum(isnan(l))/numel(l);
        
end

vid=viridis(500);
if SpatialFourier==1 && cols>1
    startTime=tf_hours-100;
    time_frac=startTime/tf_hours;
    
    if time_frac>=1
        error('Increase length of simulation or decrease startTime for spatial frequency analysis.')
    end
    
    t=T(time_frac*Nt:end)./60;
    ti=time_frac*Nt;
    [Pnew]=improvedKymo(P,rows,cols,kymoCols); %Improved selection of kymograph
    
    kymoWidth=1; %Number of cells to use for kymograph selection width
    K=cols/kymoWidth;
    
    
    if mod(cols,2)==1
        error('Make cols even')
    end
    
    figure(23)
    clf 
    t_concat=t;
    
    

    for k=2:K %Starts at k=2 because the first column in not included in the analysis (to do with the improvedKymo function)

        figure(23)

        Yraw(:,:,k)=Pnew((k-1)*rows+1:k*rows,ti:end);
        
        if k==2
            Y_RAW=Yraw(:,:,2);
            t_concat=t;
        else
            Y_RAW=[Y_RAW Yraw(:,:,k)];
            t_concat=[t_concat t];
        end    
            
        subplot(ceil(K/3),3,k-1)
        
        X=t;
        Y=1:rows;
        j=surf(X,Y,flipud(squeeze(Yraw(:,:,k))));
        title(sprintf('Kymograph %.f',k))
        ylim([Y(1) Y(end)])
        xlim([t(1) t(end)])
        xlabel('Time (hrs)')
        ylabel('Row index')
        set(j, 'LineStyle','none')
        set(gca,'YTickLabel',[]);
        set(gca,'FontSize',10)
        colormap(gca,inferno(300))
        colorbar
        view(0,90)
        axis square

        if k == 3
            figure(24)
            clf
            
            set(gcf,'renderer','painters')
            fig=gcf;

            X=t;
            Y=1:rows;
            j=surf(X,Y,flipud(squeeze(Yraw(:,:,k))));
            title(sprintf('Kymograph %.f',k))
            ylim([Y(1) Y(end)])
            xlim([t(end)-20 t(end)])
            xlabel('Time (hrs)')
            ylabel('Row index')
            set(j, 'LineStyle','none')
            set(gca,'YTickLabel',[]);
            set(gca,'FontSize',10)
            colormap(gca,inferno(300))
            colorbar
            view(0,90)
            axis square
            set(gca, 'fontsize', 13)
        end

    end
   
    %Split Y_RAW into x hour windows of expression and then average
    t_elem=(1-time_frac)*Nt*(K-1); %-1 is for the fact the first column is not included in this analysis
    split_time=0.1; %in hours (Nominal value of 2 hours)
    if split_time==0
        split_elem=1;
    else
        split_elem=round(split_time*60/dt);
    end
    
    
    M=1:split_elem:t_elem;
    M=round(M);
    II=length(M);
    
    split_t=linspace(0,tf_hours*K,II);
    
    Y_split=Y_RAW(:,M);
  
    t_split=t_concat(M);
    
    clipPos=t_split(2:end)-t_split(1:end-1);
    clipPos(clipPos<0)=-1; %End positions marked by -1
    clipPos(clipPos>0)=0;  %Fill rest with zeros
    startPos=find(clipPos==-1)+1; %Find start positions of clips
    clipPos(startPos)=1;
    clipPos=[clipPos, -1];
    clipPos(1)=1;
    clipEnd=find(clipPos==-1);
    
    YDETREND2=Y_split-repmat(mean(Y_split,1),[rows 1]); %Detrending in spatial direction by subrtracing poulation mean at each time point
    
    
    Y=fft(YDETREND2,[],1);
    winL=length(Y(:,1));
    P2 = abs(Y/winL).^2;
    P1 = P2(1:winL/2+1,:);
    P1(2:end-1,:) = 2*P1(2:end-1,:);

    Fs=1;
    f = Fs*(0:(winL/2))/winL;
 
    KymLong=YDETREND2;
    
    [FisherPer,Occurrence]=Fisher(P1, f, rows);

   

%__________________________________________________________________________

    figure(61) 
    clf
    set(gcf,'renderer','Painters')
    fig=gcf;
    fontSize=10;
    
    if GraphAppearance==1

        fig.InvertHardcopy = 'off';
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
        xlabel('Column over time')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'w')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on
        title(sprintf('Fisher G significant period detected %.0f%% of the time', 100*Occurrence))
        ylabel('Significant period (cells)')
        xlabel('Column over time')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])
%         ylim([0 10])

        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,rows],'w')
        end
        ylim([0 rows])

        drawnow

    else
        fig.InvertHardcopy = 'on';
        subplot(2, 6, [1 2 3 4 5])
        load cmap
        colormap(gca,cmap)
        P1=[P1; P1(1,:)];
        f=[f, (2*f(end)-f(end-1))];

        t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));

        P1_lines=P1;
        P1_lines(:,clipEnd)=NaN;

        h=surf(t_plot,f,P1_lines);
        xlabel('Column')
        ylabel('Frequency (1/cell)')
        title('Periodogram over time for individual columns of cells')
        xlim([ t_plot(1) t_plot(end)])
        ylim([f(1) 0.5+f(end)-f(end-1)])
        view(0,90)
        colorbar
        set(h,'LineStyle','none')
        set(gca,'FontSize',fontSize)
        pbaspect([2.3 1 1])
        colormap(gca,cmap)

        centre=(t_plot(end)-t_plot(1))./((cols-1)*2); %cols -1 for improved kymo/removed first column
        xticks(linspace(t_plot(1)+centre,t_plot(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})

        subplot(266)
        [~,maxIdx]=max(mean(P1,2));
        plot(mean(P1,2),f,'k')
        title({'Mean FT. ',sprintf('Peak period = %.1f', mean(1./f(maxIdx)))})
        ylabel('Frequency (1/cells)')
        xlabel('Power')
        pbaspect([0.5 1.4 1])
        set(gca,'FontSize',fontSize)
        ylim([0 0.5])

        subplot(2,6,[7 8 9 10 11 12 ])
        plot(split_t,FisherPer,'s','markerSize',3,'MarkerEdgeColor','none', 'MarkerFaceColor',cmap(size(cmap,1)/2,:));hold on

        title(sprintf('Fisher G significant period detected %.0f%% of the time. Mean period = %.1f', 100*Occurrence, nanmean(FisherPer)))
        ylabel('Significant period (cells)')
        xlabel('Column')
        set(gca,'FontSize',fontSize)
        xlim([split_t(1) split_t(end)])


        centre=(split_t(end)-split_t(1))./((cols-1)*2);
        xticks(linspace(split_t(1)+centre,split_t(end)-centre,(cols-1)))
        xticklabels({num2str((2:(cols))')})
        h=gca; h.XAxis.TickLength = [0 0];

        for n=1:numel(clipEnd)
            plot([split_t(clipEnd(n)), split_t(clipEnd(n))],[0,rows],'k')
        end
                ylim([0 rows])

        drawnow
    end
    
end
 
%%
[counts,edges] = histcounts(FisherPer);
percentage = 100.* counts ./length(FisherPer);

figure(62)
clf
histogram('BinEdges',edges,'BinCounts',percentage, 'facecolor', 0.4*[1 1 1], 'facealpha', 1)
title(sprintf('Mean period = %.1f cells', nanmean(FisherPer)))
ylabel('% of time')
xlabel('Significant spatial period (cells)')
axis square
set(gca, 'fontsize', 13)

ylim([0 100])
xlim([4 16])


%%


% %==========================================================================
% %%                    Spatial frequency over time
% %==========================================================================
% vid=viridis(500);
% if SpatialFourier==1
%     startTime=tf*dt/60-20;
% %     startTime = 1;
%     time_frac=startTime/tf_hours;
%     
%     if time_frac>=1
%         error('Increase length of simulation or decrease startTime for spatial frequency analysis.')
%     end
%     
%     t=T(time_frac*Nt:end)./60;
%     ti=time_frac*Nt;
%     
%     Y_RAW=(P((cols/2-1)*rows+1:cols/2*rows,ti:end)+P(cols/2*rows+1:(cols/2+1)*rows,ti:end))./2;
%     
%     kymoWidth=2; %Number of cells to use for kymograph selection width
%     K=cols/kymoWidth;
%     
%     
%     if mod(cols,2)==1
%         error('Make cols even')
%     end
%     
%     figure(23)
%     clf
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     for k=1:K
%         
%         Yraw(:,:,k)=(P(2*(k-1)*rows+1:(2*k-1)*rows,ti:end)+P((2*k-1)*rows+1:2*k*rows,ti:end))./2;
%         
%         if k==1
%             Y_RAW=Yraw(:,:,1);
%         else
%             Y_RAW=[Y_RAW Yraw(:,:,k)];
%         end
%             
%             
%         subplot(2,3,k)
%         
%         X=t;
%         Y=1:rows;
%         j=surf(X,Y,flipud(squeeze(Yraw(:,:,k))));
%         title(sprintf('Kymograph %.f',k))
%         ylim([Y(1) Y(end)])
%         xlim([t(1) t(end)])
%         xlabel('Time (hrs)')
%         ylabel('Cell index/row number')
%         set(j, 'LineStyle','none')
%         set(gca,'YTickLabel',[]);
%         set(gca,'FontSize',10)
%         colormap(gca,inferno(500))
%         colorbar
%         axis square
%         view(0,90)
%     end
%     
%     %Single Kymo
%     
%     figure(24)
%     set(gcf,'renderer','Painters')
%     clf
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     
%     
%     
% %     subplot(141)
% %     colour_index2=reshape(flipud(vecTOmat(grad(:,1),cols)),[1,cols*rows]); %Dynamic plot
% %     hex2=subplot(141);
% %     hexagons2 = patch(X_vertices,Y_vertices,colour_index2,'edgecolor','none');
% %     hexagons2 = patch(X_vertices,Y_vertices,colour_index2);
% %     set(gca,'xtick',[],'ytick',[]); 
% %     colormap(hex2,gray)
% %     colorbar; axis equal; title('Cell tracking')
% %     set(gca,'Visible','off')
% %     set(findall(gca, 'type', 'text'), 'visible', 'on')
% %     set(gca,'fontsize',10)
% %     xlim([min(min(X_vertices(:))) max(max(X_vertices(:)))])
% %     ylim([min(min(Y_vertices(:))) max(max(Y_vertices(:)))])
% 
%     
% %     subplot(1,4,[2 3 4])
%     X=t;
%     Y=1:rows;
%     plotData = flipud(squeeze(Yraw(:,:,1)));
%     if size(plotData, 2) > 6000
%         numOfSkippedElements = 10;
%         plotData = plotData(:, 1:numOfSkippedElements:end); %Crop data to reduce plotting time
%         X = X(:, 1:numOfSkippedElements:end);
%     end
%     j=surf(X,Y,plotData);
% %     title(sprintf('Kymograph %.f',k))
%     ylim([Y(1) Y(end)])
%     xlim([t(1) t(end)])
%     xlabel('Time (hrs)')
%     ylabel('Row number')
%     set(j, 'LineStyle','none')
%     set(gca,'YTickLabel',[]);
%     set(gca,'FontSize',12)
%     colormap(gca,inferno(500))
%     colorbar
%     view(0,90)
%     axis square
%     drawnow
%     
%     
%     %Split Y_RAW into x hour windows of expression and then average
%     t_elem=(1-time_frac)*Nt*K;
%     split_time=0;
%     if split_time==0
%         split_elem=1;
%     else
%         split_elem=round(split_time*60/dt);
%     end
%     
%     M=1:split_elem:t_elem;
%     M=round(M);
%     II=length(M);
%     
%     split_t=linspace(0,tf*K,II);
%     
%     Y_split=Y_RAW(:,M);
%     YDETREND2=Y_split-repmat(mean(Y_split,1),[rows 1]); %Detrending in spatial direction by subtracting poulation mean at each time point
%   
%     
% 
%     Y=fft(YDETREND2,[],1);
%     winL=length(Y(:,1));
%     P2 = abs(Y/winL).^2;
%     P1 = P2(1:winL/2+1,:);
%     P1(2:end-1,:) = 2*P1(2:end-1,:);
% 
%     Fs=1;
%     f = Fs*(0:(winL/2))/winL;
%  
%     KymLong=YDETREND2;
%     zs=zscore(KymLong); %z-scored data
% %     [PXX2,F2]=periodogram(zs,[],[],1);         %Run periodogram and boostrap intervals
%     PXX2=P1(2:end,:);
%     F2=f(2:end);
%     %Preallocation for storing stats
%     fisherG=zeros(II,1);
%     pval=zeros(II,1);
%     peakLoc=zeros(II,1);
%     peakHeight=zeros(II,1);
%     FisherPer=zeros(II,1);
%     
%     f=f(f>1/rows); %Remove whole tissue scale periodicity
%     P1=P1(f>1/rows,:);
%     for ii=1:II
% 
%         pxx=PXX2(:,ii); % periodogram
% 
%         [fisherG(ii),pval(ii),idx]=GetFisherG(pxx); % Find the peak and collect Fisher G-statistic
%         peakLoc(ii)=idx;
%         peakHeight(ii)=pxx(idx);
%         if pval(ii)<0.05
%             FisherPer(ii)=1/F2(idx);
%         else
%             FisherPer(ii)=nan;
%         end
% 
%     end
%     
%     FisherPer(FisherPer==rows)=NaN;
%     Occurrence=1-sum(isnan(FisherPer))/II;
%  
%     
%     
%     figure(61) 
%     clf
%     fig=gcf;
%     fig.InvertHardcopy = 'off';
%     
%     
%     subplot(2, 3, [1 2])
%     load cmap
%     colormap(gca,cmap)
%     t_plot=linspace(time_frac*tf/60,tf/60,numel(P1(1,:)));
%     h=surf(split_t*dt./60,f,P1);
%     xlabel('Time (hours)')
%     ylabel('Frequency (1/cell)')
%     title('Spatial Fourier transform over time')
%     xlim([split_t(1)*dt./60 split_t(end)*dt./60])
%     ylim([f(1) f(end)])
%     view(0,90)
%     colorbar
%     set(h,'LineStyle','none')
%     set(gca,'FontSize',12)
%     pbaspect([2 1 1])
%     colormap(gca,cmap)
% 
%     subplot(233)
%     [~,maxIdx]=max(mean(P1,2));
%     plot(mean(P1,2),f,'k')
% %     plot(f,P1(:,1))
%     title('Mean Fourier transform')
%     title(sprintf('Mean FT. Peak period = %.1f', mean(1./f(maxIdx))))
%     ylabel('Frequency (1/cells)')
%     xlabel('Contribution/power')
%     pbaspect([0.5 1 1])
% 
%     subplot(2,3,[4 5 6])
% %     [~,II]=max(P1,[],1);
%     plot(split_t*dt./60,FisherPer,'k.','linewidth',2)
%     title(sprintf('Significant period (Fisher G) occurence: %.1f', Occurrence))
% %     title('Significant period (Fisher G)')
%     ylabel('Significant period (cells)')
%     xlabel('Time (hours)')
%     set(gca,'FontSize',12)
%     xlim([split_t(1)*dt./60 split_t(end)*dt./60])
% %     ylim([0 10])
% drawnow
%     
%     thresh=0.6;
%     kymo=KymLong;
% 
%     for i=1:size(kymo,2)
%         vect=kymo(:,i);
%         [l(i),p(i)]=microClustDetect(vect,thresh);
%     end
%     clusterMean=nanmean(l)
%     clusterOcc=1-sum(isnan(l))/numel(l)
% %     figure(13),histogram(l,'normalization','probability'),ylabel('Frequency/occurence'),xlabel('Cluster size (radius)');
% 
% end




%==========================================================================
%%                  Phase space of mRNA and protein
%==========================================================================
if AnimatePhaseSpace==1
    idx=0;
    start=Nt/2;
    cell=1;
    for t=start:round(AnimationSpeed*Nt/400):Nt
        figure(101)
        plot(P(cell,start:t),M(cell,start:t),'--','color','w');hold on
        plot(P(cell,t),M(cell,t),'o','linewidth',6,'color',[0.8 0.2 0.2]); hold off
        xlabel('Protein count','color',[64 171 153]./255)
        ylabel('mRNA count','color',[205 102 120]./255)
        xlim([min(P(1,0.01*Nt:end)) max(P(1,0.01*Nt:end))])
        ylim([min(M(1,0.01*Nt:end)) max(M(1,0.01*Nt:end))])
        title(sprintf('Time=%.0f hours',T(t)/60))
        set(gca,'fontsize',20)
        axis square
        drawnow

%         if MakeGIF==1
%             open(vidObj);   
%         end

        if MakeGIF==1
            idx=idx+1;
            F1=getframe(gcf); %Make video with title and colour bar included
            im{idx}=frame2im(F1);
%             writeVideo(vidObj,F);
        end
    end
   
    if MakeGIF==1
        IDX=idx;
        for idx = 1:IDX
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,filename3,'gif','LoopCount',Inf,'DelayTime',1/30);
            else
                imwrite(A,map,filename3,'gif','WriteMode','append','DelayTime',1/30);
            end
        end
    end    
end
