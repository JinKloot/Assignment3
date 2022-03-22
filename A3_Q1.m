%%Use Assignment 1 Code modified with electric field 

%Jinseng Vanderkloot 101031534

%% S1 Initialization of individual electron values
clc;
clear all; 
close all;

m0 = 9.10938215e-31;            % electron mass
mn = 0.26*m0;                   % Effective mass
Temp = 300;                     % Inital Temp (K)
kb = 1.3806504e-23;             % Boltzmann constant
tmn = 0.2e-12;                  % Mean time between collision 
e = 1.602e-19;                  % Charge of an Electron (Assignment 3 for force) 

%Thermal Velocity (Question 1.A) 
vt=sqrt((2*kb*Temp)/mn);        % Sim in 2D so (2*kb*Temp), 3D is (3*kb*Temp)

% Plotting Area 
wArea = 200e-9;
lArea = 100e-9;

numElec = 30000;                 %Number of simulated Electrons 
numEPlot = 20;                  %Number of plotted Electrons 
dt = (lArea*wArea)/2;             %Typically 1/100 of region size
stepsTot = 100;                 %Total amount of steps (1000 was a long simulation) 
tTot= stepsTot*dt;              %Total Simulation time 
x = zeros(1,numElec);           %Inital X matrix          
y = zeros(1,numElec);           %Inital y matrix  
vx = zeros(1,numElec);          %Inital velocity x matrix  
vy = zeros(1,numElec);          %Inital velocity y matrix
vtot = zeros(1,numElec);        %Inital velocity matrix
VEnew = zeros(1,numElec);       %add x vector with field force velocity 
colors = rand(numElec,3);       %Color assignment for each electron                   

%Probability of Scatter 
scatOn = 1;                     %Turn Scatter on (1) or off(0)
Pscatter = 1-exp(-dt/tmn);      %Scatter Equation 
tScatter = zeros(1,numElec);    %track scatter for each particle 


%S1 Electron Random Assignments
for cnt = 1:numElec
    x(cnt)=rand()*wArea;
    y(cnt)=rand()*lArea;
    angle = (2*pi*rand());
    vx(cnt)=vt*cos(angle);   % velocity * random direction   
    vy(cnt)=vt*sin(angle);   % velocity * random direction
%     vtot(cnt)= sqrt (vx(cnt)^2)+(vy(cnt)^2);
end

%% Assignment 3 Addition of constant electric field of 0.1V
VAreaX = 0.1; %0.1V across x direction 
VAreaY = 0; %added for possibility of Y field force 
Ex = VAreaX/wArea;  % electric field along the x-axis
Ey = VAreaY/lArea;  % electric field along the y-axis
aX = (e * Ex)/mn ; % Acceleration in x due to Field 
aY = (e * Ey)/mn ; % Acceleration in y due to Field 
vxForce= (aX)*dt; % Velocity in x due to Field 
vyForce= (aY)*dt; % Velocity in y due to Field 

%The value of Ex = 50000(V/m) for 0.1V acorss the area, The force applied
%to each particle is 8e-14 N

%The acceleration due to the field in the x direction is a =
%3.381973859521560e+17 m/s^2

%% S1 Main Loop 

t=0;
intCNT = 2;
eTemp(1) = Temp;
while t < tTot
    t = t + dt; 
    
    %Store old position 
    oldx = x;
    oldy = y;
    
    %Update to new position with velocity including Field Force 
    vx(1:numElec) = vx + vxForce;
    vy(1:numElec) = vy + vxForce;
    x(1:numElec) = x(1:numElec) + (vx(1:numElec).*dt);
    y(1:numElec) = y(1:numElec) + (vy(1:numElec).*dt);
    
    vtot(1:numElec)= sqrt((vx(1:numElec).^2)+(vy(1:numElec).^2));
    
    %Apply boundary conditions 
    for check = 1:numElec
         if scatOn==1
            if Pscatter > rand()
                vx(check)=(vt/sqrt(2))*randn();
                vy(check)=(vt/sqrt(2))*randn();
                tScatter(check)= 0; %If collision, time goes to 0
            else
                tScatter(check)= tScatter(check) + dt; %track time increaing while no collision
            end
        end
        
        %If bottom contact, bounce off in opposite direction
        if (y(check)<=0)
            y(check) = 0;
            vy(check) = -vy(check);
        end
        %If top contact, bounce off in opposite directio
        if (y(check)>=lArea)
            y(check) = lArea;
            vy(check) = -vy(check);
        end
        %if left side of box, come out right side 
        if(x(check)<=0)
           x(check) = x(check) + wArea;
        end
        %if right side of box, come out left side
        if(x(check)>=wArea)
          x(check) = x(check) - wArea;
        end
    end 

    %Plot Boundary and map some electrons
    for Eplot = 1:numEPlot
        subplot (3,1,1)
        %if the electron went out of sides and back on other side, do not
        %draw line
        if abs(oldx(Eplot)-x(Eplot))<(wArea/2)
            p = plot([oldx(Eplot),x(Eplot)],[oldy(Eplot),y(Eplot)]);
        end
        p.Color=colors(Eplot,:);
        axis([0,wArea,0,lArea]);
        title('Electron Model'), xlabel('Position (m)', 'FontSize', 10), ylabel('Position (m)', 'FontSize', 10);
        
        hold on;
    end 
    pause(0.01);
    
    %Average Temprature 
    subplot (3,1,2)
    Time(:,intCNT) = t;
    allT = ((vtot(:).^2).*mn)./(2*kb);
    eTemp(:,intCNT) = mean(allT);
     
    plot(Time,eTemp,"r");
    title('Averge Temp'),xlabel('Time (s)', 'FontSize', 10), ylabel('Temp (K)', 'FontSize', 10);
    hold on;
    
    %Current over time with Field 
    Ix(:,intCNT)= e*1000000000000*mean(vx)* wArea * lArea; %Current = e * n * vd * Area (Particle density is 10e15 cm^-2)
    subplot (3,1,3)
    plot(Time,Ix,"r");
    title('Current in Area'),xlabel('Time (s)', 'FontSize', 10), ylabel('Current(A)', 'FontSize', 10);
    hold on;

    intCNT = intCNT +1;
end 
figure(2)
%Electron Temprature Map (Use Surf now to get 3D plot over the area) 
    allT = ((vtot(:).^2).*mn)./(2*kb);
    xv = linspace(min(x), max(x),100);
    yv = linspace(min(y), max(y),50);
    [X,Y] = meshgrid(xv,yv);
    ETM=griddata(x,y,allT,X,Y);
    subplot(2,1,1);
    surf(X,Y,ETM),colorbar,title("Temprature Map")
    xlabel('Width(m)', 'FontSize', 10), ylabel('Length(m)', 'FontSize', 10),zlabel('Temp(K)', 'FontSize', 10);
    axis([0, wArea, 0 lArea]);

%Electron Density Map (Use Surf now to get 3D plot over the area) 
    eMapX=linspace(0, wArea, 101); %increase by 1 for surf 
    eMapY=linspace(0, lArea, 51); %increase by 1 for surf 
    EDM=histcounts2(y,x,eMapY,eMapX);
    subplot(2,1,2)
    surf(X,Y,EDM),title('Electron Density Map');
    xlabel('Width(m)', 'FontSize', 10), ylabel('Length(m)', 'FontSize', 10),zlabel('Density (#)', 'FontSize', 10);
    
    