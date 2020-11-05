% Solving two-dimensional steady conduction equation using iterative solvers 
% Given equation is : d^2T/dx^2 + d^2T/dy^2 = 0
% Left Boundary = 400 K; Right Boundary = 800 K; Top Boundary = 600 K; Bottom Boundary = 900 K;
% nx = ny
clear all; close all; clc;
 
% Given inputs
L = 1; % Length of the domain in m
B = 1; % Width of the domain in m
 
% Nodes 
nx  = 10; % Nodes in x-direction
ny  = 10; % Nodes in y-direction
 
% Defining space variables
x = linspace(0,L,nx);
y = linspace(0,B,ny);
dx  = x(2) - x(1); %  defining grid spacing in x-direction
dy = y(2) - y(1); %  defining grid spacing in y-direction
 
%Initialization
T = ones(nx,ny);
 
% Given BCs
T(:,1) = 400;
T(:,end) = 800;
T(1,:) = 600;
T(end,:) = 900;
 
% Defining corner nodes using BCs
T(1,1) = 500; % (600+400)/2
T(end,1) = 650; % (400+900)/2
T(1,end) = 700; % (600+800)/2
T(end,end) = 850; % (900+800)/2 
 
% Creating a copy of T
Told = T; % This logic helps to write the values of T from the previous iteration
k = (2/dx^2 + 2/dy^2); % let k = (2/dx^2 + 2/dy^2)
w = 1.5; % Over-relaxation factor for SOR method
 
%Initialization of error and tolerance values 
iterative_solver = input('Solver number: Jacobi = 1 Gauss-siedel = 2 SOR = 3 \n Enter solver number:');
error = 9e9;
tol = 1e-4;
 
% Jacobi method
if iterative_solver == 1
    tic; % This command is to begin the stop watch
    jacobi_iter = 1;
        % In 'while loop' solution will be converged with the error less
        % than 1e-4
        while(error > tol)
                % For spatial discretization, another 'for loop' is used.
                % In this loop, the numerical scheme derived as shown in
                % the introduction was used. 'i-loop' is for x-direction,
                % 'j-loop' is for y-direction
                for i = 2:nx-1
                    for j = 2:ny-1
                        term1 = (1/dx^2)* (Told(i-1,j) + Told(i+1,j));
                        term2 = (1/dy^2)* (Told(i,j-1) + Told(i,j+1));
                        T(i,j) = (1/k)* (term1 + term2) ;
                    end
                end
            error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
            Told = T; % This logic helps to write the values of T from the previous step
            jacobi_iter = jacobi_iter + 1; % This is the counter for number of iterations
        end
    time = toc;  % This command is to begin the stop watch and is saved in 'time'  
    figure(1)
    contourf(x,y,T,'showText','on');
    grid on
    colormap(jet)
    title_text = sprintf('2D Steady Heat Conduction Using Jacobi Method \n Iterations = %d \n Time = %f sec', jacobi_iter, time);
    title(title_text)
    set(gca,'YDIR','reverse')
    pause(0.03)
end
 
 
% Gauss-siedal method
if iterative_solver == 2
    tic; % This command is to begin the stop watch
    GS_iter =1;
        % In 'while loop' solution will be converged with the error less
        % than 1e-4
        while(error > tol)
                % For spatial discretization, another 'for loop' is used.
                % In this loop, the numerical scheme derived as shown in
                % the introduction was used. 'i-loop' is for x-direction,
                % 'j-loop' is for y-direction
                for i = 2:nx-1
                    for j = 2:ny-1
                        term1 = (1/dx^2)* (T(i-1,j) + Told(i+1,j));
                        term2 = (1/dy^2)* (T(i,j-1) + Told(i,j+1));
                        T(i,j) = (1/k)* (term1 + term2) ;
                    end
                end
            error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
            Told = T; % This logic helps to write the values of T from the previous step
            GS_iter = GS_iter + 1; % This is the counter for number of iterations
        end
   time = toc;  % This command is to begin the stop watch and is saved in 'time' 
   figure(2)
   contourf(x,y,T,'showText','on')
   grid on
   colormap(jet)    
   title_text = sprintf('2D Steady Heat Conduction Using Gauss-Seidel Method \n Iterations = %d \n Time = %f sec', GS_iter, time);
   title(title_text)
   set(gca,'YDIR','reverse')
   pause(0.03)
end
 
% SOR method
if iterative_solver == 3
    tic; % This command is to begin the stop watch
    SOR_iter =1;
        % In 'while loop' solution will be converged with the error less
        % than 1e-4
        while(error > tol)
                % For spatial discretization, another 'for loop' is used.
                % In this loop, the numerical scheme derived as shown in
                % the introduction was used. 'i-loop' is for x-direction,
                % 'j-loop' is for y-direction
            for i = 2:nx-1
                for j = 2:ny-1
                    term1 = (1/dx^2)* (T(i-1,j) + Told(i+1,j));
                    term2 = (1/dy^2)* (T(i,j-1) + Told(i,j+1));
                    T(i,j) = (1-w)*Told(i,j) + (w/k)* (term1 + term2) ;
                end
            end
        error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
        Told = T; % This logic helps to write the values of T from the previous step
        SOR_iter = SOR_iter + 1; % This is the counter for number of iterations
        end
   time = toc;  % This command is to begin the stop watch and is saved in 'time' 
   figure(3)
   contourf(x,y,T,'showText','on')
   grid on
   colormap(jet)    
   title_text = sprintf('2D Steady Heat Conduction Using SOR Method \n Iterations = %d \n Time = %f sec', SOR_iter, time);
   title(title_text)
   set(gca,'YDIR','reverse')
   pause(0.03)
end


% Solving two-dimensional unsteady conduction equation using iterative solvers for Implicit time scheme
% Given equation is : dT/dt - alpha(d^2T/dx^2 + d^2T/dy^2) = 0
% Left Boundary = 400 K; Right Boundary = 800 K; Top Boundary = 600 K; Bottom Boundary = 900 K;
% nx = ny
clear all; close all; clc;
 
% Given inputs
L = 1; % Length of the domain in meters
B = 1; % Width of the domain in meters
dt = 0.001; % Given timestep in sec
alpha = 1.4; % Diffusivity in m^2/s
 
% Number of time steps can be calculated by nt = given time/timestep size =
% t/dt
nt = 250;
 
% Nodes 
nx  = 10; % Nodes in x-direction
ny  = 10; % Nodes in y-direction
 
% Defining space variables
x = linspace(0,L,nx);
y = linspace(0,B,ny);
dx = x(2) - x(1); %  defining grid spacing in x-direction
dy = y(2) - y(1); %  defining grid spacing in y-direction
 
%Initialization
T = ones(nx,ny);
 
% Given BCs
T(:,1) = 400;
T(:,end) = 800;
T(1,:) = 600;
T(end,:) = 900;
 
% Defining corner nodes using BCs
T(1,1) = 500; % (600+400)/2
T(end,1) = 650; % (400+900)/2
T(1,end) = 700; % (600+800)/2
T(end,end) = 850; % (900+800)/2 
 
% Creating a copy of T
Told = T; % This logic helps to write the values of T from the previous iteration
Tprev = T; % This logic helps to write the values of T from the previous timestep
 
% CFL number
k1 = alpha*dt/dx^2; % let CFL_x = k1 = alpha*dt/dx^2
k2 = alpha*dt/dy^2; % let CFL_x = k1 = alpha*dt/dx^2
cfl = k1 + k2;
% This information is given in the introduction
term1 = 1/(1 + 2*k1 + 2*k2);
term2 = k1*term1;
term3 = k2*term1;
w = 1.5; % Over-relaxation factor for SOR method
 
%Initialization of error and tolerance values 
% iterative_solver = 1;
iterative_solver = input('Solver number: Jacobi = 1 Gauss-siedel = 2 SOR = 3 \n Enter solver number:');
error = 9e9;
tol = 1e-4;
 
% Jacobi method
if iterative_solver == 1
    tic; % This command is to begin the stop watch
    jacobi_iter = 1;
        % For Time-marching, 'for loop(k)' is used 
        for k = 1:nt 
            error = 9e9;
                % In 'while loop' solution will be converged with the error less
                % than 1e-4
                while(error > tol) 
                        % For spatial discretization, another 'for loop' is used.
                        % In this loop, the numerical scheme derived as shown in
                        % the introduction was used. 'i-loop' is for x-direction,
                        % 'j-loop' is for y-direction
                        for i = 2:nx-1
                            for j = 2:ny-1
                                H = (Told(i-1,j) + Told(i+1,j));
                                V = (Told(i,j-1) + Told(i,j+1));
                                T(i,j) = (term1* Tprev(i,j) + term2*H + term3*V);
                            end
                        end
                    error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
                    Told = T; % This logic helps to write the values of T from the previous iteration
                end
            Tprev = T; % This logic helps to write the values of T from the previous timestep
            jacobi_iter = jacobi_iter + 1; % This is the counter for number of iterations
         end
    time = toc;  % This command is to begin the stop watch and is saved in 'time'   
    figure(1)
    contourf(x,y,T,'showText','on');
    grid on
    colormap(jet)
    title_text = sprintf('2D Unsteady Heat Conduction Using Jacobi Method \n Implicit Scheme \n Iterations = %d \n Time = %f sec', jacobi_iter,time);
    title(title_text)
    set(gca,'YDIR','reverse')
    pause(0.03)
end
 
% Gauss-siedal method
if iterative_solver == 2
    tic; % This command is to begin the stop watch
    GS_iter =1;
        % For Time-marching, 'for loop(k)' is used
        for k = 1:nt
            error = 9e9;
                % In 'while loop' solution will be converged with the error less
                % than 1e-4
                while(error > tol)
                        % For spatial discretization, another 'for loop' is used.
                        % In this loop, the numerical scheme derived as shown in
                        % the introduction was used. 'i-loop' is for x-direction,
                        % 'j-loop' is for y-direction
                        for i = 2:nx-1
                            for j = 2:ny-1
                                H = (T(i-1,j) + Told(i+1,j));
                                V = (T(i,j-1) + Told(i,j+1));
                                T(i,j) = (term1* Tprev(i,j) + term2*H + term3*V);
                            end
                        end
                    error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
                    Told = T; % This logic helps to write the values of T from the previous step
                end
            Tprev = T; % This logic helps to write the values of T from the previous step
            GS_iter = GS_iter + 1; % This is the counter for number of iterations
        end
   time = toc;  % This command is to begin the stop watch and is saved in 'time'     
   figure(2)
   contourf(x,y,T,'showText','on')
   grid on
   colormap(jet)    
   title_text = sprintf('2D Heat Conduction Using Gauss-Seidel Method \n Implicit Scheme \n Iterations = %d \n Time = %f sec', GS_iter,time);
   title(title_text)
   set(gca,'YDIR','reverse')
   pause(0.03)
end
 
% SOR method
if iterative_solver == 3
    tic; % This command is to begin the stop watch
    SOR_iter =1;
        % For Time-marching, 'for loop(k)' is used
        for k = 1:nt
            error = 9e9;
                % In 'while loop' solution will be converged with the error less
                % than 1e-4
                while(error > tol)
                        % For spatial discretization, another 'for loop' is used.
                        % In this loop, the numerical scheme derived as shown in
                        % the introduction was used. 'i-loop' is for x-direction,
                        % 'j-loop' is for y-direction
                        for i = 2:nx-1
                            for j = 2:ny-1
                                H = (T(i-1,j) + Told(i+1,j));
                                V = (T(i,j-1) + Told(i,j+1));
                                T(i,j) = ((1-w)*Told(i,j) + w*(term1* Tprev(i,j) + term2*H + term3*V));
                            end
                        end
                    error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
                    Told = T; % This logic helps to write the values of T from the previous step
                end
            Tprev = T; % This logic helps to write the values of T from the previous step
            SOR_iter = SOR_iter + 1; % This is the counter for number of iterations
        end
    time = toc;  % This command is to begin the stop watch and is saved in 'time'  
    figure(3)
    contourf(x,y,T,'showText','on')
    grid on
    colormap(jet)    
    title_text = sprintf('2D Heat Conduction Using SOR Method \n Implicit Scheme \n Iterations = %d \n Time = %f sec', SOR_iter,time);
    title(title_text)
    set(gca,'YDIR','reverse')
    pause(0.03)
end




% Solving two-dimensional unsteady conduction equation using iterative
% solvers for explicit and implicit time scheme to check the stability % of the solution
% Given equation is : dT/dt - alpha(d^2T/dx^2 + d^2T/dy^2) = 0
% Left Boundary = 400 K; Right Boundary = 800 K; Top Boundary = 600 K; Bottom Boundary = 900 K;
% nx = ny
clear all; close all; clc;
 
% Given inputs
L = 1; % Length of the domain in m
B = 1; % Width of the domain in m
dt = 0.001; % given timestep in sec
t = 0.25; % given time in sec
alpha = 1.4; % Diffusivity
 
% Number of time steps can be calculated by nt = given time/timestep size =
% t/dt
nt = 350;
 
% Nodes 
nx  = 10; % Nodes in x-direction
ny  = 10; % Nodes in y-direction
 
% Defining space variables
x = linspace(0,L,nx);
y = linspace(0,B,ny);
dx  = x(2) - x(1); %  defining grid spacing in x-direction
dy = y(2) - y(1); %  defining grid spacing in y-direction
 
%Initialization
T = ones(nx,ny);
 
% Given BCs
T(:,1) = 400;
T(:,end) = 800;
T(1,:) = 600;
T(end,:) = 900;
 
% Defining corner nodes using BCs
T(1,1) = 500; % (600+400)/2
T(end,1) = 650; % (400+900)/2
T(1,end) = 700; % (600+800)/2
T(end,end) = 850; % (900+800)/2 
 
% Creating a copy of T
Told = T; % This logic helps to write the values of T from the previous iteration
Tprev = T; % This logic helps to write the values of T from the previous timestep
 
% CFL number: For stable solutions, it is advisable to keep CFL <= 0.5 for explicit method.Here
% CFL = CFL_x + CFL_y
k1 = alpha*dt/dx^2; % let CFL_x = k1 = alpha*dt/dx^2
k2 = alpha*dt/dy^2; % let CFL_x = k1 = alpha*dt/dx^2
 
 
iterative_solver = input('Solver number: Jacobi = 1 Gauss-siedel = 2 SOR = 3 Explicit = 4 \n Enter solver number:');
 
% This information is given in the introduction
term1 = 1/(1 + 2*k1 + 2*k2);
term2 = k1*term1;
term3 = k2*term1;
w = 1.5; % Over-relaxation factor for SOR method
 
%Initialization of error and tolerance values 
error = 9e9;
tol = 1e-4;
 
 % Jacobi method
if iterative_solver == 1
    tic; % This command is to begin the stop watch
    jacobi_iter = 1;
        % For Time-marching, 'for loop(k)' is used 
        for k = 1:nt 
            error = 9e9;
                % In 'while loop' solution will be converged with the error less
                % than 1e-4
                while(error > tol) 
                        % For spatial discretization, another 'for loop' is used.
                        % In this loop, the numerical scheme derived as shown in
                        % the introduction was used. 'i-loop' is for x-direction,
                        % 'j-loop' is for y-direction
                        for i = 2:nx-1
                            for j = 2:ny-1
                                H = (Told(i-1,j) + Told(i+1,j));
                                V = (Told(i,j-1) + Told(i,j+1));
                                T(i,j) = (term1* Tprev(i,j) + term2*H + term3*V);
                            end
                        end
                    error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
                    Told = T; % This logic helps to write the values of T from the previous iteration
                end
            Tprev = T; % This logic helps to write the values of T from the previous timestep
            jacobi_iter = jacobi_iter + 1; % This is the counter for number of iterations
         end
    time = toc;  % This command is to begin the stop watch and is saved in 'time'   
    figure(1)
    contourf(x,y,T,'showText','on');
    grid on
    colormap(jet)
    title_text = sprintf('2D Unsteady Heat Conduction Using Jacobi Method \n Implicit Scheme \n Iterations = %d \n Time = %f sec', jacobi_iter,time);
    title(title_text)
    set(gca,'YDIR','reverse')
    pause(0.03)
end
 
% Gauss-siedal method
if iterative_solver == 2
    tic; % This command is to begin the stop watch
    GS_iter =1;
        % For Time-marching, 'for loop(k)' is used
        for k = 1:nt
            error = 9e9;
                % In 'while loop' solution will be converged with the error less
                % than 1e-4
                while(error > tol)
                        % For spatial discretization, another 'for loop' is used.
                        % In this loop, the numerical scheme derived as shown in
                        % the introduction was used. 'i-loop' is for x-direction,
                        % 'j-loop' is for y-direction
                        for i = 2:nx-1
                            for j = 2:ny-1
                                H = (T(i-1,j) + Told(i+1,j));
                                V = (T(i,j-1) + Told(i,j+1));
                                T(i,j) = (term1* Tprev(i,j) + term2*H + term3*V);
                            end
                        end
                    error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
                    Told = T; % This logic helps to write the values of T from the previous step
                end
            Tprev = T; % This logic helps to write the values of T from the previous step
            GS_iter = GS_iter + 1; % This is the counter for number of iterations
        end
   time = toc;  % This command is to begin the stop watch and is saved in 'time'     
   figure(2)
   contourf(x,y,T,'showText','on')
   grid on
   colormap(jet)    
   title_text = sprintf('2D Heat Conduction Using Gauss-Seidel Method \n Implicit Scheme \n Iterations = %d \n Time = %f sec', GS_iter,time);
   title(title_text)
   set(gca,'YDIR','reverse')
   pause(0.03)
end
 
% SOR method
if iterative_solver == 3
    tic; % This command is to begin the stop watch
    SOR_iter =1;
        % For Time-marching, 'for loop(k)' is used
        for k = 1:nt
            error = 9e9;
                % In 'while loop' solution will be converged with the error less
                % than 1e-4
                while(error > tol)
                        % For spatial discretization, another 'for loop' is used.
                        % In this loop, the numerical scheme derived as shown in
                        % the introduction was used. 'i-loop' is for x-direction,
                        % 'j-loop' is for y-direction
                        for i = 2:nx-1
                            for j = 2:ny-1
                                H = (T(i-1,j) + Told(i+1,j));
                                V = (T(i,j-1) + Told(i,j+1));
                                T(i,j) = ((1-w)*Told(i,j) + w*(term1* Tprev(i,j) + term2*H + term3*V));
                            end
                        end
                    error = max(max(abs(Told-T))); % This step continues until error is less than the given tolerance value in-order-to the solution to be converged
                    Told = T; % This logic helps to write the values of T from the previous step
                end
            Tprev = T; % This logic helps to write the values of T from the previous step
            SOR_iter = SOR_iter + 1; % This is the counter for number of iterations
        end
    time = toc;  % This command is to begin the stop watch and is saved in 'time'  
    figure(3)
    contourf(x,y,T,'showText','on')
    grid on
    colormap(jet)    
    title_text = sprintf('2D Heat Conduction Using SOR Method \n Implicit Scheme \n Iterations = %d \n Time = %f sec', SOR_iter,time);
    title(title_text)
    set(gca,'YDIR','reverse')
    pause(0.03)
end
 
 
% Explicit scheme
if iterative_solver == 4
    tic; % This command is to begin the stop watch
    Exp_Iteration =1;
        % For Time-marching, 'for loop(k)' is used
        for k = 1:nt
                % For spatial discretization, another 'for loop' is used.
                % In this loop, the numerical scheme derived as shown in
                % the introduction was used. 'i-loop' is for x-direction,
                % 'j-loop' is for y-direction
                for i = 2:nx-1
                    for j = 2:ny-1
                        term11 = k1* (Told(i-1,j) - 2*Told(i,j) + Told(i+1,j));
                        term22 = k2* (Told(i,j-1) - 2*Told(i,j) + Told(i,j+1));
                        T(i,j) = Told(i,j) + term11 + term22 ;
                    end
                end
            Told = T; % This logic helps to write the values of T from the previous step
            Exp_Iteration = Exp_Iteration + 1; % This is the counter for number of iterations    
        end
    time = toc;  % This command is to begin the stop watch and is saved in 'time'  
    figure(4)
    contourf(x,y,T,'showText','on');
    grid on
    colormap(jet)
    title_text = sprintf('2D Unsteady Heat Conduction \n Explicit Scheme \n Iterations = %d \n Time = %f sec', Exp_Iteration,time);
    title(title_text)
    set(gca,'YDIR','reverse')
    pause(0.03)
end
