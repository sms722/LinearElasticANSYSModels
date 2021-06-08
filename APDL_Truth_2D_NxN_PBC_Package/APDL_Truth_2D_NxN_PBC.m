% -----------------------------------------------------------------
% ANSYS APDL 2D NxN TRUTH CODE, with Periodic Boundary Conditions
% This code models members as cylinders with varying radii, and 3D 
% quadratic finite elements are used
% This code considers only NxN 2D trusses
% Function for calculating metamaterial properties for a 2D NxN truss
% All lengths are in [m], all stresses and moduli are in [Pa] 
% Each node has two degrees of freedom (DOF), x and y.  The (node number*2) 
%    gives the number for that node's yDOF, and (yDOF - 1) is the xDOF
% -----------------------------------------------------------------
% This code uses MATLAB to call ANSYS APDL for FEM analysis.
% First, MATLAB writes the design parameters into input file 
% (e.g. para_in.txt)
% Then, MATLAB calls ANSYS APDL to excute the APDL file 
% The APDL file reads the parameters from the input file and writes the 
% analysis result to the output file (e.g. para_out.txt)
% Finally, MATLAB reads the results from the output file
% -----------------------------------------------------------------
% CREDITS FOR ORIGINAL CODE PACKAGE:
% https://github.com/zhandawei/Call_ANSYS_in_MATLAB
% Dawei Zhan
% zhandawei{at}swjtu.edu.cn
% -----------------------------------------------------------------
% Modified/adapted for KD3M2 project by Srikar Srivatsa
% -----------------------------------------------------------------
% Sample values for testing script:
%{
clc;    
close all; 
clear;
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;4,5;5,6;5,7;5,8;5,9];
CA = [1,2;2,3;3,6;6,9;8,9;7,8;4,7;1,4;1,5;2,5;3,5;5,6;5,7;5,8];
CA = [1,2;3,6;6,9;7,8;4,7;1,4;1,5;2,5;2,9;4,5;4,8;5,8];
CA = [1,2;2,3;8,9;7,8;1,5;2,5;3,5;4,5;4,8;5,6];
sel = 0.01; E = 1816200; sidenum = 3;
rvar = (250*(10^-6)).*ones(size(CA,1),1); 
%}

function C = APDL_Truth_2D_NxN_PBC(sel,sidenum,rvar,E,CA)
    % Initialize stiffness tensor
    C = [0,0,0;0,0,0;0,0,0];
    
    % Generate vector with nodal coordinates
    NC = generateNC(sel,sidenum);
    
    % Calculate predetermined inputs
    [sortCA,sortrvar,Lvec,WPangles,edgetype] = ... 
            predetInputs(CA,NC,rvar,sidenum);
    
    % Iterate through once for each strain component
    for y = 1:1:3
        % Define strain vector: [e11, e22, e12]'
        strainvec = [0;0;0];

        % Set that component equal to a dummy value (0.01 strain), set all 
        % other values to zero
        strainvec(y) = 0.001; 
        strainvec(3) = strainvec(3)*2;
        
        % Populate Input Vector
        qsel = sel./(10^-3); qrvar = sortrvar./(10^-6); 
        qLvec = Lvec./(10^-3);
        input_vec=[qsel;E;size(CA,1);y;sidenum;qrvar;qLvec;WPangles;...
                   edgetype;sortCA(:,1);sortCA(:,2)];

        % Define paths and input/output files 
        ANSYS_path=...
 'C:\Program Files\ANSYS Inc\v211\ansys\bin\winx64\ANSYS211';
        APDL_name='Truth2D_NxN.txt';
        input_file_name='para_in.txt';
        output_file_name='para_out.txt';
        ANSYS_path=strcat('"',ANSYS_path,'"');
        disp('Current stage: ');disp(y);

        % Write input vector to para_in.txt
        %writematrix(input_vec,input_file_name);
        dlmwrite(input_file_name,input_vec,'delimiter',' ',...
           'precision','%14.6f','newline','pc');
        % Call ANSYS APDL, perform FEA
        status = system(sprintf('%s -b -p aa_r_me -i %s -o out.txt',...
            ANSYS_path,APDL_name));
        disp('Simulation status: ');disp(status);
        % Read the results from para_out.txt
        output_vec=load(output_file_name)';
        F_x = output_vec(1);
        F_y = output_vec(2);
        F_xy = output_vec(3);
        FBank(y,:) = [F_x,F_y,F_xy];

        % Calculate stress vector
        stressvec = (1/(sel*2*r)).*[F_x;F_y;F_xy];

        % Use strain and stress vectors to solve for the corresponding row 
        % of the C matrix
        Cdummy = stressvec/strainvec;
        C(:,y) = Cdummy(:,y);
    end
end

%----------%
% FUNCTION TO GENERATE NODAL COORDINATES BASED ON GRID SIZE
function NC = generateNC(sel,sidenum)
    notchvec = linspace(0,1,sidenum);
    NC = [];
    for i = 1:1:sidenum
        for j = 1:1:sidenum
            NC = [NC;notchvec(i),notchvec(j)];
        end
    end
    NC = sel.*NC;
end

% FUNCTION TO CALCULATE PREDETERMINED INPUTS
function [sortCA,sortrvar,Lvec,WPangles,edgetype] = ...
                predetInputs(CA,NC,rvar,sidenum)
    % Sort CA, rvar
    [sortCA,index] = sortrows((sort(CA'))');
    sortrvar = rvar(index);
    
    % Find member lengths, angles
    x1 = NC(sortCA(:,1),1); x2 = NC(sortCA(:,2),1);
    y1 = NC(sortCA(:,1),2); y2 = NC(sortCA(:,2),2);
    Lvec = sqrt(((x2-x1).^2)+((y2-y1).^2));
    WPangles = atan((y2-y1)./(x2-x1)); %WPangles = WPangles + 90;
    
    % Identify four sets of edge nodes
    emt1nodes = 1:1:sidenum;
    emt2nodes = ((sidenum^2)-sidenum+1):1:(sidenum^2);
    emt3nodes = 1:sidenum:((sidenum^2)-(sidenum)+1);
    emt4nodes = sidenum:sidenum:(sidenum^2);
             
    % Identify members connecting solely to edge nodes
    emt1connA = ismember(sortCA(:,1),emt1nodes);
    emt1connB = ismember(sortCA(:,2),emt1nodes);
    emt1log = (emt1connA & emt1connB);
    emt2connA = ismember(sortCA(:,1),emt2nodes);
    emt2connB = ismember(sortCA(:,2),emt2nodes);
    emt2log = 2.*(emt2connA & emt2connB);
    emt3connA = ismember(sortCA(:,1),emt3nodes);
    emt3connB = ismember(sortCA(:,2),emt3nodes);
    emt3log = 3.*(emt3connA & emt3connB);
    emt4connA = ismember(sortCA(:,1),emt4nodes);
    emt4connB = ismember(sortCA(:,2),emt4nodes);
    emt4log = 4.*(emt4connA & emt4connB);
    edgetype = emt1log + emt2log + emt3log + emt4log;
    
    % Modify work-plane angles of edge members
    WPangles(edgetype==1) = deg2rad(270);
    WPangles(edgetype==2) = deg2rad(90);
    WPangles(edgetype==3) = deg2rad(0);
    WPangles(edgetype==4) = deg2rad(180);
end




