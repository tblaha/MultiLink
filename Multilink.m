%% Multilink.m

% Author: Till Blaha; Date: 2019-05-15
% This is a proof of concept for a fast and robust solver of Multi Link
% kinematics. It is not inverse kinematics, so you cannot specify target 
% toe or travel values or the like, but you can sweep steering or push/pull
% rod inputs

%% Changelog

% 2019-05-15:   crated, very proof of concept, no meaningful suspension
%               inputs yet. Not meaningful outputs yet, but kinematics
%               solver seems to give accurate results




%%

clear
close all

addpath(genpath('./'))

% do you want to animate the sweep? it is slooww
animate = true;


%% arbitrary pickup points

pup.LRI = [100  300 110]';
pup.LFI = [-100 300 110]';
pup.URI = [100  300 250]';
pup.UFI = [-100 300 250]';
pup.SI  = [0    300 500]';
pup.TI  = [70   300 250]';

pup.SO  = [0   600 300]';
pup.TO  = [30  600 250]';
pup.LRO = [15  600 130]';
pup.LFO = [-15 600 130]';
pup.URO = [10  600 270]';
pup.UFO = [-10 600 270]';

% get suspension link vectors and their lengths
[minit, l] = getMembers(pup);

% cartesian coordinate vectors used for the iterations later (both 18x1)
% since 6 points with 3 coordinates each
Ipups0 = [pup.LRI; pup.LFI; pup.URI; pup.UFI; pup.TI; pup.SI]; % inner
Opups0 = [pup.LRO; pup.LFO; pup.URO; pup.UFO; pup.TO; pup.SO]; % outer


%% define upright vectors which contain the rigid-body geometry of the upright

% use pup.LRO as origin. This will be explained in the next chunk of code
v = zeros(3,6);
v(:,1) = zeros(3,1);
v(:,2) = pup.LFO - pup.LRO;
v(:,3) = pup.URO - pup.LRO;
v(:,4) = pup.UFO - pup.LRO;
v(:,5) = pup.TO  - pup.LRO;
v(:,6) = pup.SO  - pup.LRO;
vinit = v; % v will be changed later in iteration and we want vinit to remain unchanged

% a bit abstract: define a vector basis that is unique to the upright and
% will be used to "encode" the global wheel plane normal ([0 1 0]' when not
% steering/cambering) in a local upright coordinate system and later
% transform it back to evaluate camber/toe/etc and help with plotting.
% origin of this coordinate frame will be pup.LRO
upB = [v(:,3)/norm(v(:,3)) v(:,5)/norm(v(:,5)) cross(v(:,3), v(:,5)) / norm(cross(v(:,3), v(:,5)))];

% define the wheel normal and wheel centre and geometric contact patch in
% unperturbed state
wh_n = [0 1 0]';
wh_0 = [0 650 200]';
wh_0_LRO = wh_0 - pup.LRO;
cp_geo = [0 650 0]';
cp_geo_LRO = cp_geo - pup.LRO;

% transform to said upright coordinate system
wh_0_LRO_upB   = upB\wh_0_LRO;
wh_n_upB       = upB\wh_n;
cp_geo_LRO_upB = upB\cp_geo_LRO;



%% set up sim

% steering rack travel input
steering_range = 12:-1:-12;

% set numerical damping factor of the iterations
lambda = 1.72; % when recomputing the jacobian; this is optimal for some reason based on trial and error (1.72 for fwd diff, 3.44 for central diff)
%lambda = 1; % when reusing the jacobian (sounds faster, but is slower)

% pre-allocate outputs
wh_0_dyn   = zeros(3,  length(steering_range));
wh_n_dyn   = zeros(3,  length(steering_range));
cp_geo_dyn = zeros(3,  length(steering_range));
Opups_dyn  = zeros(18, length(steering_range));
Ipups_dyn  = zeros(18, length(steering_range));


%% simulate
tic
for j = 1:length(steering_range)
    %% inputs
    
    v=vinit; % v will be changed later, but vinit should not be modified
    m=minit; % m will be changed later, but minit should not be modified
    
    % tie rod motion
    Ipups_dyn(:,j)  = Ipups0;
    Ipups_dyn(14,j) = Ipups_dyn(14,j) + steering_range(j); 
    
%     % direct push-link extension:
%     extension = 5; %mm, for instance
%     m(:,6) = m(:,6) * (norm(m(:,6)) + extension) / norm(m(:,6));

    % rocker rotation
    % ToDo

    
    %% simulation initialization
    
    % remember all vectors for sanity checking later
    m0=m;
    v0=v;

    % initialize the iteration variables
    Opups_u = Opups0;
    Opups_l = Ipups_dyn(:,j) + m(:);
    
    % pre allocate error and command variable
    E = ones(18,1);
    jsc_prev = zeros(18,1); % jsc stands for joint space command

    %% iterate
    i = 1;
    while max(abs(E)) > 1e-2 && i < 70

        % evaluate error between current position of the outer pickup
        % points on the upright and current position of the outer pickup
        % points on the ends of the links
        % this should be zero and the end of the iterations!
        E = (Opups_u - Opups_l); 
        
        % recompute jacobian around the current member and upright
        % positions; this matrix indicates how a differential joint input
        % (rotation of the links or translation and rotation of the
        % upright) would affect the differential change in the cartesian
        % positions of the pickup points. In other terms, this is the
        % linearization of the cartesian response of the pickup points to
        % rotation of the links and translation or rotation of the upright
        J = getJacobian(m, v);

        % since d_cartesian = J*d_joint_space as described before, we can
        % issue the next commands to decrease the error as
        % joint_space_command = J^(-1)*CartesianError the negative for the
        % first 6 joints arises from the fact that those 6 are the DoFs of
        % the upright and we need to "join" the upright with the links so
        % if the error is positive for some y coordinate of some point,
        % then the upright must move left (negative y) and the links must
        % be rotated such that the y value for that coordinate is more
        % positive. This is where the method feels a bit beun and I don't
        % know why I can crank the "damping" factor that I expected to be
        % smaller than 0.5 all the way up to 1.72 and still get
        % convergence... but hey, it works
        jsc = J\E .* [-ones(6,1); ones(12,1)] * lambda;

        % apply this estimated jsc to the system using the full non-linear
        % set of equations of motion. Sounds fancier than it is, this is
        % basically just rotation matrices and translations (adding and 
        % subtracting)
        % returns new link and upright positions. Also returns new member
        % vectors m and global upright geometry vectors v.
        [Opups_l, Opups_u, m, v] = execute_jsc(jsc, Ipups_dyn(:,j), Opups_u, m, v); % recompute jacobian

        % it seemed to me that computing the jacobian everytime would be
        % expensive, and I could get away with not doing it and keeping the
        % jacobian computed around the initial conditions but it was slower
        % since it needed a lot more iterations and the execute_jsc
        % function is also quite expensive. so leave the commented
        %[Opups_r, Opups_u, m, v] = execute_jsc(jsc+jsc_prev, Ipups_dyn(:,j), Opups_u, m0, v0); % reuse jacobian
        %jsc_prev = jsc_prev+jsc;

        % sanity check that the members didn't change lengths all of a
        % sudded or the upright geometry is distorted
        assert(norm(vecnorm(m) - vecnorm(m0)) < 1e-10)
        assert(norm(vecnorm(v) - vecnorm(v0)) < 1e-10)

        i=i+1;
    end
    % save the converged pickup points. The choice of _l or _u is arbitrary
    % and actually the mean between the two could make this entire thing
    % more precise? Or does that actually violate the asserts?
    if max(abs(E)) < 1e-2
        Opups_dyn(:,j) = (Opups_l);
    else
        warning(strcat("Probable lock-up at ", num2str(steering_range(j))) )
        Opups_dyn(:,j) = zeros(size(Opups_l));
    end
    
    
    
    
    %% post-processing -- compute wheelplane

    % recompute the upright basis transformation matrix with the final upright
    % position vectors
    upB_dyn = [v(:,3)/norm(v(:,3)) v(:,5)/norm(v(:,5)) cross(v(:,3), v(:,5)) / norm(cross(v(:,3), v(:,5)))];
    
    % apply the wheel centre position and wheel plane normal given in the
    % upright basis to that matrix to get it in global coordinates
    % wheel plane geometric centre
    wh_0_LRO_dyn    = upB_dyn * wh_0_LRO_upB;
    wh_0_dyn(:,j)   = Opups_dyn(1:3,j) + wh_0_LRO_dyn;
    
    % wheel plane normal
    wh_n_dyn(:,j)   = upB_dyn * wh_n_upB;
    
    % tire geometric contact patch
    cp_geo_LRO_dyn  = upB_dyn * cp_geo_LRO_upB;
    cp_geo_dyn(:,j) = Opups_dyn(1:3,j) + cp_geo_LRO_dyn;
    
    
    %% post-processing -- "measure" toe, camber, etc
    
    % todo
    
end
toc







%% animate

if animate
    filename = 'Plots/MultiLinkSim.gif';
    mkdir("Plots")

    % indices for forwards-backwards animation
    js = [1:length(steering_range) length(steering_range)-1:-1:2];
    
    % time per frame, a bit arbitrary...
    dt = abs(mean(diff(steering_range)))/10;

    % go
    tic
    h = figure('Position', [100 100 800 800], 'Visible', 'on');
    ax = axes;
    for j = js
        
        % get the centre of the wheel at the inner side and the outer side
        % which will be used to plot a cylinder (this is a 3x2 matrix)
        width_tire = 160; % arbitrary
        wheelCs = wh_0_dyn(:,j) + wh_n_dyn(:,j) .* width_tire .* [0.5 -0.5];

        % use this weird geom3d toolbox to create what they call "edges"
        E = createEdge3d(reshape(Ipups_dyn(:,j), 3, 6)', reshape(Opups_dyn(:,j), 3, 6)');

        % clear figure (there must be a better way of doing this... 
        clf
        hold on
            % draw links
            drawEdge3d(E, 'color', 'k', 'LineWidth', 2)
            % draw "tire"
            drawCylinder([wheelCs(:)' 150],  'open', 'EdgeColor', 'k', 'FaceAlpha', 0.5)
        hold off
        axis equal
        view([43 38])
        %view([0 90])
        xlim([-200 200])
        ylim([270 800])
        zlim([0 400])

        % needed for animation
        drawnow

        % save to gif
        frame = getframe(h); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if j == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,  'DelayTime',dt); 
        else 
            imwrite(imind,cm,filename,'gif','WriteMode','append',  'DelayTime',dt); 
        end 
    end
    toc
end

