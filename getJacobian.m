function J = getJacobian(m, v)
%GETJACOBIAN linearizes relation of joint space input and cartesian output

    % Author: Till Blaha; Date: 2019-05-15
    
    % compute jacobian around the current member and upright positions;
    % this matrix indicates how a differential joint input (rotation of the
    % links or translation and rotation of the upright) would affect the
    % differential change in the cartesian positions of the pickup points.
    % In other terms, this is the linearization of the cartesian response
    % of the pickup points to rotation of the links and translation or
    % rotation of the upright.
    
    % jacobian will be 18x18 matrix:
    % 3 translation of rigid body upright ...
    %   + 3 rotations of rigid body  upright ...
    %   + 6 links * 2 rotations around inner pup per link = 18 inputs
    % 6 outer pups * 3 cartesian coordinates per pup = 18 outputs
    
    %%% format %%%
    % column meanings
    %   [x y z]_upright | [theta delta gamma]_upright | [gamma theta]_LR | [gamma theta]_LF | etc...
    % row meanings
    %   [x; y; z]_LRO | [x; y; z]_LFO | etc...
    J = zeros(18);
    
    %% 1. translation of upright
    
    % simple: derivative of linear motion in the global x-y-z frame is the
    % exact same for each pickup point on the upright. So these first 3
    % columns are the identity matrix repeated for as many rows as
    % necessary
    J(:,1:3) = repmat(eye(3), 6, 1);
    
    
    %% 2. rotation of the upright
    % ...around the LRO by choice of the vectors v
    
    % the linearization of rotational motion around a fixed axis is simply
    % cross(2*pi*axis_unit_vector, vector_of_remote_point_to_rotation_centre)
    
    % the cross function in matlab only works with matrices if we have the
    % exact same number of vectors in either matrix. The Kronecker product
    % is used to triple the occurence of the vectors v since we need to do
    % 3 cross products with different axes for each point on the upright.
    % example kron([1 2; 3 4], [1 1]) == [1 1 2 2; 3 3 4 4]
    up_rots = cross(2*pi*repmat(eye(3), 1, 6), kron(v, [1 1 1]));
    
    % So far, I couldn't figure out a nice way to reshape these cross
    % products in batches of 3... so this feels a bit beun but works and
    % should be fast
    J(:,4:6) = [up_rots(:,1:3); up_rots(:,4:6); up_rots(:,7:9); up_rots(:,10:12); up_rots(:,13:15); up_rots(:,16:18)];
    
    
    %% 3. rotation of all the links

    % again the linearization of rotational motion around a fixed axis is
    % simply cross(2*pi*axis_unit_vector, vector_of_remote_point_to_rotation_centre)
    % also, for the suspension links we only need 2 rotations as indicated
    % above: the 3rd one is rotation around the links own axis, which is
    % meaningless in this simulation. Also, the order is a bit arbitrary,
    % but well...
    diagonal = cross(2*pi*repmat([[0 0 1]' [1 0 0]'], 1, 6), kron(m, [1 1]));
    % diagonal is a 12x3 matrix (2 3-component vectors per pickup point)
    
    % this time, each pair of cross products above "reacts" to a different
    % angle, so I used the block diagonal function to arrange it
    % accordingly in the matrix
    J(:,7:18) = blkdiag(diagonal(:,1:2), diagonal(:,3:4), diagonal(:,5:6), diagonal(:,7:8), diagonal(:,9:10), diagonal(:,11:12));
    
    

end