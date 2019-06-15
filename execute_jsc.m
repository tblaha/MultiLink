function [Opups_r, Opups_u, m, v] = execute_jsc(jsc, Ipups, Opups_u, m, v)
%EXECUTE_JSC applies a joint space command to a multilink suspension

    % Author: Till Blaha; Date: 2019-05-15
    
    % This function is split in two parts: 
    % 1. apply the rotations and translations to the upright computing the
    % next positions of the outer pickup locations on the upright.
    % 2. apply the rotations of the suspension links to get the next
    % positions of the outer ends (or outer pickup points) of the links 
    
    %% 1. upright
    
    % translate the 6 outer pups using the first three (already cartesian)
    % joint space variables
    Opups_u = Opups_u + repmat(jsc(1:3), 6, 1);

    % rotate 6 outer pups around LRO. Rotation sequency is weird in this
    % code and doesn't really matter since we will evaluate camber and toe
    % differently later anyway, but this is the way I have it written down
    rmat = roty(jsc(5)*180/pi)*rotx(jsc(4)*180/pi)*rotz(jsc(6)*180/pi);
    
    % execute the rotation by having the 3x3 rotation matrices on the
    % diagonal of a big 12x12 rotation matrix. See 2. for a better
    % explanation of this since it is basically the same
    vvec_new = blkdiag(rmat, rmat, rmat, rmat, rmat, rmat) * v(:);
    
    % add the new vectors to the origin (RLO which is 1:3)
    Opups_u = repmat(Opups_u(1:3), 6, 1) + vvec_new;

    % assembly new upright geometry vectors (this will be sanity checked
    % later such that it can't violate the rigid upright geometry by a
    % mistake in the above)
    v = [vvec_new(1:3) vvec_new(4:6) vvec_new(7:9) vvec_new(10:12) vvec_new(13:15) vvec_new(16:18)];
    
    
    %% 2. suspension links
    
    % rotate all the links. This is done with a big rotation matrix:
    % m_new = big_rot_mat * m(:) instead of 6 separate 3x3 matrices for
    % slight speed up
    
    % preallocate the 12x12 matrix (6 links and 2 rotations about the inner
    % pickup points each)
    links_rmats = zeros(12);
    for i = 1:6
        % rotation sequency for each link is first x, then z, but it
        % doesn't really matter, since the Jacobian linearization doesn't
        % take this into account anyway and the resulting joint commands
        % are not used in the code later.
        % sorry for the indexing mess... this basically assigned 3x3 blocks
        % along the diagonal of the 12x12 matrix using the odd joint
        % commands 7,9,11,13,15,17 for z rotation and the even commands
        % 8,10,12,14,16,18 for the x rotation of the respective link
        links_rmats(i*3-2:i*3, i*3-2:i*3) = rotz(jsc(i*2-1+6)*180/pi)*rotx(jsc(i*2+6)*180/pi);
    end
    
    % execute the rotation to the suspension members 
    mvec_new = links_rmats * m(:);
    
    % add the new suspension member vectors to their respective inner
    % pickup points
    Opups_r = Ipups + mvec_new;
    
    % assembly new suspension link vectors (this will be sanity checked
    % later such that it can't violate the link lengths by a mistake in the
    % above)
    m = [mvec_new(1:3) mvec_new(4:6) mvec_new(7:9) mvec_new(10:12) mvec_new(13:15) mvec_new(16:18)];
    

end