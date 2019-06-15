function [m, l] = getMembers(pup)

    m = zeros(3, 6);

    m(:,1) = pup.LRO - pup.LRI; 
    m(:,2) = pup.LFO - pup.LFI; 
    m(:,3) = pup.URO - pup.URI; 
    m(:,4) = pup.UFO - pup.UFI; 
    m(:,5) = pup.TO  - pup.TI;  
    m(:,6) = pup.SO  - pup.SI;  
    
    l = vecnorm(m);


end