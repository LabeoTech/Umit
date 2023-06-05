 function MI =  mutual_information(counts)
    % This function calculates the Mutual Information from the joint
    % histogram of two images.
    
    pxy = counts / sum(counts(:));
    px = sum(pxy, 2); % marginal for x over y
    py = sum(pxy, 1); % marginal for y over x
    px_py = px * py; % Multiply marginals
    
    % Only non-zero pxy values contribute to the sum
    nzs = pxy > 0;
    
    % Now we can do the calculation using the pxy, px_py 2D arrays
    MI = sum(pxy(nzs) .* log(pxy(nzs) ./ px_py(nzs)));   
end