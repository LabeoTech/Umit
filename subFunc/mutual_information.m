 function MI =  mutual_information(counts)
    % This function calculates the Mutual Information from the joint
    % histogram of two images.
    
    pXY = counts./sum(counts(:));
    pX = sum(pXY,2, 'omitnan');
    pY = sum(pXY,1, 'omitnan');
    
    MI = sum(pXY.*log(pXY./(pX*pY)),'all','omitnan');    
end