function [results] = t2FC(xyData1, xyData2, varargin)
    % Syntax: [results] = tSquaredFourierCoefs(xyData,testMu,alphaVal)
    %
    % Returns the results of running Hotelling's t-squared test that the mean
    % of the 2D data in xyData is the same as the mean specified in testMu at
    % the specified alphaVal (0-1).
    %
    % Based on Anderson (1984) An Introduction to Multivariate Statistical 
    % Analysis, Wiley
    %
    % In:
    %   xyData - n x 2 x q matrix containing n data samples
    %            if q = 1, data will be tested against zero    
    %            if q = 2, a paired t-test will be performed,
    %            testing xyData(:,:,1) against xyData(:,:,2). 
    %            2 is the maximum length of the third dimension. 
    %            Function assumes that the 2D data in xyData(:,:,1)
    %            and the optional xyData(:,:,2) are Fourier 
    %            coefficients organized with the real coefficients 
    %            in the 1st column and the imaginary
    %            coefficients in the 2nd column. Rows = samples.
    % 
    % <options>:
    %   testMu - 2-element vector specifying the mean to test against ([0,0]) 
    % 
    %   alphaVal - scalar indicating the alpha value for testing ([0.05])
    %
    %   pdStyle - do PowerDiva style testing (true/[false])
    %
    % Out:
    %
    %   results - struct containing the following fields:
    %             alpha, tSqrdCritical, tSqrd, pVal, H (0 if can't reject null; 1 if
    %             rejected the null hypothesis)

    opt	= ParseArgs(varargin,...
        'testMu'		, [0,0], ...
        'alphaVal'		, 0.05	, ...
        'pdStyle',        false ...
        );

    dims = size(xyData1);
    if dims(1) < 2
        error('input data must contain at least 2 samples');
    end
    if dims(2) ~= 2
        error('input data must be a matrix of 2D row samples');
    end
    if length(dims) < 3 % if no third dimension
        xyData1(:,:,2) = zeros(size(xyData1));
        xyData2(:,:,2) = zeros(size(xyData2));
    elseif dims(3) > 2
         error('length of third dimension of input data may not exceed two')
    else
    end

    if length(opt.testMu) ~= 2
        error('testMu should be a 2D vector');
    end

    
    % remove NaNs
    d1 = xyData1(:,:,1) - xyData1(:,:,2);  
    d2 = xyData2(:,:,1) - xyData2(:,:,2);
    
    dData1 = d1(~any(isnan(d1),2),:);
    dData2 = d2(~any(isnan(d2),2),:);
    
    n1 = size(dData1, 1);
    n2 = size(dData2, 1);
    
    try
        [sMu1, sCov1] = eigFourierCoefs(dData1);
        [sMu2 ,sCov2] = eigFourierCoefs(dData2);
    catch
        fprintf('The covariance matrix of xyData could not be calculated, probably your data do not contain >1 sample.');
    end

    results.alpha = opt.alphaVal;
    p = 2;
    df1 = p;
    df2 = n1 + n2 - p - 1;    
    
    t0Sqrd = (((n1 + n2 - 2)*p)/df2)*finv( 1 - opt.alphaVal, df1, df2);
  
    results.tSqrdCritical = t0Sqrd;
    
    %% if covariance matrices are the same:
    if (sCov1 == sCov2)
        S_pool = ( (n1 - 1)*sCov1 + (n2 - 1)*sCov2)/(n1 + n2 - 2);
        diff_S_pool = S_pool*(1/n1 + 1/n2);      
        diff_mu = sMu1 - sMu2;
        
        invCovMat = inv(diff_S_pool);    
        results.tSqrd = (diff_mu - opt.testMu)*invCovMat*(diff_mu - opt.testMu)'; 
        tSqrdF = results.tSqrd*(df2/(df1*(df2 + 1))); %F-approx
        results.pVal = 1 - fcdf(tSqrdF, df1, df2);  % compute p-value
    else % same for now
        S_pool = ( (n1 - 1)*sCov1 + (n2 - 1)*sCov2)/(n1 + n2 - 2);
        diff_S_pool = S_pool*(1/n1 + 1/n2);      
        diff_mu = sMu1 - sMu2;
        
        invCovMat = inv(diff_S_pool);    
        results.tSqrd = (diff_mu - opt.testMu)*invCovMat*(diff_mu - opt.testMu)'; 
        tSqrdF = results.tSqrd*(df2/(df1*(df2 + 1))); %F-approx
        results.pVal = 1 - fcdf(tSqrdF, df1, df2);  % compute p-value
    end
end

