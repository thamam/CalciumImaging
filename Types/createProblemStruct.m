function opts =  createProblemStruct()

opts = struct('tmax',4,'x_params',[8,2,2],'N_samp',50,'N_trial',5 ,...
    'rateOffset',2, 'cosntx' ,false,'overDisp', false,...
    'threshVal' ,0,'min_dt',2e-3 ,'pmax', Inf,'plotOpt', false,...
    'n_params' ,[0.5, 0.5, 1] );

end