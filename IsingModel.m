function stats = IsingModel(varargin)
% IsingModel Binary Ising Model 
%   stats = IsingModel(spikes) fits a binary Ising model to the numeric array spikes.
%       spikes should be m x n, where m is the number of bins and n is the
%       number of dimensions (neurons), and should contain only 1s and 0s.
%       stats is a structure with fields:
%           'entropy' - vector of length 3 containing:
%                       { 'S1' (1st order MaxEnt Model entropy),
%                         'S2' (2nd order MaxEnt Model / Ising entropy),
%                         'Sn' entropy of the empirical distribution }
%           'info2' = S1 - S2 
%           'infoN' = S1 - Sn
%           'performance' = info2/infoN (Ising model multi-information fraction)
%           'J_lagrange' - Ising model fit param matrix
%           'fitstats' - output of Matlab's optimization routine
%           'logL' - vector of length 4 containing:
%                       { log likelihood of uniform model,
%                         log likelihood of 1st order MEM,
%                         log likelihood of 2nd order / Ising,
%                         entropy of empirical distribution = -Sn }
%           'probabilities' - 1 x 3 vector of probability distributions:
%                       { 1st order MaxEnt Model,
%                         2nd order MaxEnt Model / Ising,
%                         empirical distribution }
%           'nmissing' - number of states for which no samples were observed
%           'Nc' - critical number of neurons needed to avoid the
%                  perturbative regime (ncells>>Nc or you are in trouble, see ref. 1)
%           'ncells' - number of cells
%           'correlations' - observed correlation matrix
%           'nsamples' - number of samples
%
%   stats = IsingModel() generates simulated spike trains from 10
%       uncorrelated neurons and then fits an Ising model.
%
%   stats = IsingModel() generates simulated spike trains from 10
%       WEAKLY CORRELATED neurons and then fits an Ising model.
%       Correlated spike trains are generated using the dichotomized
%       Gaussian method (see ref. 2.  For this to work, you must have the 
%       dichotomized Gaussian code in your path. Get it at: http://bit.ly/pop-spike.
%       You need folders 'util' and 'dich_gauss').

%   stats = IsingModel(...,'verbose') plots the log likelihood surface as its
%       minimum is searched.  NB: Slow!

% Adrian Bondy, 2015
% (initially written for CSH Computational Vision summer course, summer 2014)
    
    global state
    p=inputParser;
    p.KeepUnmatched=true;
    p.addOptional('spikes',[],@(x)validateattributes(x,{'numeric','logical'},{'ndims',2,'binary'}));
    if ismethod(p,'addParameter')
        method=@addParameter;
    else
        method=@addParamValue;
    end
    method(p,'ncells',10,@(x)validateattributes(x,{'numeric'},{'scalar','>',1}));
    method(p,'nbins',10^5,@(x)validateattributes(x,{'numeric'},{'scalar'}));    
    method(p,'verbose',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    method(p,'display',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    method(p,'performance',false,@(x)validateattributes(x,{'logical'},{'scalar'}));
    method(p,'corrspikes',false,@(x)validateattributes(x,{'logical'},{'scalar'}));    
    p.parse(varargin{:});
    state=p.Results;
    state.count=0;
    state=rmfield(state,'spikes');
    spikes=p.Results.spikes;
    if ~ismember('spikes',p.UsingDefaults) % user input spikes
        [state.nbins,state.ncells]=size(spikes);
        if state.ncells<2
            error('Dimension 2 of spikes must be at least 2.');
        end
    else
        % simulate spikes
        if state.corrspikes 
            %% simulate correlated spike trains
            goodcov=0;                
            while ~goodcov
                try
                    myCov=normrnd(0.05,0.1,state.ncells,state.ncells);
                    mnvar=normrnd(0.5,0.1,1,state.ncells);               
                    myCov(logical(eye(state.ncells)))=mnvar;
                    myCov=nearestSPD(myCov);
                    spikes=sampleDichGauss01(mnvar,myCov,state.nbins);
                    goodcov=1;
                catch
                end
            end
        else
            %% simulated independent spikes
            fprintf('Simulating w/ %d independent spike trains, firing prob. is 0.5.\n',state.ncells);            
            spikes=double(dec2bin(randsample(2^state.ncells,state.nbins,true)-1)=='1'); % generate spike train                        
        end
    end
    %% list all 2^ncells words and calculate empirical distribution
    state.allW=double(dec2bin(0:2^state.ncells-1)=='1');    
    nstates=size(state.allW,1);
    [~,spkid] = ismember(spikes,state.allW,'rows');
    [counts,i] = Counts(spkid);
    state.freq=zeros(2^state.ncells,1);
    state.freq(i)=counts;
    state.freq=max(eps,state.freq);
    state.freq=state.freq(:);
    state.freq = state.freq./sum(state.freq); % empirical distribution
    nmissing=sum(state.freq==eps);
    if nmissing
        fprintf('%d of 2^%d states never observed.  (More data may be needed to estimate the state probabilities).',nmissing,state.ncells);
    end
    logL = sum(state.freq.*log(state.freq)); % entropy of the empirical distribution   
    Sn = -logL;
    %% the independent model
    empirical_p=mean(spikes);
    p1 = prod( bsxfun(@times,(1-state.allW),(1-empirical_p)) + bsxfun(@times,state.allW,empirical_p) , 2 ); % probability distribution of the independent model
    S1 = -sum(p1.*log(p1)); % entropy of the independent model
    logL1 = sum(state.freq.*log(p1)); % log likelihood of the data under the independent model
    %% fit the Ising model
    J_0=cov(spikes);
    %J_0=eye(state.ncells); % initializing w/ an identity matrix seems to work best
    fprintf('Fitting Ising model.\n');tic;
    options=optimset('TolFun',10^-7,'Display','off','TolX',10^-7,'MaxFunEvals',10^7,'MaxIter',10^7);    
    if any(strncmpi(toolboxes,'Opt',3))
        optimfun=@fminunc;
        warning('off','optim:fminunc:SwitchingMethod');
    else
        optimfun=@fminsearch;
    end
    if state.display
        figure;hold on;
    end    
    for i=nstates:-1:1 
        state.statemats(:,:,i)=state.allW(i,:)'*state.allW(i,:);
    end
    state.zrdiag=repmat(~eye(state.ncells),[1,1,nstates]);        
    [J,~,~,output] = optimfun(@IsingObjective,J_0(:),options); % find the J that maximizes the log likelihood of the Ising model
    %% evaluate the Ising model
    [neglogL2,p2] = IsingObjective(J(:)); 
    logL2=-neglogL2;
    J=reshape(J,sqrt(numel(J)),sqrt(numel(J)));
    dummy_logL = sum ( state.freq * log(1/nstates) ); % log likelihood of uniform model
    S2 = -sum(p2.*log(p2)); % Ising model entropy
    %% calculate model performance
    In=S1-Sn;
    I2=S1-S2;
    performance=I2/In;
    Nc=1./mean(mean(spikes,2));
    correlations=corr(spikes);
    stats=struct('performance',performance,'info2',I2,...
        'infoN',In,'entropy',[S1 S2 Sn],'J_lagrange',J,'fitstats',output,...
        'logL',[dummy_logL logL1 logL2 logL],'probabilities',[p1(:) p2(:) state.freq(:)],...
        'nmissing',nmissing,'Nc',Nc,'ncells',state.ncells,'correlations',correlations,'nsamples',state.nbins);      
    if state.verbose
        fprintf('Took %s seconds to evaluate function %d times.\n',myround(toc,3),output.funcCount);
        fprintf('Entropy of empirical distribution is %g.\n',myround(Sn,3));
        fprintf('Entropy of maximum entropy Ising model is %g.\n',myround(S2,3));  
        % entropy figure
        figure;bar([S1 S2 Sn]);
        set(gca,'ylim',[Sn-(S2-Sn) S1+(S1-S2)],'xticklabel',{'1st order model' 'Ising model' 'Observed'});
        box off;
        ylabel('Entropy');
        title({['performance of Ising model is ',num2str(round(performance*1000)/10),'%']});    
        fprintf('Performance of the 2nd order model is %g%%\n.',round(performance*1000)/10);
        fprintf('Nc is %s.\n',num2str(myround(Nc,2)));
        % Log Likelihood Figure (turned this off for now because it really
        % doesn't add anything beyond the entropy figure -- entropy and log
        % likelihood are extremely closely related)
        %logL_percentages = ( stats.logL - stats.logL(1) ) ;
        %logL_percentages=logL_percentages / logL_percentages(end);
        %bar(logL_percentages(2:end));
        %set(gca,'xticklabel',{'1st order model' 'Ising model' 'Observed'});       
        %ylabel('Relative Log Likelihood');
        %box off
    end
end

function [neglogL,p] = IsingObjective(J)
    % calculates the negative log-likelihood of an Ising model defined by parameter matrix J
    % for data in the global variable 'state'
    global state
    J=reshape(J,state.ncells,state.ncells);
    diagJ=diag(J)';
    energysum=squeeze(sum(sum(bsxfun(@times,state.statemats,J) .* state.zrdiag))); 
    E = sum(bsxfun(@times,diagJ,state.allW),2) + energysum/2;  % Ising model state energies
    z = sum ( exp ( E ) ) ; % Ising model partition function
    p =  ( 1/z .* exp ( E ) ) ; % state probabilities
    neglogL = - sum ( state.freq .* log(p) ); % negative log likelihood
    if state.display
        state.count=state.count+1;
        scatter(state.count,neglogL,'.k');
        drawnow;
    end  
end

%% References
% 1. Roudi, Y., Nirenberg, S., & Latham, P. E. (2009).
% Pairwise maximum entropy models for studying large biological systems: when they can work and when they can't.
% PLoS Computational Biology, 5(5), e1000380. doi:10.1371/journal.pcbi.1000380

% 2. Macke, J. H., Berens, P., Ecker, A. S., Tolias, A. S., & Bethge, M. (2009).
% Generating spike trains with specified correlation coefficients. 
% Neural Computation, 21(2), 397?423. doi:10.1162/neco.2008.02-08-713
