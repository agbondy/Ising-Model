%% IsingModelTutorial
% Adrian Bondy, 2014-2015
% originally written at CSH Comp. Vision summer school, summer 2014

%% set up
% make sure subfolder 'tools' is in your path, 
% as well as the dichotomized gaussian m-files
% (get the latter at : http://bit.ly/pop-spike) and see Ref. 5
%addpath('/.../mackelab-pop_spike-28a1584cb671/util');
%addpath('/.../mackelab-pop_spike-28a1584cb671/dich_gauss');
clear cell spikes1 spikes2 spikes3

%% Intro
% The Ising model was originally developed in the field of statistical physics
% to describe how local interactions between nearby atoms could explain
% macro-scale properties of magnets. In neuroscience, it has been adapted to 
% model large-scale statistical properties of neuronal populations using
% only pairwise measurements.
%%
% The Ising model is a so-called "Maximum
% Entropy" model (MEM).  Generally, an MEM is the most "unstructured" distribution 
% over states of a system given certain constraints.
% The Ising Model is the special case where the system is binary and the only constraints
% are the pairwise interactions of the system.
%%
% What do we mean by "states" in the brain?
% Consider a set of three neurons each of which,
% in any given time bin, can either fire a spike or not.
% We can think of this as a binary physical system which can exist in one of 2^3 states:
% [0 0 0]
% [1 0 0]
%   ...
% [1 1 1] 
%
% The number of possible states grows exponentially with population size,
% and high-order interactions reflecting complex network connectivity can 
% create highly skewed distributions over these states. For instance, 
% say neuron 1 only spikes if it receives synchronous input from neurons 2 and 3.
%%
% So what can the Ising model tell us?  Well, it describes the most unstructed distribution of
% states if we confine ourselves to only
% those distribution that obey a set of 2nd order or pairwise statistics,
% i.e. a known covariance matrix. In some cases, the maximum entropy distribution
% may closely match the observed distribution of states,
% thus accurately describing a system with a paramter set that grows quadratically, rather
% than exponentially, with N.
% (See refs. [1-3] for experimental results).
%%
% On the flip side, we may be interested in seeing how badly the Ising model fails, 
% as a way of testing for the presence of higher order interactions that might otherwise be hard to reveal.  
% For instance, we'd expect that the Ising model might have some trouble explaining the
% example above, since knowing the correlations between neuron 1 and the
% other two neurons won't capture its dependence on synchronous input 
% (although, as we'll see, even these higher-order interactions can be
% explained relatively well with an Ising Model, which creates interpretive
% problems.)
%%
% NB: Ising models should be applied with care.  There are a number of
% situations where they can lead to seriously misleading conclusions, particularly when
% firing rates are low, small populations are used, and when the number of
% samples is limited. See refs. 4 and 6 for detailed treatments of these issues.

%%
% We'll consider three simulations involving a population of six neurons.
ncells=6;
nbins=10^5;

%% Example 1: The Independent Case  
% Let's generate six independent spike trains, with random mean firing rates.
rate = normrnd(0.2,0.05,1,ncells);
clear spikes
for cell=1:ncells
    spikes1(:,cell) = rand(nbins,1)<=rate(cell);
end

%%
% Because the spike trains are independent, we know that the probability of observing any state is just the product of the marginal probabilities.,
% (Remember, given independent x and y, P(x|y) = P(x).)  So an Ising model isn't going to help us much here.  The first-order
% statistics (empirical mean rate) are sufficient to fully describe the probability distribution over states.

%% 
% But let's show this empirically.
% First let's make a list of all the states the system can take.  
% 'dec2bin' is a handy function for this:

allstates=double(dec2bin(0:2^ncells-1)=='1');  

%% 
% And now let's calculate the empirical distribution over states (the observed frequencies).

decspikes=bin2dec(num2str(spikes1)); % this line takes a while to run
[counts,i] = Counts(decspikes);
freq=zeros(2^ncells,1);
freq(i+1)=counts;
freq=max(eps,freq);
freq=freq(:);
p = freq./sum(freq); % empirical distribution

%%
% Now we can use the observed means to estimate the
% probability of observing each state, without knowing the actual set of
% states we observed, just by multiplying the marginal distributions.

rates=mean(spikes1);
p1 = prod( bsxfun(@times,(1-allstates),(1-rates)) + bsxfun(@times,allstates,rates) , 2 ); % probability distribution of the independent model

%%
% Now let's plot the probabilities of each state under the independent
% model against the observed frequencies:

figure;
h=scatter(p(p>10^-7),p1(p>10^-7),25,'o','filled');
minp=min([p(p>10^-7) ; p1(p>10^-7)]);
set(gca,'yscale','log','xscale','log','xlim',[minp 1e0],'ylim',[minp 1e0]);
refline(1,0);
ylabel('Predicted Probabilities (Independent)');
xlabel('Empirical probabilities');
title('Example 1: Independent Model Performance');
legend(h,{'Independent'});

%% 
% This model does pretty well, as you can see.  In fact, it's the 1st order
% Maximum Entropy Model for our population, as it gives the most
% likely distribution of states constrained only on the empirical mean rate
% of each neuron (i.e. the 1st order statistics).

%%
% Question 1: Why doesn't this scatter plot lie EXACTLY along the identity line? 
%%
% Answer: Sampling error.

%% Example 2: 2nd Order Data
% Now what if we create a set of CORRELATED spike trains?
% There's a convenient and fast way to generate sample
% spike trains that are constrained by whatever pairwise correlations you
% want but which have *almost* no excess higher order correlations, using a
% dichotomized multivariate Gaussian distribution (see ref 5).
%%
% First, let's make a random covariance matrix, with a diagonal ridge.
success=false;
while ~success
    myCov=normrnd(0.05,0.1,ncells);
    myCov=triu(myCov);
    myCov=myCov*myCov';
    ridge=eye(ncells)*0.2;
    myCov=myCov+ridge;
    %%
    % now generate an example set of spike trains using the dichotomized
    % Gaussian method:
    try
        spikes2=sampleDichGauss01(rate,myCov,nbins);
        success=true;
    catch
        success=false;
    end
end

%% 
% compare the covariance matrix that we generated to the empirical covariance matrix of the spike trains we sampled to see this method works.
figure;
subplot(1,3,1);imagesc(Cov2Corr(myCov)); makeCorrelationImage;
title('Desired Correlation Matrix');
subplot(1,3,2);imagesc(corr(spikes2)); makeCorrelationImage;
title({'Observed Correlation Matrix' 'of Simulated Spike Trains'});
subplot(1,3,3);imagesc(corr(spikes2) - Cov2Corr(myCov)); makeCorrelationImage;
title('Difference');
set(gcf,'Position',[944 1180 1499 318]);
matchclim(gcf);
%% Independent Model
% First let's see how well the 1st order MEM does on our new set of spikes
decspikes=bin2dec(num2str(spikes2));
[counts,i] = Counts(decspikes);
freq=zeros(2^ncells,1);
freq(i+1)=counts;
freq=max(eps,freq);
freq=freq(:);
p = freq./sum(freq); % empirical distribution
rates=mean(spikes2);
p1 = prod( bsxfun(@times,(1-allstates),(1-rates)) + bsxfun(@times,allstates,rates) , 2 ); 
figure;
h=scatter(p(p>10^-7),p1(p>10^-7),25,'o','filled');
minp=min([p(p>10^-7) ; p1(p1>10^-7)]);
set(gca,'yscale','log','xscale','log','xlim',[minp 1e0],'ylim',[minp 1e0]);
refline(1,0);
ylabel('Predicted Probabilities (Independent)');
xlabel('Observed Probabilities');
title('Example 2: Independent Model Performance');
legend(h,'Independent');


%%
% The first order model does not do a good job at predicting the state
% probabilities we observe.  How could it, when our spike trains are not
% independent?
%% Pairwise Model
% OK. Now let's fit a second order (Ising) model to the same set of correlated spikes.
% In words, this means we'll use only the covariance matrix to constrain the set of probability distributions
% we'll consider and then we'll pick the one with the most entropy. 
% I wrote a function "IsingModel" which finds this distribution using gradient descent.  
% It can be used on any set of spike train data, not just in this
% tutorial (see the help).
stats = IsingModel(spikes2,'verbose',true);

%%
% A good way to see how well the Ising model did is to compare its
% performance to the performance of the independent model.  We can do this
% by comparing the entropies of both models relative to the entropy of the observed data.
% Remember, maximum entropy models try to
% generate the model with  the most entropy given the
% constraints.  So the observed entropy is necessarily a lower bound.
% We'll calculate performance as follows:

S_1 = stats.entropy(1); % entropy of the first order model
S_2 = stats.entropy(2); % entropy of the Ising model
S_observed = stats.entropy(3); % entropy of the observed distribution
performance = (S_1 - S_2) / ( S_1 - S_observed);

%%
% This metric (sometimes called the multi-information fraction) is equivalent to: D_KL(p_observed || p_independent) / D_KL(p_observed || p_pairwise) 
%%
% The multi-information fraction is bounded on [0,1] where 1 indicates perfect performance and
% 0 indicates no improvement over the independent model.  Check the value
% of 'performance' to see how we did.

%% 
% plot the results
figure;
h(1)=scatter(p(p>10^-7),p1(p>10^-7),25,'o','filled');
set(gca,'yscale','log','xscale','log');
refline(1,0);
ylabel('Predicted Probabilities (Independent)');
xlabel('Empirical probabilities');
title('Example 1: Independent Model Performance');
hold on;
h(2)=scatter(p(p>10^-7),stats.probabilities(p>10^-7,2),25,'bo','filled');
minp=min([p(p>10^-7) ; stats.probabilities(stats.probabilities>10^-7)]);
set(gca,'yscale','log','xscale','log','xlim',[minp 1e0],'ylim',[minp 1e0]);
ylabel('Predicted Probabilities (Pairwise)');
refline(1,0);
xlabel('Observed Probabilities');
title('Example 2 Ising Model Performance');
legend(h,{'Independent','Pairwise'});

%%
% OK, so the Ising model does a pretty darn good job at matching the
% empirically observed distribution of state probabilities.  But why
% doesn't it do perfectly?  After all, we know for a fact that there are no
% correlations beyond pairwise, since we designed the spikes this way.
% The answer is that sampling error effectively introduces higher order
% correlations. Unfortunately, the bias this introduces is always in 
% the direction of underestimating the performance of the Ising model, 
% since undersampling makes it
% look like our population has more structure (i.e. less entropy) than it
% really does.  (If this isn't intuitive, consider that some states are
% really, really unlikely given our covariance matrix.  With limited
% samples, we might never observe these states.  So our empirical
% distribution is 0 at this point.  To even get the math to work, we need to bump this up to some tiny number.  
% And even then, the entropy of this has gone waaay down.
% An Ising model will always assign p>0 to that point and thus deviate from a perfect description of the data).
% The only foolproof way around this is to get more samples.  But recent
% work has at least quantitatively investigated the issue of limited sampling for maximum entropy models (see ref 5) 
% and offers a correction term which may be a good solution in some contexts.


%% Example 3: Single Common Input
% In this example, we'll create spike trains that share a common input
% added to which is a bit of independent noise.  This will create a global
% fluctuation in the population.

clear spikes3
latentAmplitude=0.2;
for i=1:1000
    interval = (nbins/1000)*(i-1)+1:i*nbins/1000 ;               
    spikes3(interval,:) = rand(nbins/1000,ncells)>min(0.99,normrnd(0.45,latentAmplitude,1,1));
end

%% 
% Let's generate a raster plot of the spikes.
figure;
plotSpikeRaster(spikes3(1:1000,:)');
ylabel('Cell No');
xlabel('Time Bin');
title('Example 3: Spiking Population with Common Input');
set(gca,'ytick',1:ncells);

%% 
% You should be able to see the correlated up and down states clearly.

%% 
% Now let's see how the Ising model does at characterizing the observed distribution of states.

stats = IsingModel(spikes3,'verbose',true);
figure;
h(1)=scatter(p(p>10^-7),stats.probabilities(p>10^-7,1),25,'o','filled');hold on
h(2)=scatter(p(p>10^-7),stats.probabilities(p>10^-7,2),25,'ob','filled');
minp=min([p(p>10^-7) ; stats.probabilities(stats.probabilities>10^-7)]);
set(gca,'yscale','log','xscale','log','xlim',[minp 1e0],'ylim',[minp 1e0]);
ylabel('Predicted Probabilities');
refline(1,0);
xlabel('Observed Probabilities');
legend(h,{'Independent','Pairwise'});
S_1 = stats.entropy(1); % entropy of the first order model
S_2 = stats.entropy(2); % entropy of the Ising model
S_observed = stats.entropy(3); % entropy of the observed distribution
performance = (S_1 - S_2) / ( S_1 - S_observed);

%%
% The Ising model does reasonably well.  Perhaps counterintuitively,
% common input generates a correlation structure which can be described
% extremely well by second order interactions.  Pairwise statistics don't
% imply pairwise inputs.  This emphasizes that we are very separated from
% mechanism here, an important point to consider.

%%
% this example illustrates another disappointing feature of the Ising
% model.  It doesn't know anything about temporal dynamics, just state
% probabilities.  So if we sample from the Ising model fit to our sample
% population with a common input, we don't see the up and down states
% anymore.
sample_idx = randsample(2^ncells,nbins,true,stats.probabilities(:,2));
sample_spikes = logical(allstates(sample_idx,:));
figure;
plotSpikeRaster(sample_spikes(1:1000,:)');
ylabel('Cell No');
xlabel('Time Bin');
title('Example 3: Ising Model for Spiking Population with Common Input');
set(gca,'ytick',1:ncells);

%% 
%but see that these have the same covariance matrix
figure;
subplot(1,3,1);
imagesc(corr(spikes3));makeCorrelationImage;
title('"Real" simulated data with common input');
subplot(1,3,2);
imagesc(corr(sample_spikes));makeCorrelationImage;
title('Ising Model Fit');
subplot(1,3,3);
imagesc(corr(spikes3)-corr(sample_spikes));makeCorrelationImage;
title('Difference');
matchclim(gcf,'clim',[0 0.5]);
set(gcf,'Position',[944 1180 1499 318]);

%% Example 4: Higher Order Correlations
% In this example, we'll create spike trains with third-order structure:
% i.e. where the probabilities of spiking in one neuron depends on the product of the activity of two other neurons.
rate = 0.45;
clear spikes4 cells spikes_tmp
cells=1:ncells;
cells(2,:)=circshift(cells,[0,1]);
for cell=1:ncells
    spikes_tmp(:,cell) = rand(nbins,1)>normrnd(rate,0.2,1,1);
end
for cell=1:ncells
    spikes4(:,cell)=spikes_tmp(:,cells(1,cell))~=spikes_tmp(:,cells(2,cell));
end

%% 
% Let's generate a raster plot of the spikes.
figure;
plotSpikeRaster(spikes4(1:1000,:)');
set(gca,'ytick',1:ncells);
ylabel('Cell No');
xlabel('Time Bin');
title('Example 4: Spiking Population with Higher Order Correlations');

%% Pairwise Model
% Now let's see how the Ising model does at characterizing the observed distribution of states.

stats = IsingModel(spikes4,'verbose',true);
figure;
h(1)=scatter(p(p>10^-7),stats.probabilities(p>10^-7,1),25,'o','filled');hold on
h(2)=scatter(p(p>10^-7),stats.probabilities(p>10^-7,2),25,'ob','filled');
minp=min([p(p>10^-7) ; stats.probabilities(stats.probabilities>10^-7)]);
set(gca,'yscale','log','xscale','log','xlim',[minp 1e0],'ylim',[minp 1e0]);
ylabel('Predicted Probabilities');
refline(1,0);
legend(h,{'Independent','Pairwise'});
xlabel('Observed Probabilities');
S_1 = stats.entropy(1); % entropy of the first order model
S_2 = stats.entropy(2); % entropy of the Ising model
S_observed = stats.entropy(3); % entropy of the observed distribution
performance = (S_1 - S_2) / ( S_1 - S_observed);

%%
% Given the way we constructed the spikes and the scatter plot comparing the probability distributions,
% are you at all surprised at the performance value?  One counterintuitive
% observation: it's hard to generate third-order correlations without also
% generating some 2nd order correlations.

figure;imagesc(corr(spikes4));makeCorrelationImage;title('Example 4 Correlation Matrix');

%%
% Something to think about: the ratio of the number of parameters in an Ising model
% versus the number of parameters in the full model is n^2/2^n.  This ratio
% depends on n!  Which means you are GUARANTEED to get a better approximation
% with an Ising model applied to a smaller population.  Of course, what we really want is
% to be able to extrapolate our findings to the case of large n (i.e. the
% number of neurons in the brain).  It turns out this is a complicated
% affair!! (See ref 4).


%% Conclusions and Explorations
% Question: How should we interpret high performance from
% an Ising model in real data? 
%%
% Question: What insight does it give us (or not give us) about
% patterns of correlated activity in neurons?
%%
% Try rerunning this simulation several time and seeing if the results are
% different just due to sampling variation.  Also try running it with
% different numbers of cells and seeing what happens.

%% References:
% # Schneidman, E., Berry, M. J., Segev, R., & Bialek, W. (2006). 
% Weak pairwise correlations imply strongly correlated network states in a neural population. 
% Nature, 440(7087), 1007?12. doi:10.1038/nature04701
% # Shlens, J., Field, G. D., Gauthier, J. L., Grivich, M. I., Petrusca, D., Sher, A., ? Chichilnisky, E. J. (2006).
% The structure of multi-neuron firing patterns in primate retina. 
% The Journal of Neuroscience, 26(32), 8254?66. doi:10.1523/JNEUROSCI.1282-06.2006
% # Tang, A., Jackson, D., Hobbs, J., Chen, W., Smith, J. L., Patel, H., & Beggs, J. M. (2008).
% A maximum entropy model applied to spatial and temporal correlations from cortical networks in vitro.
% The Journal of Neuroscience, 28(2), 505?18. doi:10.1523/JNEUROSCI.3359-07.2008
% # Roudi, Y., Nirenberg, S., & Latham, P. E. (2009).
% Pairwise maximum entropy models for studying large biological systems: when they can work and when they can't.
% PLoS Computational Biology, 5(5), e1000380. doi:10.1371/journal.pcbi.1000380
% # Macke, J. H., Berens, P., Ecker, A. S., Tolias, A. S., & Bethge, M. (2009).
% Generating spike trains with specified correlation coefficients. 
% Neural Computation, 21(2), 397?423. doi:10.1162/neco.2008.02-08-713
% # Macke, J., Murray, I., & Latham, P. (2013).
% Estimation Bias in Maximum Entropy Models. 
% Entropy, 3109?3129. doi:10.3390/e15083209