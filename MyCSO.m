function [ x,bestever,F ] = MyCSO( FitnessFcn,numberOfVariables,lb,ub,maxfe )
%MYCSO Summary of this function goes here
%   Detailed explanation goes here
% addpath(genpath(pwd));


% global initial_flag

%d: dimensionality
d = numberOfVariables;
%maxfe: maximal number of fitness evaluations
% maxfe = d*5000;
%runnum: the number of trial runs
% runnum = 1;

lu = single([lb;ub]);

% The frist six benchmark functions in  CEC'08  test suite.
% Function 7 is excluded because of the following error thrown by the test suite:
% 'Undefined function 'FastFractal' for input arguments of type 'char'.'

n = d;
initial_flag = 0;


%phi setting (the only parameter in CSO, please SET PROPERLY)
if(d >= 2000)
    phi = 0.2;
elseif(d >= 1000)
    phi = 0.1;
elseif(d >=500)
    phi = 0.05;
else
    phi = 0;
end
% population size setting
if(d >= 5000)
    population = 1500;
elseif(d >= 2000)
    population = 1000;
elseif(d >= 1000)
    population = 500;
elseif(d >= 100)
    population = 100;
else
    population = 50;
end


% initialization
% F = zeros(population,maxfe);
XRRmin = repmat(lu(1, :), population, 1);
XRRmax = repmat(lu(2, :), population, 1);
rand('seed', sum(100 * clock));
p = XRRmin + (XRRmax - XRRmin) .* rand(population, d);
clear('XRRmin','XRRmax');
fitness = FitnessFcn(p);
v = zeros(population,d,'single');
bestever = 1e200;

FES = population;
F(:,1) = [FES;fitness];
gen = 1;


tic;
% main loop
while(FES < maxfe)
    
    
    % generate random pairs
    rlist = randperm(population);
    rpairs = [rlist(1:ceil(population/2)); rlist(floor(population/2) + 1:population)]';
    
    % calculate the center position
    center = ones(ceil(population/2),1)*mean(p);
    
    % do pairwise competitions
    mask = (fitness(rpairs(:,1)) > fitness(rpairs(:,2)));
    losers = mask.*rpairs(:,1) + ~mask.*rpairs(:,2);
    winners = ~mask.*rpairs(:,1) + mask.*rpairs(:,2);
    
    
    %random matrix
    randco1 = rand(ceil(population/2), d,'single');
    randco2 = rand(ceil(population/2), d,'single');
    randco3 = rand(ceil(population/2), d,'single');
    
    % losers learn from winners
    v(losers,:) = randco1.*v(losers,:) ...,
        + randco2.*(p(winners,:) - p(losers,:)) ...,
        + phi*randco3.*(center - p(losers,:));
    p(losers,:) = p(losers,:) + v(losers,:);
    
    % boundary control
    for i = 1:ceil(population/2)
        p(losers(i),:) = max(p(losers(i),:), lu(1,:));
        p(losers(i),:) = min(p(losers(i),:), lu(2,:));
    end
    
    
    % fitness evaluation
    fitness(losers,:) = FitnessFcn(p(losers,:));
    bestever = min(bestever, min(fitness));
    
    
    %         fprintf('Best fitness: %e\n', bestever);
    FES = FES + ceil(population/2);
    F(:,gen+1) = [FES;fitness];
    gen = gen + 1;
end
[~,ind] = min(fitness);
x = p(ind,:);
% results(funcid, runnum) = bestever;
% fprintf('Run No.%d Done!\n', run);
% disp(['CPU time: ',num2str(toc)]);





end

