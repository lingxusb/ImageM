function [ fpara, prob ] = Fitequi( sinten,snum )
%FITEQUI fit a equidistant distribution for the input [x, density]
%   sinten: is the spot intensity
%   snum: is the hist result of the spot intensity
%   fpara(1): density of a single transcript
%   fpara(2): the variation for each peak
%   prob: percentage for different transcripts

%% cost function for two different distributions
function r = cost(fpara, sinten, snum)
    %get matrix of activation
    for j = 1:floor(max(sinten)/fpara(1))
        sample = normrnd(fpara(1)*j,fpara(2),1,10000);
        m(j,:) = hist(sample,sinten);
    end
    sample_num =  abs(mldivide(m',snum'));
    y = SampleDis(sinten, fpara,floor(sample_num*100000/sum(sample_num)));
    r = norm(snum/sum(snum)-y);
end

%% optimization
optims=optimset('fminsearch');
optims.MaxIter = 1000;
optims.MaxFunEvals = 1000;
fval_prev = 1000;
[fpara,prob] = fminsearch(@(fpara) cost(fpara, sinten, snum),[5000 1000])

%% plot fitting results
for j = 1:floor(max(sinten)/fpara(1))
   sample = normrnd(fpara(1)*j,fpara(2),1,10000);
   m(j,:) = hist(sample,sinten);
end
sn =  abs(mldivide(m',snum'));
y_pre = SampleDis(sinten, fpara,floor(sn*100000/sum(sn)));
plot(sinten,snum/sum(snum));
hold on
plot(sinten,y_pre,'r');
end

%% get sampled distribution
%   sinten: x-axis
%   fpara: parameter that constrain the distribution
function y_pre = SampleDis(sinten, fpara, sample_num)
    peak_num = floor(max(sinten)/fpara(1));
    sample = [];
    for j = 1:peak_num
        sample = [sample, normrnd(fpara(1)*j,fpara(2),1,sample_num(j))];
    end
    y_pre = hist(sample,sinten)/sum(sample_num);
end



