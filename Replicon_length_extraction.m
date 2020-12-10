load ../mat' files'/dnaseq_sonication_async_smooth.mat
loglambda=estimateLambda(async_reads);
n = length(async_strains);
L = zeros(n,1);
for i = 1:n
    L(i) = power(10,mean(loglambda(i)));
    
end
Asynch_Lambda = [async_strains,num2cell(L)];
save ../mat' files'/Asynch_Lambda Asynch_Lambda async_strains async_reads;
