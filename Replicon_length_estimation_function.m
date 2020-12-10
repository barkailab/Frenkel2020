function [ loglambda1, projected_data] = estimateLambda(data)
%ESTIMATELAMBDA Summary of this function goes here

%after lambda_prediction_parameters is set - the next line  can be deleted -
%saveParameters();




load ../mat' files'/lambda_prediction_parameters.mat
c1=(M1\data)';
projected_data = M1(:,1:2)*c1(:,1:2)';
loglambda1=alpha_to_lambda(1)*(c1(:,2)./c1(:,1))+alpha_to_lambda(2);
end

function [] = saveParameters()
load ../mat' files'/simulations.mat
[u,s,v] = svd(simulated_data,0);
M1=u*s;
alpha_to_lambda=polyfit(v(:,2)./v(:,1),log10(lambda(sim_type))',1);
save lambda_prediction_parameters M1 alpha_to_lambda X
end
