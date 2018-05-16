clear;
clc;

%% Import Real Data
dat = importdata("luca_alldata.txt");
label = importdata("luca_ally.txt");
X = dat.data';
Y = label;
[n,p] = size(X);

X_train = X(1:round(2*n/3),:);
X_test = X(round(2*n/3)+1:end,:);
Y_train = Y(1:round(2*n/3),:);
Y_test = Y(round(2*n/3)+1:end,:);
test_size = size(X_test,1);

%%  Logistic + SCAD %%
col = size(X_train,2);
row = size(X_train,1);
beta = zeros(col,1);

%  Calculating the beta_zero  %
temp = sum(Y_train)/row;
beta_zero = log(temp/(1-temp));

% Inputting X, Y, beta_int and lambda %
beta_int = [beta_zero;beta];
x0 = ones(row,1);
X = [x0,X_train];
Y = Y_train;

% Setting lambda
lambda_max = norm(X'*Y,'inf'); % according to the https://github.com/yangziyi1990/SparseGDLibrary.git
lambda_min = lambda_max * 0.001;
m = 10;
for i = 1:m
    Lambda1(i) = lambda_max * (lambda_min/lambda_max)^(i/m);
    lambda = Lambda1(i);
    beta = Logistic_SCAD_func(X,Y,beta_int,lambda);   
    beta_path(:,i) = beta;
    fprintf('iteration times:%d\n',i);
end

%% Result (without CV)
for i = 1:m
    beta_i = beta_path(:,i);
    beta_zero = beta_i(1);
    beta = beta_i(2:end);
    l = beta_zero + X_test * beta;
    prob = exp(l)./(1 + exp(l)); 
    
    for j = 1:test_size
        if prob(j) > 0.5
            Y_validation(j) = 1;
        else
            Y_validation(j) = 0;
        end
    end
    
    error=Y_validation'-Y_test;
    error_number_testing=length(nonzeros(error));
    beta_non_zero=length(nonzeros(beta_i));
    
    fprintf('beta path: %d\n', i);
    fprintf('The error number of testing data: %d\n', error_number_testing);
    fprintf('The number of nonzero beta: %d\n', beta_non_zero);  
    
    [accurancy,sensitivity,specificity]=performance(Y_test,Y_validation');
    fprintf('The accurancy of testing data (SCAD): %f\n' ,accurancy);
    fprintf('The sensitivity of testing data (SCAD): %f\n' ,sensitivity);
    fprintf('The specificity of testing data (SCAD): %f\n' ,specificity);
    fprintf('\n');
end

%% Cross Validation
% [Opt,Mse]=CV_SCAD_logistic(X,Y,Lambda1);
% beta_opt=beta_path(:,Opt);

%% Result (with CV)
% fprintf('Beta: %f\n' ,beta_opt(1:11));
% 
% beta_zero=beta_opt(1); 
% beta=beta_opt(2:end); 
% l = beta_zero + X_test * beta;
% prob=exp(l)./(1 + exp(l)); 
% for i=1:test_size
%     if prob(i)>0.5
%         Y_validation(i)=1;
%     else
%         Y_validation(i)=0;
%     end
% end
% 
% error=Y_validation'-Y_test;
% error_number_testing=length(nonzeros(error));
% beta_non_zero=length(nonzeros(beta_opt));
% fprintf('The error number of testing data: %d/n', error_number_testing);
% fprintf('The number of nonzero beta: %d/n', beta_non_zero);

%% Performance
% [accurancy,sensitivity,specificity]=performance(Y_test,Y_validation');
% fprintf('The accurancy of testing data (SCAD): %f\n' ,accurancy);
% fprintf('The sensitivity of testing data (SCAD): %f\n' ,sensitivity);
% fprintf('The specificity of testing data (SCAD): %f\n' ,specificity);
