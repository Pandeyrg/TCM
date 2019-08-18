
clear;

N_WORLD_FEATURES = 5;
N_ITEMS = 10;
ENCODING_TIME = 500;
TEST_TIME = 20;

% continuing with the optimal scheduling structure.
schedule = zeros(N_ITEMS , 2);
for i=6:10
    schedule(i,1) = 490 + i;
    schedule(i,2) = i;
end

x = 99;
schedule(1,1)=1;
schedule(1,2)=1;
for i = 1:5
    j = 6-i;
    schedule(j,1)= 496 - x*i ;
    schedule(j,2)= j ;
end

%disp(median(diff(schedule(: , 1))));

schedule_load = ENCODING_TIME/median(diff(schedule(:,1)));                  % variable important for parts 2,3 of assignment

encoding = zeros( N_ITEMS, N_WORLD_FEATURES+1 );

world_m = [1 2 1 2 3];              % can generate randomly for yourself
world_var = 1;
%delta = 0.05;                       
beta_param = 0.001;                 
m = 1;
s_mean = 0.1;
l_mean = 10;
% simulating encoding 
X = zeros(ENCODING_TIME , 1); % to store the delta and hence learn its model later.

for time = 1:ENCODING_TIME
    a = rand;
    if a< 0.5
        delta  = normrnd(s_mean , world_var);
    else
        delta = normrnd(l_mean , world_var);
    end
     X(time , 1 ) = delta;    
    
    world_m = world_m + delta;
    world = normrnd(world_m, world_var);
    % any item I want to encode in memory, I encode in association with the
    % state of the world at that time.
    if(m<(N_ITEMS+1))
        if(time==schedule(m,1))
            encoding(m,:) =    [world ,m ];                                             % encode into the encoding vector
            m =  m + 1;
        end                                %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               end;  
    end
end

%%%%%%% E M Algorithm %%%%%%%%%%%
X = transpose(X);

[lab , model , llh] = mixGaussEm(X,2);
t =length(model.mu);
%disp(model)

if t == 2  
    mu1 = model.mu(1,1);
    mu2 = model.mu(1 , 2);
    sigm1 =model.Sigma(1,1,1);
    sigm2 = model.Sigma(1,1,2);
    p = model.w(1,1);  % probability of using the first gaussian for sampling as learned by EM alogrithm.
else
    % incase EM algorithm thinks one distribution can alone make for the
    % GMM  which might occur depending on the delta's  used in the past
    % during encoding.
    mu = model.mu(1,1);
    sigm = model.Sigma(1,1,1);
end


%disp(mu1);
%disp(mu2);

%disp(sigm1)
%disp(sigm2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while(time<ENCODING_TIME+TEST_TIME)
    if t ==2
        a = rand;
        if a < p
            delta = normrnd(  mu1 , sigm1*sigm1 );
        else
            delta = normrnd( mu2 , sigm2*sigm2 );
        end
    else
        delta = normrnd(mu , sigm*sigm);
    end
    
    world_m = world_m + delta ;
    cue = normrnd(world_m, world_var);
                                                                                    % model world evolution                                                                                      
    for m = 1:N_ITEMS
        soa(m) = dot(cue,encoding(m, 1 : N_WORLD_FEATURES))/dot(cue,cue); % finding association strengths
    end
    
                                                                            
    try
        out(time-ENCODING_TIME+1) = find(drawFromADist(soa));
    catch ME
        out(time - ENCODING_TIME+1) = 0;
    end
    
    time = time + 1;       
end
out(ENCODING_TIME + TEST_TIME) = 0;
success = length(unique(out)) - 1 ;  % success is number of unique retrievals
disp("When is delta is from GMM whose parameters is unknown to the retrieving agents");
disp(success);


