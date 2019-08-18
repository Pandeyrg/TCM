clear;

N_WORLD_FEATURES = 5;
N_ITEMS = 10;
ENCODING_TIME = 500;
TEST_TIME = 20;


%%%%% delta is sampled from a mixture of two gaussians from here
%%%%% onwards.

%%% generating the most accuarate scheduling structure with respect to retrieval rate.
schedule = zeros(N_ITEMS , 2);
for i = 1:10
    schedule(i,1)= 490 + i ;
    schedule(i,2)= i ;
end

%disp(median(diff(schedule(: , 1))));

schedule_load = ENCODING_TIME/median(diff(schedule(:,1)));                  % variable important for parts 2,3 of assignment

encoding = zeros( N_ITEMS, N_WORLD_FEATURES+1 );

world_m = [1 2 1 2 3];              
world_var = 1;
%delta = 0.05;                       
beta_param = 0.001;                 
m = 1;
s_mean = 0.1;
l_mean = 10;
% simulating encoding 

for time = 1:ENCODING_TIME
    a = rand;
    if a > 0.5
        delta = normrnd(s_mean , world_var);
    else
        delta = normrnd(l_mean, world_var);
    end
    
      
    world_m = world_m + delta;
    world = normrnd(world_m, world_var);
    % any item I want to encode in memory, I encode in association with the
    % state of the world at that time.
    if(m<(N_ITEMS+1))
        if(time==schedule(m,1))
            encoding(m,:) =    [world , m ];                                             % encode into the encoding vector
            m =  m + 1;
        end                                %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               end;  
    end
end


while(time<ENCODING_TIME+TEST_TIME)
% the state of the world is the retrieval cue
    a = rand;
    if a > 0.5
        delta = normrnd(s_mean, world_var);
    else
        delta = normrnd(l_mean, world_var);
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
success = length(unique(out)) - 1 ;  % success is number of distinct retrievals
disp("Success Rate = ");
disp(success);
disp("ENCODING LOAD = ")
disp(schedule_load);

disp(" This shows a trade off between retrieval rate and encoding load ");
disp(" This scheduling structure has a very good retrieval rate, but alongwith that we have high encoding load also which needs to be minimised");



