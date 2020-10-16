clear
close all

%% [IMPORTANT] For this code to run, please download shapecomm.de toolbox from:
% https://www.shapecomm.de/downloads/shapecomm-webdm-v0.4.zip

M = 4; 
L = 20;
r = 1.5;

% Reference CCDM
dm = shapecomm.webdm(M, L, r);

% Input bits
bits_in = randi([0,1], dm.k,1);
% bits_in = zeros(dm.k, 1); % test
% bits_in = ones(dm.k, 1); % test

% bits_in = [   1     0     1     0     1     0     1     1     1     0   ];

% Distribution of the symbols  
dm.type_seq.'

% Verify encoder output with shapecomm's system
verify = dm.encode(bits_in)


%% Encoder 
L = L;
k = dm.k;
candidateList = dm.type_seq;

% p_a of each symbol in the alphabet 
p_a = candidateList/sum(candidateList);

% cmf of p_a
p = cumsum(p_a);

% Initialise limits for the interval 
ll = 0; % lower limit
uu = 1; % upper limit

% Reset 
n = 0;

% Pre-allocate symbols out 
symbols_out = zeros(L,1);

% Plot out the encoder operation
% 
z = 0:k+L;
u_hist = NaN*zeros(k+L+1,1);            u_hist(1) = ll;
v_hist = NaN*zeros(k+L+1,1);            v_hist(1) = uu;
p_hist = NaN*zeros(k+L+1,length(p));    p_hist(1,:) = p_a;

figure
hold on
h_p = area(z,p_hist);
h_u = plot(z,u_hist,'-+');
h_v = plot(z,v_hist,'-+');
legend('cmf(a1)', 'cmf(a2)', 'cmf(a3)' , 'cmf(a4)', 'lower limit', 'upper limit')


m = 0;

% Encoder loop
while n<L
    if m<k
        % Load next bit to encode
        m = m+1;
        b = bits_in(m);
    else
        b = 0; % append .5s if ran out of bits to encode
    end
    
    % Update limits for the loaded bit
    q = uu - ll;
    ll = ll + b*q/2; %lower
    uu = uu - (1-b)*q/2; % upper 
    
    % Update figure
    u_hist(1+m+n) = ll;      h_u.YData = u_hist;
    v_hist(1+m+n) = uu;      h_v.YData = v_hist;

    p_hist(1+m+n,:) = p_a;  for i = 1:length(p);h_p(i).YData = p_hist(:,i);end
    drawnow
    
    % Ppick the first candidate for which [p_lower,p_upper) \in [u,v)
    % first pick the first candidate for which p>=u
	candidate = find([0;p(1:end-1)]<=ll,1,'last');
    
    %check if p+r < v
    while p(candidate)>=uu
        % we have found a candidate completely within [u,v)
        
        %outpout symbol
        n = n+1; % increment symbol output by 1
        symbols_out(n) = candidate;
        
        % Catch any errors by comparing our input to shapecomm's
        if verify(n) ~= candidate
%             keyboard
        end
        
        %adjust and rescale u and v
        if candidate>1
            ll = ll - p(candidate-1);
            uu = uu - p(candidate-1);
        end
        ll = ll/p_a(candidate);
        uu = uu/p_a(candidate);
        
        % Draw one symbol from the candidate pool
        candidateList(candidate) = candidateList(candidate) - 1;
        assert(all(candidateList>=0), 'Picked a symbol from an empty pool');
        
        % Recalculate p_a and p (i.e. cmf) based on the updated
        % pool of the candidates
        p_a = candidateList/sum(candidateList);
        p = cumsum(p_a);
        
        % Update candidate
        candidate = find([0;p(1:end-1)]<=ll,1,'last');
        
        %update figure
        u_hist(1+m+n) = ll;      h_u.YData = u_hist;
        v_hist(1+m+n) = uu;      h_v.YData = v_hist;
        p_hist(1+m+n,:) = p_a;  for i = 1:length(p);h_p(i).YData = p_hist(:,i);end
        drawnow
    end
end


symbols_out


%% Decoder 

symbols_in = symbols_out;
candidateList = dm.type_seq;

% Initialise lower and upper limits for the interval
ll = 0;
uu = 1;

% p_a of each symbol in the alphabet 
p_a = candidateList/sum(candidateList);

% cmf of p_a
p = [0; cumsum(p_a)];

% bit index
m = 0;

% Preallocate bits
bits_out = zeros(k,1);

% Decoder loop
for n = 1:L
    if m>=k
        warn('ended early')
        break
    end
    candidate = symbols_in(n);
    
    q = uu - ll;
    uu = ll + q*p(candidate+1);
    ll = ll + q*p(candidate);
    
    % Update candidate list and recalculate p_a and p accordingly
    candidateList(candidate) = candidateList(candidate)-1;
    p_a = candidateList/sum(candidateList);
    p = [0; cumsum(p_a)];
    
    % Recalculate intervals
    while (uu<0.5||ll>0.5)&&(m<k)
        m = m+1;
        fprintf('%d ',bits_in(m))
        if uu<0.5
            bits_out(m) = 0;
            ll = 2*ll;
            uu = 2*uu;
        elseif ll>=0.5 % upper interval is [.5, 1)
            bits_out(m) = 1;
            ll = 2*ll - 1;
            uu = 2*uu - 1;
        else
            error('No!')
        end
        fprintf('%d\n',bits_out(m))
    end
end


