% SABVA JAY DILIPBHAI
% 202101224

tic
% Given H Matrix Here
H = load("H1.mat");

% Calcuation Of N, U, K
n = size(H.H, 2);
u = size(H.H, 1);
k = n - u;

% Degree Of Check Node and Variable Node
DC = sum(H.H(1, :));
DV = sum(H.H(:, 1));

% Valid Codeword
codeword = zeros(1, n);

% TannerGraph
CN = zeros(u, 2 * DC); % we take value corresponding to variable node
VN = zeros(n, DV + 1); % last row will have msg and all other will be connection

% TannerGraph Filler

% VN Array
for vn = 1: n
    idx = 1;
    for cn = 1: u
        if H.H(cn, vn) == 1 % if thery are connected
            VN(vn, idx) = cn;
            idx = idx + 1;
        end
    end
end

% CN Array
for cn = 1: u
    idx = 1;
    for vn = 1: n
        if H.H(cn, vn) == 1
            CN(cn, idx) = vn;
            idx = idx + 1;
        end
    end
end

% Monte Carlo Simulation
% Probability Array Because of algo-conv
p = 0: 0.02: 1;
tmax = 25; % max iteration allowed
Nsim = 1000; % No of simulation
Psuccess = zeros(size(p)); % store success for each p and all of Nsim
inx = 1; %index
algocon = zeros(length(p), tmax); % algo convergence 2d array for each probabilty there erasure per iteration

% Monte Carlo Loop
for Pe = p
    success = 0;
    erasure_sum = zeros(1, tmax); % for algo-conv
    for L = 1: Nsim

        % BEC Channel
        revcode = zeros(1, n); % Inital Recieved Codeword
        for i = 1: n
            r = rand;
            if r <= Pe
                revcode(i) = 9; % erasure
            else
                revcode(i) = codeword(i); % as it is
            end
        end

        % Magic Happens - Black Box
        % I am the heart of the This Entire Code Please Dont Touch Me
        % Erasures I am coming To Remove You

        decoded = ones(1, n) * 9; % decoded codeword initialized with erasures
        erasure_perit = zeros(1, tmax); % erasure per iterations

        % Loading Each Variable Node
        for vn = 1: n
            VN(vn, DV + 1) = revcode(vn);
        end

        % For t == 0 Tell Each CN that I am this 
        for cn = 1: u
            for vn = DC + 1: 2 * DC
                CN(cn, vn) = VN(CN(cn, vn - DC), DV + 1); % Don't be confuse just visulize how things work 
            end
        end

        % real magic start from here
        for t = 1: tmax

            % erasure_perit
            for vn = 1: n
                if VN(vn, DV + 1) == 9
                    erasure_perit(t) = erasure_perit(t) + 1; % Yo I found Erasure
                end
            end

            % Now CN To VN For which (SUM % 2)
            noOfCNCorrect = 0; % count how many CN we corrected or it is if == u then our job is done
            
            for cn = 1: u
                % Find How Many Erasure and Where it is
                erasure = find(CN(cn, DC + 1: 2 * DC) == 9);
                
                if size(erasure) == 0 % no erasure this CN has so this CN is correct
                    noOfCNCorrect = noOfCNCorrect + 1;
                elseif size(erasure) == 1 % if it has one erasure that can be solved
                    % which VN has erasure
                    vnind = CN(cn, erasure(1));
                    sum = 0; % sum mod 2
                    for vn = DC + 1: 2 * DC
                        if(CN(cn, vn) == 1 || CN(cn, vn) == 0)
                            sum = sum + CN(cn, vn);
                            sum = mod(sum, 2);
                        end
                    end
                    if sum == 0
                        VN(vnind, DV + 1) = 0;
                    else
                        VN(vnind, DV + 1) = 1;
                    end
                end

                % check whether all CN are corrected in this iteration then
                % exit this iteration loop
                if noOfCNCorrect == u
                    break;
                end
            end
            % VN To CN
            % Now it is time to Tell all CN That Which VN is updated
            for cn = 1: u
                for vn = DC + 1: 2 * DC
                    CN(cn, vn) = VN(CN(cn, vn - DC), DV + 1);
                end
            end
        end

        % after completing all iterations
        decoded(1, :) = VN(:,DV + 1);

        if decoded == codeword
            success = success + 1;
        else
            % Not Able To Decode This (I will try in next sim)
        end
        erasure_sum = erasure_sum + erasure_perit;
    end
    erasure_avg = erasure_sum ./ Nsim;
    algocon(inx, :) = erasure_avg;
    Psuccess(inx) = success / Nsim;
    inx = inx + 1;
end

% Plot 
% Probability Of Successful Decoding
plot(p,Psuccess);
xlabel('Error Probability Of BEC Channel');
ylabel('Probability Of Successful Decoding');
title('BEC Hard Decision - 202101224');
grid on;


% Algorithmatic Convergence 
figure;
hold on;
for i = 1: length(p)
    plot(algocon(i,:),'DisplayName',['p=' num2str(p(i))]);
end
legend('show');
xticks(0:10:25);
xlim([1,25]);
xlabel('Iteration Index');
ylabel('Average Number Of Erasure');
title('Algorithmatic Convergence (BEC Hard Decision) - 202101224');
toc
