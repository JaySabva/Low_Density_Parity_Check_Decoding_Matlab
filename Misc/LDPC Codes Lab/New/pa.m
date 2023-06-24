% SABVA JAY DILIPBHAI
% 202101224

% Soft Decision Final
clear all
tic
% Given H Matrix Here
H = load("Hmatrix2.mat");

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
VN = zeros(n, 2 * DV + 1); % last row will have msg and all other will be connection

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
p = [0.3, 0.4, 0.5, 0.51, 0.52];
tmax = 50; % max iteration allowed
Nsim = 100; % No of simulation
Psuccess = zeros(size(p)); % store success for each p and all of Nsim
inx = 1; %index
algocon = zeros(length(p), tmax); % algo convergence 2d array for each probabilty there erasure per iteration

% Monte Carlo Loop
for Pe = p
    disp(Pe);
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
            VN(vn, 2 * DV + 1) = revcode(vn);
        end

        % For t == 0 Tell Each CN that I am this and this are the
        % probabilites
        for cn = 1: u
             for vn = 1: DC
                 if VN(CN(cn, vn), 2 * DV + 1) == 0
                     CN(cn, vn + DC) = 0;
                 elseif VN(CN(cn, vn), 2 * DV + 1) == 1
                     CN(cn, vn + DC) = 1;
                 elseif VN(CN(cn, vn), 2 * DV + 1) == 9
                     CN(cn, vn + DC) = 0.5;
                 end
             end
        end 

        decodedPrev = VN(:, 2 * DV + 1);
        decodedCur = ones(1, n) * 9;

        % real magic start from here
        for t = 1: tmax

            % erasure_perit
            for vn = 1: n
                if decodedPrev(n) == 9
                    erasure_perit(t) = erasure_perit(t) + 1; % Yo I found Erasure
                end
            end

            % Now CN To VN For which (SUM % 2)
            % Bernard Equation
            for cn = 1:u
               for vn = 1:DC
                   % Count Probability Of * All Except This VN
                   prob0 = 1;
                   for vnvalue = DC + 1: 2 * DC
                       if vn + DC ~= vnvalue
                           prob0 = prob0 * (1 - 2 * CN(cn, vnvalue));
                       end
                   end
                   prob0 = (1 + prob0) / 2;
                   % Update That CN - > VN in hard decision there is only
                   % one to update so no need to check which cn is sending
                   % this but in soft we need to store which cn is sending
                   % this value
                   for vn1 = 1: DV
                       if VN(CN(cn,vn), vn1) == cn
                           VN(CN(cn, vn), vn1 + DV) = 1 - prob0;
                           break;
                       end
                   end
               end
            end

            % VN To CN
            % Now it is time to Tell all CN That Which VN is updated
            % See Slide's Of The Professor Yash Vasavada's It is Briliant. 
            for vn = 1:n
                for cn = 1:DV
                   if(VN(vn, 2 * DV + 1) == 0)
                        prob = 0;
                    elseif (VN(vn, 2 * DV + 1) == 1)
                        prob = 1;
                    elseif (VN(vn, 2 * DV + 1) == 9)
                        prob = 0.5;
                    end
                    prob0 = 1 - prob;
                    prob1 = prob;
                    for cn1 = DV + 1: 2 * DV
                        if cn + DV ~= cn1
                            prob0 = prob0 * (1 - VN(vn, cn1));
                            prob1 = prob1 * VN(vn, cn1);
                        end
                    end
                    for cni = 1: DC
                        if CN(VN(vn,cn), cni) == vn
                            CN(VN(vn, cn), cni + DC) = prob1 / prob1 + prob0;
                        end
                    end
                end
            end 

            % Now it's time to estimate the codeword
            % are you able to find the similarity between VN To CN and this
            for vn = 1: n
                if(VN(vn, 2 * DV + 1) == 0)
                    lamda1 = -Inf;
                elseif (VN(vn, 2 * DV + 1) == 1)
                    lamda1 = Inf;
                elseif (VN(vn, 2 * DV + 1) == 9)
                    lamda1 = 0;
                end
                for cn = DV + 1 : 2 * DV 
                    lamda1 = lamda1 + log(VN(vn, cn)) / (1 - VN(vn, cn));
                end
                if lamda1 == Inf
                    decodedCur(1, vn) = 1;
                elseif lamda1 == 0
                    decodedCur(1, vn) = 9;
                elseif lamda1 == -Inf
                    decodedCur(1, vn) = 0;
                end
            end
            % Check Is This decodedCur having erasure
            flag = 0;
            if decodedCur(1, :) == 9
                flag = 1;
            end
            % if it does not have erasure then procced in this block
            if flag == 0
                % no change in prev and cur then there no need to procced
                % further
                if H.H * transpose(decodedCur) == codeword
                    disp("YES");
                    break;
                end
                %if decodedCur == decodedPrev
                 %   break;
                %end
            end
            decodedPrev = decodedCur;
        end
        % after completing all iterations
        decoded = decodedCur;

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
title('BEC Soft Decision - 202101224');
grid on;


% Algorithmatic Convergence 
figure;
hold on;
for i = 1: length(p)
    plot(algocon(i,:) / n,'DisplayName',['p=' num2str(p(i))]);
end
legend('show');
xlim([1,50]);
xlabel('Iteration Index');
ylabel('Error Probability');
title('Algorithmatic Convergence (BEC Soft Decision) - 202101224');
toc
