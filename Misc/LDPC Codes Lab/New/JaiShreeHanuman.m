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

% Probability Array Because of algo-conv
p = 0: 0.02: 1;
p = [0.3, 0.4 ,0.5, 0.51, 0.52];
tmax = 45; % max iteration allowed
Nsim = 1000; % No of simulation
Psuccess = zeros(size(p)); % store success for each p and all of Nsim
inx = 1; %index
algocon = zeros(length(p), tmax); % algo convergence 2d array for each probabilty there erasure per iteration

for Pe = p
    disp(Pe);
    success = 0;
    erasure_sum = zeros(1, tmax);

    for L = 1: Nsim

        % BEC Channel
        revcode = zeros(1, n); % Inital Recieved Codeword
        for i = 1: n
            r = rand;
            if r <= Pe
                revcode(1, i) = 9; % erasure
            else
                revcode(1, i) = codeword(1, i); % as it is
            end
        end

        for vn = 1: n
            VN(vn, 2 * DV + 1) = revcode(1, vn);
        end

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
        decodedPrev = revcode(1,:);
        erasure_perit = zeros(1, tmax);
        decodedCur = ones(1,n) * 9;
        for t = 1: tmax
            for vn= 1:n
                if decodedPrev(1, vn) == 9
                    erasure_perit(1, t) = erasure_perit(1, t) + 1;
                end
            end

            for cn = 1:u
                for vn = 1:DC
                    prob0 = 1;
                    for vnvalue = DC + 1: 2 * DC
                        if vn + DC ~= vnvalue
                            prob0 = prob0 * (1 - 2 * CN(cn, vnvalue));
                        end
                    end
                    prob0 = (1 + prob0) / 2;
                    for vn1 = 1: DV
                        if VN(CN(cn,vn), vn1) == cn
                            VN(CN(cn, vn), vn1 + DV) = 1 - prob0;
                            break;
                        end
                    end
                end
            end

            for vn = 1: n
                if VN(vn, 2 * DV + 1) == 0
                    prob1 = 0;
                elseif VN(vn, 2 * DV + 1) == 1
                    prob1 = 1;
                else
                    prob1 = 0.5;
                end
                prob0 = 1 - prob1;
                for cn = 1: DV
                    p1 = 1;
                    p0 = 1;
                    for cn1 = DV + 1:2 * DV
                        if cn + DV ~= cn1
                            p1 = p1 * VN(vn, cn1);
                            p0 = p0 * (1 - VN(vn, cn1));
                        end
                    end
                    p1 = prob1 * p1;
                    p0 = prob0 * p0;
                    for cni = 1: DC
                        if CN(VN(vn,cn), cni) == vn
                            CN(VN(vn, cn), cni + DC) = p1 / (p1 + p0);
                        end
                    end
                end
            end

            for vn = 1: n
                if VN(vn, 2 * DV + 1) == 0
                    prob1 = 0;
                elseif VN(vn, 2 * DV + 1) == 1
                    prob1 = 1;
                else
                    prob1 = 0.5;
                end
                prob0 = 1 - prob1;
                for cn1 = DV + 1:2 * DV
                    prob1 = prob1 * VN(vn, cn1);
                    prob0 = prob0 * (1 - VN(vn, cn1));
                end
                if prob1 / prob0 > 1
                    decodedCur(1, vn) = 1;
                elseif prob1 == prob0
                    decodedCur(1, vn) = 9;
                elseif prob1/ prob0 < 1
                    decodedCur(1, vn) = 0;
                end
            end

            erasure = find(decodedCur(1, :) == 9);
            if size(erasure, 2) == 0
                break;
            end
            decodedPrev = decodedCur;
        end
        decoded = decodedCur;
        if decoded == codeword
            success = success + 1;
        else
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

