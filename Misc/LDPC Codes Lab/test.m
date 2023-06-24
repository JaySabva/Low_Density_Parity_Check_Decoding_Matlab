clear all
H = load("H1.mat");

n = size(H.H, 2);
u = size(H.H, 1);
k = n - u;

% Degree Of Check Node and Variable Node
DC = sum(H.H(1, :));
DV = sum(H.H(:, 1));

codeword = zeros(1, n);

CN = zeros(u, 2 * DC);
VN = zeros(n, 2 * DV + 1);

for vn = 1: n
    idx = 1;
    for cn = 1: u
        if H.H(cn, vn) == 1 % if thery are connected
            VN(vn, idx) = cn;
            idx = idx + 1;
        end
    end
end

for cn = 1: u
    idx = 1;
    for vn = 1: n
        if H.H(cn, vn) == 1
            CN(cn, idx) = vn;
            idx = idx + 1;
        end
    end
end

p = 0:0.02:1;
tmax = 25;
Nsim = 100;
Psuccess = zeros(size(p));
inx = 1;
algocon = zeros(length(p), tmax);

for Pe = p
    success = 0;
    erasure_sum = zeros(1, tmax);
    for L = 1: Nsim
        revcode = zeros(1, n);
        for i = 1:n
            r = rand;
            if r <= Pe
                revcode(i) = 9;
            else
                revcode(i) = codeword(i);
            end
        end
        decoded = ones(1, n) * 9;
        erasure_perit = zeros(1, tmax);

        for vn = 1: n
            VN(vn, 2 * DV + 1) = revcode(vn);
        end

        for cn = 1: u
            for vn = 1: DC
                CN(cn, vn + DC) =  VN(CN(cn, vn), 2 * DV + 1);
            end
        end

        for t = 1: tmax
            for vn = 1:n
                if VN(vn, 2 * DV + 1) == 9
                    erasure_perit(t) = erasure_perit(t) + 1;
                end
            end

            noOfCNCorrect = 0;
            for cn = 1: u

                erasure = find(CN(cn, DC + 1: 2* DC) == 9);

                if size(erasure) == 0
                    noOfCNCorrect = noOfCNCorrect + 1;
                elseif size(erasure) == 1

                    vn = CN(cn, erasure(1));
                    sum = 0;
                    for vnval = DC + 1: 2 * DC
                        if(CN(cn, vnval) == 1 || CN(cn, vnval) == 0)
                            sum = sum + CN(cn, vnval);
                            sum = mod(sum, 2);
                        end
                    end
                    for vni = 1:DV
                        if VN(vn, vni) == cn
                            if sum == 0
                                VN(vn, DV + 1) = 0;
                            elseif sum == 1
                                VN(vn, DV + 1) = 1;
                            end
                        end
                    end
                    if noOfCNCorrect == u
                        break;
                    end
                else
                    for vnx = erasure
                        for vni = 1:DV
                            if VN(CN(cn, vnx), vni) == cn
                                VN(CN(cn, vnx), vni + DV) = 9;
                            end
                        end
                    end
                end
            end
            for vn = 1:n
                for cn = 1:DV
                    cnt0 = 0;
                    cnt1 = 0;
                    cnt9 = 0;
                    for cn1 = DV + 1: 2 * DV
                        if cn + DV ~= cn1
                            if(VN(vn, cn1)) == 0
                                cnt0 = cnt0 + 1;
                            elseif VN(vn, cn1) == 1
                                cnt1 = cnt1 + 1;
                            else 
                                cnt9 = cnt9 + 1;
                            end
                        end
                    end    
                   for cni = 1:DC
                       if CN(VN(vn, cn), cni) == vn
                           if cnt9 == DV - 1
                                CN(VN(vn, cn), cni + DC) = 9;
                           elseif cnt1 > cnt0
                               CN(VN(vn, cn), cni + DC) = 1;
                           else
                               CN(VN(vn, cn), cni + DC) = 0;
                           end
                       end
                   end
                end
            end

             for vn = 1:n
                    cnt0 = 0;
                    cnt1 = 0;
                    cnt9 = 0;
                for cn = DV + 1: 2 * DV 

                    if VN(vn, cn) == 9
                        cnt9 = cnt9 + 1;
                    elseif VN(vn, cn) == 1
                        cnt1 = cnt1 + 1;
                    else 
                        cnt0 = cnt0 + 1;
                    end
                end
                if cnt9 == DV
                    VN(vn, 2 * DV + 1) = 9;
                elseif cnt1 > cnt0
                    VN(vn, 2 * DV + 1) = 1;
                else 
                    VN(vn, 2 * DV + 1) = 0;
                end
             end
             decoded(1, :) = VN(:,2 * DV + 1);
             if H.H * transpose(decoded) == codeword
                 break;
             end
        end

        % Say My Name - LDPC + TanneGraph (LDPC On Erasure Gone)
        if H.H * transpose(decoded) == codeword % decoded == codeword
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
    plot(algocon(i,:) / n,'DisplayName',['p=' num2str(p(i))]);
end
legend('show');
xlim([1,25]);
xlabel('Iteration Index');
ylabel('Error Probability');    
title('Algorithmatic Convergence (BEC Hard Decision) - 202101224');
grid on;