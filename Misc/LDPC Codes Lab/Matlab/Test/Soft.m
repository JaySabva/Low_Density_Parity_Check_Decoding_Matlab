%SABVA JAY DILIPBHAI
%202101224
%Soft Decoding

%From H Matrix Count N, u, k
H = load("H1.mat");
u = size(H.H, 1);
n = size(H.H, 2);
k = n - u;

%Find the Degree Of Check node, and Degree Of Variable Node
DC = sum(H.H(1,:));
DV = sum(H.H(:,1));

%Valid CodeWord
codeword = zeros(1, n);

%VN and CN Tanner Graph
VN = zeros(n, 2* DV + 1);
CN = zeros(u, 2*DC);

%We need to fill vn and cn so that we can decide which vn is connected to
%which vn and vice versa

%VN Array
for vn = 1:n
    idx = 1;
    for cn = 1:u
        if H.H(cn, vn) == 1
            VN(vn, idx) = cn;
            idx = idx + 1;
        end
    end
end

%CN Array
for cn = 1:u
    idx = 1;
    for vn = 1:n
        if H.H(cn, vn) == 1
            CN(cn, idx) = vn;
            idx = idx + 1;
        end
    end
end

p = 0: 0.02:1;
tmax = 25;
Nsim = 1000;
successar = zeros(size(p));
inx = 1;
algocon = zeros(length(p), tmax);
for Pe = p
    success = 0;
    erasuresum = zeros(1, tmax);
    for L = 1:Nsim

        %BEC Channel
        revcode = zeros(1, n);
        for i = 1:n
            r = rand;
            if (r <= Pe)
                revcode(i) = 9;
            else
                revcode(i) = codeword(i);
            end
        end
        %I am going to decode it at any cost
        [decoded,erasureit] = TannerGraph(VN, CN, revcode, DC, DV, n, u, tmax);
        
        %Say My Name I am TannerGraph Decoding I can decode anything
        if decoded == codeword
            success = success + 1;
        else
            %Next Time I will try best
        end
        erasuresum = erasuresum + erasureit;
    end  
    erasureavg = erasuresum / Nsim;
    algocon(inx,:) = erasureavg;
    successar(inx) = (success / Nsim);

    inx = inx + 1;
end

%Probability Of Successfull Decoding
plot(p,successar);
xlabel('Error Proability Of BEC Channel');
ylabel('Success Probability Of Decoding');
grid on;

%{
%Algorithmatic Convergence
figure;
hold on;
ylim([-1, 10]);
xticks(0:3:50);
for i = 1:length(p)
    plot(1:tmax, algocon(i,:));
end
xlabel('Iteration Index');
ylabel('Average Number of Erasures');
title('Convergence of Tanner Graph Decoding');
%}

%Dont Touch If it work Please-----------<*|*>
function [decoded ,erasureit] = TannerGraph(VN, CN, revcode, DC, DV, n, u,tmax)
    %e iteration ma ketli error che
    erasureit = zeros(1, tmax);
    %load each variable node
    for vn = 1: n
        VN(vn, 2 * DV + 1) = revcode(vn);
    end

    %Loading Each CN With VN For 0th iteration
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
    decodedPrev = VN(:,2 * DV + 1);
    decodedCur = ones(n) * 9;
    %itreration
    for t = 1: tmax

        %count how many erasure
        for er = 1: n
            if  decodedPrev(n)== 9
                erasureit(t) = erasureit(t) + 1;
            end
        end

        %Belief Propagation 
        % CN TO VN
        for cn = 1:u
            for vn = 1:DC
                prob = 0.5;
                for vnvalue = DC + 1: 2 * DC
                    if vn + DC ~= vnvalue
                        prob = prob * (1 - 2 * CN(cn, vnvalue));
                    end
                end
                prob = 0.5 + prob;
                for vn1 = 1: DV
                    if VN(CN(cn,vn), vn1) == cn
                        VN(CN(cn, vn), vn1 + DV) = 1 - prob;
                    end
                end
            end
        end

        %VN TO CN
        for vn = 1:n
            for cn = 1:DV
                if(VN(vn, 2 * DV + 1) == 0)
                    lamda1 = 0;
                elseif (VN(vn, 2 * DV + 1) == 1)
                    lamda1 = 1;
                elseif (VN(vn, 2 * DV + 1) == 9)
                    lamda1 = 0.5;     
                end
                for cn1 = DV + 1: 2 * DV
                    if cn + DV ~= cn1
                        lamda1 = lamda1 * (VN(vn, cn1)) / (1 - VN(vn, cn1));
                    end
                end
                for cni = 1: DC
                    if CN(VN(vn,cn), cni) == vn
                        CN(VN(vn, cn), cni + DC) = lamda1;
                    end
                end
            end
        end 

        %Decoded
        for vn = 1: n
            if(VN(vn, 2 * DV + 1) == 0)
                lamda = 0;
            elseif (VN(vn, 2 * DV + 1) == 1)
                lamda = 1;
            elseif (VN(vn, 2 * DV + 1) == 9)
                lamda = 0.5;
            end
            for cn = DV + 1 : 2 * DV 
                lamda = lamda * (VN(vn, cn)) / (1 - VN(vn, cn));
            end
            if lamda > 0.5
                cHat(vn) = 1;
            elseif lamda == 0.5
                cHat(vn) = 9;
            elseif lamda < 0.5
                cHat(vn) = 0;
            end
        end
        erasu = find(cHat == 9);
        if(size(erasu) == 0)
        if(cHat == cHat_prev)
            break;
        end
        end
        cHat_prev = cHat;
    end
    %after completing all iterations
    decoded = cHat;
end

