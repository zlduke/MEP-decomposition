function [v,record_MSE] = get_MUAP(S, v_init, Nitr)
%   modified at 20200506 from convStudy_adaptiveMUAP_v2.m


[Nemg,Lemg] = size(S);
Lv = length(v_init);
N_offset = Lemg - Lv + 1;
record_MSE = zeros(Nitr,1);

%   Part 1: get lin-decomp for each of the EMG measurement
%       optimization: min|T(v)x - s|, x>0 for each s in S (i.e., EMG_scope)
%       output      : X (Nidx_offset * Nemg)

options = optimoptions('quadprog','Display','off');

v = v_init;
v = v/range(v);
v = v-mean(v);

for itr = 1:Nitr %roughly takes 20 iterations

    Tv = v2T(v, N_offset);
    X = zeros(N_offset,Nemg);
    for i = 1:Nemg
        s = S(i,:);
%         s = s-mean(s);
        H = Tv'*Tv;         % scaled by 1/2
        f = -Tv'*s';        % scaled by 1/2
        A = -eye(N_offset);
        b = zeros(N_offset,1);
        x = quadprog(H,f,A,b,[],[],[],[],[],options);
        X(:,i) = x;
    end
    
    record_MSE(itr) = mean(mean((Tv*X - S').^2));

    %   Part 2: estimate actual MUAP from lin-decomp results
    %       optimization: min|T(x)v - s| for each x in X (i.e., MUAP weights)
    %       output      : v_opt for each MUAP (test first)

    TT = 0;
    Ts = 0;
    for i = 1:Nemg
        s = S(i,:);     % 1 x LENemg
        x = X(:,i);     % Nidx_offset x Nemg
        Tx = v2T(x, Lv); % (Nidx_offset+LENv-1) x (LENv)
        TT = TT + Tx'*Tx;
        Ts = Ts + Tx'*s';
    end
    new_v = sqrtm(TT)\Ts;
    new_v = new_v/range(new_v);
    new_v = new_v-mean(new_v);

    v = new_v;
end
