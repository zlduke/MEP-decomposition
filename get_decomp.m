function [X,S_est,Vpp_est,Ipk_est,Ipk_raw] = get_decomp(S, v)
%   S: raw signal to be detected
%   v: template signal

[Nemg,Lemg] = size(S);
N_offset = Lemg - length(v) + 1;
options = optimoptions('quadprog','Display','off');

v = v/range(v);     %   normalize
% v = v-mean(v);    %   better to not zero-mean it 

Tv = v2T(v, N_offset);
X = zeros(N_offset,Nemg);

% ppm = ParforProgressbar(Nemg);

for i = 1:Nemg
    s = S(i,:);
%     s = s-mean(s);
    H = Tv'*Tv;   % scaled by 1/2
    f = -Tv'*s';      % scaled by 1/2
    A = -eye(N_offset);
    b = zeros(N_offset,1);
    x = quadprog(H,f,A,b,[],[],[],[],[],options);
    X(:,i) = x;
%     ppm.increment();
end

S_est = (Tv*X)';
Vpp_est = range(S_est,2);
Ipk_est = Vpp_est*0;
Ipk_raw = Vpp_est*0;
for i = 1: Nemg
    Ipk_est(i) = median(find(S_est(i,:) == max(S_est(i,:))));
    Ipk_raw(i) = median(find(S(i,:) == max(S(i,:))));
end