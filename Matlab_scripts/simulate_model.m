function batch=simulate_all(t,Y,Vh, Kh,a, Qlag, NQlag, epsilon_NQ, epsilon_Q, epsilon,epsilonAG, reusablility_ratio, dQ, dNQ, b)
% here we consider Q and NQ are cells that are programmed to become Q and
% NQ but at the beginning they are the same
G =Y(1);           %glucose (S.cerevisiae limiting nutrient)
Q =Y(2);           %quiescent cells
NQ_Q =Y(3);        %non-quiescent cells: descendants of Q
NQ_NQ =Y(4);       %non-quiescent cells: descendants of NQ
Qprep = Y(5);      %non-quiescent cells preparing to become Q <do not proliferate>
D = Y(6);          %dead cells
A = Y(7);          %aminoacid reservuar
%%%%%%%%%%%%%%%%%
% handle the lags
if t >= Qlag
    JgQ  = Vh*G./(Kh+G); 
else 
    JgQ = 0;
end

if t >= NQlag
    JgNQ = Vh*G./(Kh + G);
else 
    JgNQ = 0;
end

%%%%%%%%%%%%%%%%%%%%%%
% rate of switching from NQ_NQ and NQ_Q to Q
epsilon_of_glucose = switch_rate_of_glucose(G);
% we introduce Qprep to slow down the dynamics


Gdot     = -JgQ.*NQ_Q - JgNQ.*NQ_NQ + epsilonAG*A; 
Qdot     = epsilon*Qprep -dQ*Q;
NQ_Qdot  = a*JgQ.*NQ_Q - epsilon_of_glucose*epsilon_Q*NQ_Q - dNQ*NQ_Q;
NQ_NQdot = a*JgNQ.*NQ_NQ - epsilon_of_glucose*epsilon_NQ*NQ_NQ - dNQ*NQ_NQ;
Qprepdot = b*(epsilon_of_glucose*epsilon_Q*NQ_Q + epsilon_of_glucose*epsilon_NQ*NQ_NQ) - epsilon*Qprep -dNQ*Qprep;
Ddot     = dNQ*(NQ_NQ + NQ_Q + Qprep) + dQ*Q; % no need to track this
Adot     = reusablility_ratio*(dNQ*(NQ_NQ + NQ_Q + Qprep) + dQ*Q) - epsilonAG*A ; 
batch  =[Gdot;Qdot; NQ_Qdot;NQ_NQdot;Qprepdot;Ddot;Adot];
end 

