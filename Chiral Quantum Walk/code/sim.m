function [qqq, ttt] = sim()
% Simulates a quantum walk in an exciton number preserving circuit.

% Ville Bergholm 2012


global qit

X = qit.sx;
Y = qit.sy;

dim = [2 2 2];


% symmetric and antisymmetric exciton number preserving two-qubit Hamiltonians
HS = 0.5*(kron(X, X) +kron(Y, Y)); % symmetric, real
HA = 0.5*(kron(X, Y) -kron(Y, X)); % antisymmetric, imaginary

% mixed symmetry case
theta = 0.3145 * pi;
H = cos(theta) * HS +sin(theta) * HA;

% single-exciton subspace
ind = [5 3 2];

% 2q Hamiltonian
gate_dim = {[2 2], [2 2]};
H = lmap(H, gate_dim);

% full Hamiltonians
h = gate.two(H, [1 2], dim).data;
H12 = full(h(ind, ind));

h = gate.two(H, [2 3], dim).data;
H23 = full(h(ind, ind));

h = gate.two(H, [1 3], dim).data;
H13 = full(h(ind, ind));

H = {H12, H23, H13};


initial = 0;
s0 = state(initial, 3);

desc = sprintf('Initial state: %d,  theta = %g pi', initial, theta/pi)


if true
    % optimize the angles
    % define the optimization problem
    problem.objective = @(x) goal_func(x);
    problem.x0 = abs(0 + randn(1, 3))

    % try to minimise objective function to zero
    [x, cost, exitflag, output] = fminunc(problem.objective, problem.x0);
    t_opt = x.^2;
else
    % all gates use the same time/angle
    t_opt = pi * [1 1 1];
end

N = 200
ttt = linspace(0, 2, N);
res = zeros(3, N);

for k=1:N
    t = ttt(k) * t_opt;
    U = circuit(H, t);
    s = u_propagate(s0, U);
    res(:, k) = s.data;
end
t_opt/pi

figure();
total_t = ttt * sum(t_opt);
plot(total_t/pi, abs(res).^2);
legend('1', '2', '3')
xlabel('total gate time / \pi')
ylabel('probability ampl.')
title(desc)


function err = goal_func(x)
  
  x = x.^2; % no negative times
  U = circuit(H, x);
  s = u_propagate(s0, U);
  temp = s.data;
  err = 1-abs(temp(2))^2;
end
end






function U = circuit(H, t)

  U1 = expm(-1i * t(1) * H{1});
  U2 = expm(-1i * t(2) * H{2});
  U3 = expm(-1i * t(3) * H{3});

  % palindromic gate sequence
  U = U1*U2*U3 * U3*U2*U1;
end
