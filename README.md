# DFS

```
A = [1 1 1 0;
    2 -3 0 1];
b = [1;3];
C = [2 3 0 0];
[n_constraints , n_vars] = size(A);
Zmax = -inf;
Xopt = [];

comb = nchoosek(1:n_vars, n_constraints);

for i = 1:size(comb,1)
    B = A(:, comb(i,:));
    if det(B) ~= 0
        xB = B \b;
        X = zeros(n_vars, 1);
        X(comb(i,:)) = xB;

        if all(X>= 0)
            Z = C*X;
            fprintf('BFS: %s  Z = %.2f\n', mat2str(X'), Z);
            if Z > Zmax
                Zmax = Z;
                Xopt = X;
            end
        end
    end
end
fprintf('\n Optimal Solution:\n')
disp(Xopt');
fprintf('Maximum Z = %.2f\n', Zmax)
```

# Graph

```
clc
clear
format short

% -x_1 + 3x_2 <= 10
%  x_1 +  x_2 <= 6
%  x_1 -  x_2 <= 2

A = [-1  3;
      1  1;
      1 -1];
B = [10;
     6;
     2];

x1 = 0:0.01:max(B);
x2_1 = (B(1) - A(1,1)*x1) / A(1,2);
x2_2 = (B(2) - A(2,1)*x1) / A(2,2);
x2_3 = (B(3) - A(3,1)*x1) / A(3,2);

x2_1 = max(0, x2_1);
x2_2 = max(0, x2_2);
x2_3 = max(0, x2_3);

figure
plot(x1, x2_1, 'r', x1, x2_2, 'b', x1, x2_3, 'g', 'LineWidth', 1.5)
hold on
grid on
xlabel('x_1')
ylabel('x_2')
title('Feasible Region and Optimal Solutions')

% --- Find all candidate corner points ---
pt = [0 0];   % include origin

% Pairwise intersections of constraint lines
for i = 1:size(A,1)
    for j = i+1:size(A,1)
        M = [A(i,:); A(j,:)];
        if abs(det(M)) > 1e-10   % check non-singular
            X = M \ [B(i); B(j)];
            pt = [pt; X'];
        end
    end
end

% Intersections of each constraint with x1=0 and x2=0 axes
for i = 1:size(A,1)
    if abs(A(i,2)) > 1e-10
        pt = [pt; 0, B(i)/A(i,2)];
    end
    if abs(A(i,1)) > 1e-10
        pt = [pt; B(i)/A(i,1), 0];
    end
end

% --- Filter feasible points ---
points = unique(pt, 'rows');
feasible = [];
for i = 1:size(points,1)
    x = points(i,:)';
    if all(A*x <= B + 1e-6) && all(x >= -1e-6)
        feasible = [feasible; max(0, x')];   % snap near-zero to 0
    end
end

P = unique(feasible, 'rows');

% --- Plot feasible region ---
if size(P,1) >= 3
    K = convhull(P(:,1), P(:,2));
    fill(P(K,1), P(K,2), 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'k')
end
plot(P(:,1), P(:,2), 'ko', 'MarkerFaceColor', 'k')

% --- Objective function ---
C = [1 5];     % Z = x1 + 5*x2
Z = P * C';

% Display corner point table
Corner_Values = array2table([P Z], 'VariableNames', {'x1','x2','Z'});
disp('Corner points and objective values:')
disp(Corner_Values)

[maxZ, idxMax] = max(Z);
[minZ, idxMin] = min(Z);

fprintf('Maximum Z = %.4f at (%.4f, %.4f)\n', maxZ, P(idxMax,1), P(idxMax,2));
fprintf('Minimum Z = %.4f at (%.4f, %.4f)\n', minZ, P(idxMin,1), P(idxMin,2));

plot(P(idxMax,1), P(idxMax,2), 'r*', 'MarkerSize', 12, 'LineWidth', 2)
plot(P(idxMin,1), P(idxMin,2), 'b*', 'MarkerSize', 12, 'LineWidth', 2)

legend('-x_1 + 3x_2 = 10', ...
       'x_1 + x_2 = 6', ...
       'x_1 - x_2 = 2', ...
       'Feasible Region', ...
       'Corner Points', ...
       'Maximum Z', ...
       'Minimum Z', ...
       'Location', 'best')
xlim([0 max(B)+1])
ylim([0 max(B)+1])
hold off
```
