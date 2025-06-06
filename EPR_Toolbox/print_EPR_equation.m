function print_EPR_equation (X,best_par,pop)

% display the EPR equation
if size(X,2) == 1
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 2
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 3
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 4
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 5
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 6
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 7
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 8
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 9
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 10
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)]\n', [best_par(end) pop(end,:)])
end


% display the EPR equation
if size(X,2) == 11
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 12
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 13
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 14
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 15
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)*x15^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)*x15^(%3.1f)]\n', [best_par(end) pop(end,:)])
end

% display the EPR equation
if size(X,2) == 16
    fprintf('%6.12f +\n', best_par(1));
    for i = 1:size(best_par,1)-2
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)*x15^(%3.1f)*x16^(%3.1f)] +\n', [best_par(i+1) pop(i,:)]);
    end
    % display the last term
    fprintf('%6.12f * [x1^(%3.1f)*x2^(%3.1f)*x3^(%3.1f)*x4^(%3.1f)*x5^(%3.1f)*x6^(%3.1f)*x7^(%3.1f)*x8^(%3.1f)*x9^(%3.1f)*x10^(%3.1f)*x11^(%3.1f)*x12^(%3.1f)*x13^(%3.1f)*x14^(%3.1f)*x15^(%3.1f)*x16^(%3.1f)]\n', [best_par(end) pop(end,:)])
end