function [timing] = lmeasure(filename,mode)
A=RBread(filename);
switch mode
    case 0
        [p count] = analyze (A,'sym',3);
        tic
        L=lchol(A(p,p));
        timing=toc;
    case 1
        [p count] = analyze (A,'sym',3);
        tic
        L=lchol(A);
        timing=toc;
    case 2
        tic
        L=lchol(A);
        timing=toc;
    case 3
        [p count] = analyze (A,'sym',3);
        AT=A(p,p);
        t = zeros(1,10);
        for i = 1:10
            tic
            L=lchol(AT);
            t(i)=toc;
        end
        timing=mean(t);
    otherwise
        disp("Unknown option.")
end
end

