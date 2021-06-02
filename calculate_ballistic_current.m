function current_ballistic = calculate_ballistic_current(bias,Tmatrix)
    e = 1.602e-19;
    h = 6.626e-34;
    sizes = size(Tmatrix);
    n = sizes(1);
    m = sizes(2);
    numfermis = sizes(3);
    current = 0;
    Energy = linspace(0,bias,numfermis);
    for i = 1:n
        for j = 1:m
            Tmn = Tmatrix(i,j,:);
            Tmn = reshape(Tmn,[1,numfermis]);
            current = current+(2*e/h)*trapz(Energy,Tmn);
        end
    end
    current_ballistic = current;
    