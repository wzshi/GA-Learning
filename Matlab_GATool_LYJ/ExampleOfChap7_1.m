figure(1);
% sym variable;
fplot('variable.*sin(10*pi*variable)+2.0', [-1,2]);
NIND = 40;        %Number of individuals
MAXGEN = 50;      %Maximum number of generations
PRECI = 20;       %Precision of variables
GGAP = 0.9;       %Generation gap

trace = zeros(2,MAXGEN);
FieldD = [20;-1;2;1;0;1;1];         %Build field descriptor
Chrom = crtbp(NIND, PRECI);         %
gen = 0;
variable = bs2rv(Chrom, FieldD);
ObjV = variable.*sin(10*pi*variable)+2.0;

while gen<MAXGEN,
    FitnV = ranking(-ObjV);         %Assign fitness values
    SelCh = select('sus', Chrom, FitnV, GGAP);
    SelCh = recombin('xovsp', SelCh, 0.7);
    SelCh = mut(SelCh);
    variable = bs2rv(SelCh, FieldD);
    ObjVSel = variable.*sin(10*pi*variable)+2.0;
    [Chrom ObjV] = reins(Chrom, SelCh, 1, 1, -ObjV, -ObjVSel);
    ObjV = -ObjV;
    gen = gen+1;
    
    [Y,I] = max(ObjV), hold on;
    variable = bs2rv(Chrom, FieldD);
    plot(variable(I), Y, 'ro');
    trace(1,gen) = max(ObjV);
    trace(2,gen) = sum(ObjV)/length(ObjV);
end

variable = bs2rv(Chrom, FieldD);
hold on, grid;
plot(variable', ObjV', 'b*');
figure(2);
plot(trace(1,:)');
hold on;
plot(trace(2,:)', '-.');
grid;
legend('解的变化','种群均值的变化');
    