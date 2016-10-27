%定义遗传算法参数
NIND = 40;
MAXGEN = 2;
NVAR = 20;
PRECI = 20;
GGAP = 0.9;
trace = zeros(MAXGEN,2);
%建立区域描述器（Build field descriptor）
FieldD = [rep([PRECI],[1,NVAR]);rep([-512;512],[1,NVAR]);rep([1;0;1;1],[1,NVAR])];
Chrom = crtbp(NIND,NVAR*PRECI);
gen = 0;
ObjV = objfun1(bs2rv(Chrom,FieldD));
figure(1);
while gen<MAXGEN
    FitnV = ranking(ObjV);
    SelCh = select('sus',Chrom,FitnV,GGAP);
    SelCh = recombin('xovsp',SelCh,0.7);
    SelCh = mut(SelCh);
    ObjVSel = objfun1(bs2rv(SelCh,FieldD));
    [Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
    gen = gen+1;
    plot(ObjV,'*');hold on;grid on;
    trace(gen,1) = min(ObjV);
    trace(gen,2) = sum(ObjV)/length(ObjV);
end
figure(2);
plot(trace(:,1));
hold on;
plot(trace(:,2),'r');
grid on;
legend('解的变化','种群均值的变化');