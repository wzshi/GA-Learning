%**********************************
%多目标优化问题
%
%min f1 = (x1^2+x2^2)/4
%min f2 = x1(1-x2)+10
%s.t. 1<=x1<=4,1<=x2<=2
%**********************************

NIND = 100;
MAXGEN = 50;
NVAR = 2;
PRECI = 20;
GGAP = 0.9;
trace1 = zeros(MAXGEN,2);
trace2 = zeros(MAXGEN,2);
trace3 = zeros(MAXGEN,2);
FieldD = [rep([PRECI],[1,NVAR]);[1,1;4,2];rep([1;0;1;1],[1,NVAR])];
Chrom = crtbp(NIND,NVAR*PRECI);
v = bs2rv(Chrom,FieldD);
gen = 1;
while gen<=MAXGEN
    [NIND,N] = size(Chrom);
    M = fix(NIND/2);
    [ObjV1,ans] = ExampleOfChap7_12_f(v(1:M,:));
    FitnV1 = ranking(ObjV1);
    SelCh1 = select('sus',Chrom(1:M,:),FitnV1,GGAP);
    [ans,ObjV2] = ExampleOfChap7_12_f(v(M+1:NIND,:));
    FitnV2 = ranking(ObjV2);
    SelCh2 = select('sus',Chrom(M+1:NIND,:),FitnV2,GGAP);
    SelCh = [SelCh1;SelCh2];
    SelCh = recombin('xovsp',SelCh,0.7);
    Chrom = mut(SelCh);
    v = bs2rv(Chrom,FieldD);
    [f1,f2] = ExampleOfChap7_12_f(v);
    trace1(gen,1) = min(f1);
    trace1(gen,2) = sum(f1)/length(f1);
    trace2(gen,1) = min(f2);
    trace2(gen,2) = sum(f2)/length(f2);
    trace3(gen,1) = min(f1+f2);
    trace3(gen,2) = sum(f1)/length(f1)+sum(f2)/length(f2);
    gen = gen+1;
end
figure(1);clf;
plot(trace1(:,1));hold on;
plot(trace1(:,2),'-.');
plot(trace1(:,1),'.');
plot(trace1(:,2),'.');grid on;
legend('解的变化','种群均值的变化');
xlabel('迭代次数');ylabel('目标函数值');

figure(2);clf;
plot(trace2(:,1));hold on;
plot(trace2(:,2),'-.');
plot(trace2(:,1),'.');
plot(trace2(:,2),'.');grid on;
legend('解的变化','种群均值的变化');
xlabel('迭代次数');ylabel('目标函数值');

figure(3);clf;
plot(trace3(:,1));hold on;
plot(trace3(:,2),'-.');
plot(trace3(:,1),'.');
plot(trace3(:,2),'.');grid on;
legend('解的变化','种群均值的变化');
xlabel('迭代次数');ylabel('目标函数值');

figure(4);clf;
plot(f1);
hold on;
plot(f2);
grid on;



    
    