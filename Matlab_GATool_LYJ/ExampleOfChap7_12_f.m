function [ f1 f2 ] = ExampleOfChap7_12_f( x )
%EXAMPLEOFCHAP7_12_F1 此处显示有关此函数的摘要
%   此处显示详细说明
f1 = x(:,1).*x(:,1)/4+x(:,2).*x(:,2)/4;
f2 = x(:,1).*(1-x(:,2))+10;
end

