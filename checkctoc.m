clear all
close all

C = rand(6,6);
C = (C+transpose(C))/2;

c4 = C2toc4(C);
c2 = c4toC2(c4);

c2-C

