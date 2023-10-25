clear all; echo off; close all force; clc; %clear another variables
k=25
q=1000
tp1=150
tp2=50     
r1=0.5
r2=1
dr=.1/2
r=r1:dr*2:r2

n=length(r)

A=[2*(r(1)-dr)+r(1)+dr -(r(1)+dr) 0 0 0; -(r(2)-dr) (r(2)+dr+r(2)-dr) -(r(2)+dr) 0 0;0 -(r(3)-dr) (r(3)+dr+r(3)-dr) -(r(3)+dr) 0;0 0 -(r(4)-dr) (r(4)+dr+r(4)-dr) -(r(4)+dr); 0 0 0 -(r(5)-dr) 2*(r(5)+dr)+r(5)-dr],
B=[2*(r(1)-dr)*tp1+q/k*dr^2, q/k*dr^2, q/k*dr^2 , q/k*dr^2 ,2*(r(5)+dr)*tp2+q/k*dr^2]
C=inv(A)*B'
plot(C)
