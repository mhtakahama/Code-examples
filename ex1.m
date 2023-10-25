clear all; echo off; close all force; clc; %clear another variables
k=25
q=1000
tp1=423
tp2=323
x=.2

A=[3 -1 0 0 0; -1 2 -1 0 0;0 -1 2 -1 0;0 0 -1 2 -1; 0 0 0 -1 3],
B=[2*tp1, 0, 0 ,0 ,2*tp2]
C=inv(A)*B'-273
