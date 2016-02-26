function [ output_args ] = p_sigmf( x,a,c )
%P_SIGMF Summary of this function goes here
%   Detailed explanation goes here

    output_args = 1/(1+exp(-a*(x-c)));
end

