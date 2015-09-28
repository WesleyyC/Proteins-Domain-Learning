function [ tf ] = converge( M1,M2,e )
    tf = abs(sum(sum(M1-M2)))<e;
end