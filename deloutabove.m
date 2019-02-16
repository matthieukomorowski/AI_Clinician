function [ reformat ] = deloutabove( reformat, var, thres )
%DELOUTABOVE delete values above the given threshold, for column 'var'

ii=reformat(:,var)>thres;
reformat(ii,var)=NaN;


end

