function [ reformat ] = deloutbelow( reformat,var,thres )
%DELOUTABOVE delete values below the given threshold, for column 'var'

ii=reformat(:,var)<thres;
reformat(ii,var)=NaN;

end

