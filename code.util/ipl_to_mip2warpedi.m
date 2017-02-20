function mip2warped_omni_i = ipl_to_mip2warpedi(ipl_depth)

% i=919 is 28% 
% i=744 is 62% 

mip2warped_omni_i = (ipl_depth-0.28) * (744-919) / (0.62-0.28) + 919;

end
