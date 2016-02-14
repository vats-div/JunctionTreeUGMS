function CM = cov2corr(CM)

Sdm = diag(diag(CM).^(-0.5));
CM = Sdm * CM * Sdm;

end