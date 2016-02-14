function CM = mycorr(X)

CM = cov(X);
Sdm = diag(diag(CM).^(-0.5));
CM = Sdm * CM * Sdm;

end