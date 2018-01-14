function out = gpeFun(psi_fun,potential,model_pars)
dens    = abs(psi_fun).^2;
out     = 1i*((-U_now*dens-model_pars.K3_im*dens.^2-potential).*psi_fun);
end