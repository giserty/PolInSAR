function [fai0,p] = line_fai(PDH,PDL,Opt1,Opt2,Opt3,HH,HV,VV,HHpVV,HHmVV,n_pol,m,n)
if n_pol == 2
    [fai0,p] = PD_fai0(PDH,PDL,m,n);
elseif n_pol == 5
    [fai0,p] = PDSVD_fai0(PDH,PDL,Opt1,Opt2,Opt3,m,n);
else
    [fai0,p] = PDPauli_fai0(PDH,PDL,HH,HV,VV,HHpVV,HHmVV,m,n);
end
end