function Qit=stretch_out(Ditd,Dita,E,Efermi,Vt)
Q_SI=1.6e-19;
Qit=ones(1,length(Efermi));
for i=1:length(Efermi)
    fermi=1./(1+exp((E-Efermi(i))/Vt));
    funct=Ditd.*(1-fermi)-Dita.*(fermi);
    Qit(i)=Q_SI*trapz(E,funct);
end
end