function [p,v,t] = calc_ref(tra,ss,se,vs,as,dt)
    if((se-ss)/(vs^2/as) > 1)
        t1 = 0:dt:vs/as-dt;
        t2 = t1(end):dt:(se-ss)/vs-dt;
        t3 = t2(end):dt:vs/as+t2(end);
        vs1 = as*t1;
        vs2 = vs*ones(size(t2));
        vs3 = vs - as*(t3-t2(end));
        s1 = 0.5*as*t1.^2 +ss;
        s2 = vs*(t2 -t1(end)) + s1(end);
        s3 = vs*(t3-t2(end)) - 0.5*as*(t3-t2(end)).^2 + s2(end);
        ta = [t1 t2 t3];
        vsa =[vs1 vs2 vs3];
        sa = [s1 s2 s3];
    else
        T = 2*sqrt((se-ss)/as);
        t1 = 0:dt:T/2-dt;
        t2 = T/2:dt:T;
        vs1 = as*t1;
        vs2 = as*T/2 - as*(t2-T/2);
        s1 = 0.5*as*t1.^2 + ss;
        s2 = se - 0.5*as*(t2-T).^2;
        ta = [t1 t2];
        vsa = [vs1 vs2];
        sa = [s1 s2];
    end
    t =ta;
    p = tra.p(sa);
    v = tra.v(sa,vsa);
end
%0 < t < t1
%ds = as*t
%s = 0.5*as*t^2 + ss
%ds = vs = as*t1
%-> t1 = vs/as
%-> s1 = 0.5*vs^2/as + ss
%
%t1 < t <t2
%ds = vs
%s = vs*(t-t1) + s1
%s2 = vs*t2 - 0.5*vs^2/as + ss
%
%t2 < t <t3
%ds = vs - as*(t-t2)
%s = vs*(t-t2) - 0.5*as*(t-t2)^2 + s2
%ds = vs - as*(t-t2) = 0
%-> t3 = vs/as + t2
%-> s4 = vs*t2 + ss =se
%-> t2 = (se-ss)/vs

%0 < t <T/2
%ds = as*t
%s = 0.5*as*t^2 + ss
%sm = as*T^2/8 + ss
%sm = as*T/2
%
%T/2 < t < T
%ds = as*T/2 - as*(t-T/2)
%s = se - 0.5*as*(t-T)^2
%dsm = se - as*T^2/8
%as*T^2/8 + ss = se - as*T^2/8
%as*T^2/4 = se-ss
%T = 2*sqrt((se-ss)/as);

