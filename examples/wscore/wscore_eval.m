function [result] = evalinmodel (p,ts,xs,expr)
%  wscore

   s0 = p(1);
   vt = p(2);
   we = p(3);
   wr = p(4);
   wp = p(5);
   wq = p(6);
   thetar = p(7);
   thetax = p(8);
   nr = p(9);
   nx = p(10);
   gmax = p(11);
   Kp = p(12);
   vm = p(13);
   Km = p(14);
   ns = p(15);
   kq = p(16);
   Kt = p(17);
   nq = p(18);
   cl = p(19);
   k_cm = p(20);
   b = p(21);
   dm = p(22);
   kb = p(23);
   ku = p(24);
   M = p(25);

[npoints,~] = size(xs);
result = zeros(npoints,1);

for asdf = 1:npoints
  x = xs(asdf,:);
  t = ts(asdf);
     si = x(1);
     a = x(2);
     r = x(3);
     mm = x(4);
     mp = x(5);
     mq = x(6);
     mr = x(7);
     mt = x(8);
     rmm = x(9);
     rmp = x(10);
     rmq = x(11);
     rmr = x(12);
     rmt = x(13);
     zmm = x(14);
     zmr = x(15);
     zmp = x(16);
     zmq = x(17);
     zmt = x(18);
     em = x(19);
     et = x(20);
     ep = x(21);
     Q = x(22);

    f = cl*k_cm;
    Kg = gmax/Kp;
    gam = gmax*a/(Kg+a);
    ttrate = (rmq+rmr+rmp+rmt+rmm)*gam;
    lam = ttrate/M;
    Mdyn = (nr*(r+rmr+rmp+rmt+rmm+rmq+zmr+zmp+zmt+zmm+zmq)+nx*(ep+Q+et+em));
    Rmass = (nr*(rmr+zmr)+nx*(mp+mq+mt+mm+rmp+rmt+rmm+rmq+zmp+zmt+zmm+zmq))/3;
    fr = nr*(r+rmr+rmp+rmt+rmm+rmq+zmr+zmp+zmt+zmm+zmq)/Mdyn;
    nuimp = et*vt*s0/(Kt+s0);
    nucat = em*vm*si/(Km+si);

     dx_si = 0 +1*(nuimp) -1*(nucat) -1*(lam*si);%   si
     dx_a = 0 +1*(ns*nucat) -1*(ttrate) -1*(lam*a);%   a
     dx_r = 0 -1*(kb*r*mm) +1*(ku*rmm) -1*(kb*r*mp) +1*(ku*rmp) -1*(kb*r*mq) +1*(ku*rmq) -1*(kb*r*mr) +1*(ku*rmr) -1*(kb*r*mt) +1*(ku*rmt) +1*((gam/nx)*rmm) +1*((gam/nx)*rmp) +1*((gam/nx)*rmq) +2*((gam/nr)*rmr) +1*((gam/nx)*rmt) -1*(lam*r);%   r
     dx_mm = 0 +1*(we*a/(thetax+a)) -1*(kb*r*mm) +1*(ku*rmm) +1*((gam/nx)*rmm) -1*(dm*mm) -1*(lam*mm);%   mm
     dx_mp = 0 +1*(wp*a/(thetax+a)) -1*(kb*r*mp) +1*(ku*rmp) +1*((gam/nx)*rmp) -1*(dm*mp) -1*(lam*mp);%   mp
     dx_mq = 0 +1*((wq*a/(thetax+a))/(1+(Q/kq)^nq)) -1*(kb*r*mq) +1*(ku*rmq) +1*((gam/nx)*rmq) -1*(dm*mq) -1*(lam*mq);%   mq
     dx_mr = 0 +1*(wr*a/(thetar+a)) -1*(kb*r*mr) +1*(ku*rmr) +1*((gam/nr)*rmr) -1*(dm*mr) -1*(lam*mr);%   mr
     dx_mt = 0 +1*(we*a/(thetax+a)) -1*(kb*r*mt) +1*(ku*rmt) +1*((gam/nx)*rmt) -1*(dm*mt) -1*(lam*mt);%   mt
     dx_rmm = 0 +1*(kb*r*mm) -1*(ku*rmm) -1*(f*rmm) +1*(b*zmm) -1*((gam/nx)*rmm) -1*(lam*rmm);%   rmm
     dx_rmp = 0 +1*(kb*r*mp) -1*(ku*rmp) -1*(f*rmp) +1*(b*zmp) -1*((gam/nx)*rmp) -1*(lam*rmp);%   rmp
     dx_rmq = 0 +1*(kb*r*mq) -1*(ku*rmq) -1*(f*rmq) +1*(b*zmq) -1*((gam/nx)*rmq) -1*(lam*rmq);%   rmq
     dx_rmr = 0 +1*(kb*r*mr) -1*(ku*rmr) -1*(f*rmr) +1*(b*zmr) -1*((gam/nr)*rmr) -1*(lam*rmr);%   rmr
     dx_rmt = 0 +1*(kb*r*mt) -1*(ku*rmt) -1*(f*rmt) +1*(b*zmt) -1*((gam/nx)*rmt) -1*(lam*rmt);%   rmt
     dx_zmm = 0 +1*(f*rmm) -1*(b*zmm) -1*(lam*zmm);%   zmm
     dx_zmr = 0 +1*(f*rmr) -1*(b*zmr) -1*(lam*zmr);%   zmr
     dx_zmp = 0 +1*(f*rmp) -1*(b*zmp) -1*(lam*zmp);%   zmp
     dx_zmq = 0 +1*(f*rmq) -1*(b*zmq) -1*(lam*zmq);%   zmq
     dx_zmt = 0 +1*(f*rmt) -1*(b*zmt) -1*(lam*zmt);%   zmt
     dx_em = 0 +1*((gam/nx)*rmm) -1*(lam*em);%   em
     dx_et = 0 +1*((gam/nx)*rmt) -1*(lam*et);%   et
     dx_ep = 0 +1*((gam/nx)*rmp) -1*(lam*ep);%   ep
     dx_Q = 0 +1*((gam/nx)*rmq) -1*(lam*Q);%   Q
  result(asdf) = eval(expr);
  end
end

