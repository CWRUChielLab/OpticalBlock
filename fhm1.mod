TITLE FHm1 channel
: Frankenhaeuser - Huxley channels for Xenopus L
: standard temperature for FH model is 20 degC
:
: This version adds local temperature, and is from the code from the
: following paper:
:
: Mou, Zongxia, I.F. Triantis, V.M. Woods, C. Toumazou, and K. Nikolic. 2012.
: "A Simulation Study of the Combined Thermoelectric Extracellular Stimulation
: of the Sciatic Nerve of the Xenopus Laevis: The Localized Transient Heat
: Block." IEEE Transactions on Biomedical Engineering 59 (6): 1758-69.
: doi:10.1109/TBME.2012.2194146.
:
: Copyright (c) 2012, Zongxia Mou and Konstantin Nikolic
: All rights reserved.
:
: Redistribution and use in source and binary forms, with or without
: modification, are permitted provided that the following conditions are
: met:
:
:     * Redistributions of source code must retain the above copyright
:       notice, this list of conditions and the following disclaimer.
:     * Redistributions in binary form must reproduce the above copyright
:       notice, this list of conditions and the following disclaimer in
:       the documentation and/or other materials provided with the distribution
:
: THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
: AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
: IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
: ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
: LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
: CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
: SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
: INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
: CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
: ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
: POSSIBILITY OF SUCH DAMAGE.

NEURON {
    SUFFIX fhm1
    USEION na READ nai, nao WRITE ina
    USEION k READ ki, ko WRITE ik
    NONSPECIFIC_CURRENT il, ip
    RANGE pnabar, pkbar, ppbar, gl, el, il, ip
    GLOBAL inf,tau
        RANGE localtemp
}


UNITS {
    (molar) = (/liter)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {
    v (mV)
    localtemp (degC) : 20
    pnabar=8e-3 (cm/s)
    ppbar=.54e-3 (cm/s)
    pkbar=1.2e-3 (cm/s)
    nai (mM) : 13.74
    nao (mM) : 114.5
    ki (mM) : 120
    ko (mM) : 2.5
    gl=30.3e-3 (mho/cm2)
    el = -69.74 (mV)
}
STATE {
    m h n p
}
ASSIGNED {
    ina (mA/cm2)
    ik (mA/cm2)
    ip (mA/cm2)
    il (mA/cm2)
    inf[4]
    tau[4] (ms)
}

INITIAL {
    mhnp(v*1(/mV))
    m = inf[0]
    h = inf[1]
    n = inf[2]
    p = inf[3]
}

BREAKPOINT {
    LOCAL ghkna, q10ch
    SOLVE states METHOD cnexp
    ghkna = ghk(v, nai, nao)
    q10ch = 1.4^((localtemp - 20)/10)
    ina = q10ch*pnabar*m*m*h*ghkna
    ip = q10ch*ppbar*p*p*ghkna
    ik = q10ch*pkbar*n*n*ghk(v, ki, ko)
    il = q10ch*gl*(v - el)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
    :assume a single charge
    LOCAL z, eci, eco
    z = (1e-3)*FARADAY*v/(R*(localtemp+273.15))
    eco = co*efun(z)
    eci = ci*efun(-z)
    ghk = (.001)*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
        efun = 1 - z/2
    }else{
        efun = z/(exp(z) - 1)
    }
}

DERIVATIVE states {    : exact when v held constant
    mhnp(v*1(/mV))
    m' = (inf[0] - m)/tau[0]
    h' = (inf[1] - h)/tau[1]
    n' = (inf[2] - n)/tau[2]
    p' = (inf[3] - p)/tau[3]
}

UNITSOFF
FUNCTION alp(v,i) { LOCAL a,b,c,q10m,q10h,q10n,q10p :rest = -70  order m,h,n,p
    v = v+70
    if (i==0) {
        q10m=1.8^((localtemp - 20)/10)
           a=.36 b=22. c=3.
        alp = q10m*a*expM1(b - v, c)
    }else if (i==1){
           q10h = 2.8^((localtemp - 20)/10)
        a=.1 b=-10. c=6.
        alp = q10h*a*expM1(v - b, c)
    }else if (i==2){
        q10n = 3.2^((localtemp - 20)/10)
           a=.02 b= 35. c=10.
        alp = q10n*a*expM1(b - v, c)
    }else{
        q10p = 3^((localtemp - 20)/10)
        a=.006 b= 40. c=10.
        alp = q10p*a*expM1(b - v , c)
    }
}

FUNCTION bet(v,i) { LOCAL a,b,c,q10m,q10h,q10n,q10p :rest = -70  order m,h,n,p
    v = v+70
    if (i==0) {
        q10m=1.7^((localtemp - 20)/10)
           a=.4  b= 13.  c=20.
        bet = q10m*a*expM1(v - b, c)
    }else if (i==1){
        q10h = 2.9^((localtemp - 20)/10)
           a=4.5  b= 45.  c=10.
        bet = q10h*a/(exp((b - v)/c) + 1)
    }else if (i==2){
           q10n = 2.8^((localtemp - 20)/10)
           a=.05  b= 10.  c=10.
        bet = q10n*a*expM1(v - b, c)
    }else{
        q10p = 3^((localtemp - 20)/10)
        a=.09 b= -25. c=20.
        bet = q10p*a*expM1(v - b, c)
    }
}

FUNCTION expM1(x,y) {
    if (fabs(x/y) < 1e-6) {
        expM1 = y*(1 - x/y/2)
    }else{
        expM1 = x/(exp(x/y) - 1)
    }
}

PROCEDURE mhnp(v) {LOCAL a, b :rest = -70
    TABLE inf, tau DEPEND localtemp FROM -100 TO 100 WITH 200
    FROM i=0 TO 3 {
        a = alp(v,i)  b=bet(v,i)
        tau[i] = 1/(a + b)
        inf[i] = a/(a + b)
    }
}
UNITSON
