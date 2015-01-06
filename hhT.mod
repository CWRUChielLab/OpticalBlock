TITLE hhT.mod   squid sodium, potassium, and leak channels

COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium,
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
}

? interface
NEURON {
        SUFFIX hhT
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk, m_alpha_q10, m_beta_q10
        RANGE n_alpha_q10, n_beta_q10, h_alpha_q10, h_beta_q10, localtemp
        GLOBAL minf, hinf, ninf, mtau, htau, ntau
        THREADSAFE : assigned GLOBALs will be per thread
}

PARAMETER {
        localtemp = 20 (degC)
        gnabar = .12 (S/cm2)        <0,1e9>
        gkbar = .036 (S/cm2)        <0,1e9>
        gl = .0003 (S/cm2)          <0,1e9>
        el = -54.3 (mV)
        m_alpha_q10 = 1.8
        m_beta_q10 = 1.8
        n_alpha_q10 = 3
        n_beta_q10 = 3
        h_alpha_q10 = 3
        h_beta_q10 = 3
}

STATE {
        m h n
}

ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

        gna (S/cm2)
        gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        minf hinf ninf
        mtau (ms) htau (ms) ntau (ms)
}

? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
        ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
        ik = gk*(v - ek)
        il = gl*(v - el)
}


INITIAL {
        rates(v)
        m = minf
        h = hinf
        n = ninf
}

? states
DERIVATIVE states {
        rates(v)
        m' =  (minf-m)/mtau
        h' = (hinf-h)/htau
        n' = (ninf-n)/ntau
}


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
        LOCAL  m_alpha, m_beta, m_sum, m_alpha_q10_term, m_beta_q10_term, h_alpha, h_beta, h_sum, h_alpha_q10_term, h_beta_q10_term, n_alpha, n_beta, n_sum, n_alpha_q10_term, n_beta_q10_term
        TABLE minf, mtau, hinf, htau, ninf, ntau DEPEND localtemp FROM -100 TO 600 WITH 700

UNITSOFF
        : The magic numbers are experimental values from Hodgekin and Huxley 1952d;
        : these will all need to be updated for difference reference values.
                :"m" sodium activation system
        m_alpha_q10_term = m_alpha_q10^((localtemp - 6.3)/10)
        m_beta_q10_term = m_beta_q10^((localtemp - 6.3)/10)
        m_alpha = m_alpha_q10_term * .1 * vtrap(-(v+40),10)
        m_beta =  m_beta_q10_term * 4 * exp(-(v+65)/18)
        m_sum = m_alpha + m_beta
        mtau = 1/m_sum
        minf = m_alpha/m_sum
                :"h" sodium inactivation system
        h_alpha_q10_term = h_alpha_q10^((localtemp - 6.3)/10)
        h_beta_q10_term = h_beta_q10^((localtemp - 6.3)/10)
        h_alpha = h_alpha_q10_term * .07 * exp(-(v+65)/20)
        h_beta =  h_beta_q10_term * 1 / (exp(-(v+35)/10) + 1)
        h_sum = h_alpha + h_beta
        htau = 1/h_sum
        hinf = h_alpha/h_sum
                :"n" potassium activation system
        n_alpha_q10_term = n_alpha_q10^((localtemp - 6.3)/10)
        n_beta_q10_term = n_beta_q10^((localtemp - 6.3)/10)
        n_alpha = n_alpha_q10_term * .01*vtrap(-(v+55),10)
        n_beta =  n_beta_q10_term * .125*exp(-(v+65)/80)
        n_sum = n_alpha + n_beta
        ntau = 1/n_sum
        ninf = n_alpha/n_sum
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

UNITSON
