/* asl_models.cc Kinetic curve models for ASL

 Michael Chappell - IBME & FMRIB Analysis Group

 Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "asl_models.h"
#include "fabber_core/fwdmodel.h"

#include <algorithm>

#include "fwdmodel_asl_2compartment.h"
#include "fwdmodel_asl_grase.h"
#include "fwdmodel_asl_multiphase.h"
#include "fwdmodel_asl_quasar.h"
#include "fwdmodel_asl_rest.h"
#include "fwdmodel_asl_satrecov.h"
#include "fwdmodel_asl_satrecovdualfa.h"
#include "fwdmodel_asl_turboquasar.h"
#include "fwdmodel_asl_multite.h"
#include "fwdmodel_asl_velocityselective.h"

extern "C" {
int CALL get_num_models() { return 8; }
const char *CALL get_model_name(int index)
{
    switch (index)
    {
    case 0:
        return "aslrest";
        break;
    case 1:
        return "buxton";
        break;
    case 2:
        return "asl_multiphase";
        break;
    case 3:
        return "quasar";
        break;
    case 4:
        return "turboquasar";
        break;
    case 5:
        return "asl_2comp";
        break;
    case 6:
        return "satrecov";
        break;
    case 7:
        return "multite";
        break;
    case 8:
        return "velocity_selective";
        break;
    default:
        return NULL;
    }
}

NewInstanceFptr CALL get_new_instance_func(const char *name)
{
    if (string(name) == "aslrest")
    {
        return ASLFwdModel::NewInstance;
    }
    if (string(name) == "buxton")
    {
        return GraseFwdModel::NewInstance;
    }
    if (string(name) == "asl_multiphase")
    {
        return MultiPhaseASLFwdModel::NewInstance;
    }
    if (string(name) == "quasar")
    {
        return QuasarFwdModel::NewInstance;
    }
    if (string(name) == "turboquasar")
    {
        return TurboQuasarFwdModel::NewInstance;
    }
    if (string(name) == "asl_2comp")
    {
        return ASL2CompartmentModel::NewInstance;
    }
    if (string(name) == "satrecov")
    {
        return SatrecovFwdModel::NewInstance;
    }
    if (string(name) == "multite")
    {
        return multiTEFwdModel::NewInstance;
    }
    if (string(name) == "velocity_selective")
    {
        return VelocitySelectiveFwdModel::NewInstance;
    }
    else
    {
        return NULL;
    }
}
}

namespace OXASL
{
// --- Kinetic curve functions ---

// Wrapper around the error function to screen out very high and low values to avoid 'erfc underflow' errors
static double myerf(double x)
{
    return MISCMATHS::erf(max(min(x, 5.0), -5.0));
}

// Arterial

double AIFModel_nodisp::kcblood(const double ti, const double deltblood, const double taub,
    const double T_1b, const bool casl, const ColumnVector dispparam) const
{
    // Non dispersed arterial curve
    //
    // To avoid problems with the discontinuous gradient at ti=deltblood
    // and ti=deltblood+taub, we smooth the transition at these points
    // using a Gaussian convolved step function. The sigma value could
    // be exposed as a parameter (small value = less smoothing). This is
    // similar to the effect of Gaussian dispersion, but can be computed
    // without numerical integration
    double kcblood;
    if (casl)
        kcblood = 2 * exp(-deltblood / T_1b);
    else
        kcblood = 2 * exp(-ti / T_1b);

    if (ti < (deltblood + taub / 2))
    {
        // If deltblood is smaller than the lead in scale, we could 'lose' some
        // of the bolus, so reduce degree of lead in as deltblood -> 0. We
        // don't really need it in this case anyway since there will be no
        // gradient discontinuity
        double leadscale = min(deltblood, m_leadscale);
        if (leadscale > 0)
        {
            kcblood *= 0.5 * (1 + myerf((ti - deltblood) / leadscale));
        }
        else if (ti < deltblood)
        {
            kcblood = 0;
        }
    }
    else
    {
        // 'lead out' rather than immediate drop to zero
        kcblood *= 0.5 * (1 + myerf(-(ti - deltblood - taub) / m_leadscale));
    }

    return kcblood;
}

// NOTE: for cASL the version here is an over simplificaiton (just changing the
// decay term and leaving the rest alone) since it ignores the fact that some
// blood will be more delayed than the rest due to dispersion
double AIFModel_gammadisp::kcblood(const double ti, const double deltblood, const double taub,
    const double T_1b, const bool casl, const ColumnVector dispparam) const
{
    // Gamma dispersed arterial curve (pASL)
    double kcblood = 0.0;

    // extract dispersion parameters - stored as the log of the true value
    double s;
    double p;
    s = dispparam(1);
    s = exp(s);
    double sp = (dispparam.Row(2)).AsScalar();
    sp = exp(sp);
    if (sp > 10)
        sp = 10;
    p = sp / s;

    double k = 1 + p * s;

    double sprime; // Sharpness parameter modified to include T1 decay effects

    if (casl)
    {
        sprime = s + 1/T_1b;
    }

    if (ti < deltblood)
    {
        kcblood = 0.0;
    }
    else if (ti >= deltblood && ti <= (deltblood + taub))
    {
        if (casl)
        {
            kcblood = 2 * exp(-deltblood / T_1b) * pow(s/sprime,k); // Original T1 decay term modulated by additional dispersion term
            kcblood *= (1 - igamc(k, sprime * (ti - deltblood))); // Modified incomplete gamma function call to account for dispersion-T1 interaction
        }
        else // PASL
        {
            kcblood = 2 * exp(-ti / T_1b);
            // This part is specific to the gamma dispersion model
            kcblood *= (1 - igamc(k, s * (ti - deltblood)));
        }
        
    }
    else //(ti > deltblood + taub)
    {
        if (casl)
        {
            kcblood = 2 * exp(-deltblood / T_1b) * pow(s/sprime,k); // Original T1 decay term modulated by additional dispersion term
            kcblood *= (igamc(k, sprime * (ti - deltblood - taub)) - igamc(k, sprime * (ti - deltblood))); // Modified incomplete gamma function call to account for dispersion-T1 interaction
        }
        else // PASL
        {
            kcblood = 2 * exp(-ti / T_1b);
            // This part is specific to the gamma dispersion model
            kcblood *= (igamc(k, s * (ti - deltblood - taub)) - igamc(k, s * (ti - deltblood)));    
        }        
    }

    return kcblood;
}

double AIFModel_gvf::kcblood(const double ti, const double deltblood, const double taub,
    const double T_1b, const bool casl, const ColumnVector dispparam) const
{
    // GVF AIF shape
    double kcblood = 0.0;

    // extract dispersion parameters
    double s;
    double p;
    s = (dispparam.Row(1)).AsScalar();
    s = exp(s);
    double sp = (dispparam.Row(2)).AsScalar();
    sp = exp(sp);
    if (sp > 10)
        sp = 10;
    p = sp / s;

    // gamma variate arterial curve
    // NOTE:    taub does not directly affect the shape, jsut scale of this KC -
    // see below
    //          NOT a good idea to use when inferring bolus duration.

    if (ti < deltblood)
    {
        kcblood = 0.0;
    }
    else // if(ti >= deltblood) && ti <= (deltblood + taub))
    {
        if (casl)
            kcblood = 2 * exp(-deltblood / T_1b);
        else
            kcblood = 2 * exp(-ti / T_1b);

        kcblood *= gvf(ti - deltblood, s, p);
    }
    // we do not have bolus duration within the GVF AIF - the duration is
    // 'built' into the function shape

    kcblood /= taub; // the 'original' bolus duration scales the magtiude of the
                     // KC becuase the area under the KC is preserved under
                     // dispersion (apart from T1 decay)
    return kcblood;
}

// NOTE: for cASL the version here is an over simplificaiton (just changing the
// decay term and leaving the rest alone) since it ignores the fact that some
// blood will be more delayed than the rest due to dispersion
double AIFModel_gaussdisp::kcblood(const double ti, const double deltblood, const double taub,
    const double T_1b, const bool casl, const ColumnVector dispparam) const
{
    // Gaussian dispersion arterial curve
    // after Hrabe & Lewis, MRM, 2004 for pASL
    // my derrivation based on the spatial gauss solution for cASL
    double kcblood = 0.0;
    double sqrt2 = sqrt(2);

    double sig1;
    sig1 = exp((dispparam.Row(1)).AsScalar());
    double erf1;
    double erf2;

    if (casl)
    {
        double a = 1 / (sig1 * sig1);
        double b = 1 / T_1b;

        double Q = (2 * a * (deltblood - ti) - b) / (2 * a);
        double S = (a * (deltblood - ti) * (deltblood - ti) + b * ti) / a;

        kcblood = 2 * exp(-deltblood / T_1b) * exp(Q * Q - S);

        erf1 = sqrt(a) * (taub + Q);
        erf2 = sqrt(a) * Q;
    }
    else
    {
        // we will only have diserpsion SD (leading edge)
        // assume trailing edge is related to leading edge as Hrabe did
        double sig2;
        if (deltblood > 0)
        {
            sig2 = sig1 * sqrt((deltblood + taub) / deltblood);
        }
        else
            sig2 = sig1;

        kcblood = 2 * exp(-ti / T_1b);

        erf1 = (ti - deltblood) / (sqrt2 * sig1);
        erf2 = (ti - deltblood - taub) / (sqrt2 * sig2);
    }

    kcblood *= 0.5 * (myerf(erf1) - myerf(erf2));
    return kcblood;
}

double AIFModel_spatialgaussdisp_alternate::integral(double t, double k, double A, double B, double C) const
{
    // This implements the indefinite integral from wolfram
    double erf1 = (C - A*t)/sqrt(t);
    double erf2 = (C + A*t)/sqrt(t);
    double P1 = -myerf(erf1);
    double P2 = myerf(erf2) - 1;
    // Check P2 because when P2 -> 0 4AC -> inf and we get overflow
    if (P2 > 1e-6)
    {
        P2 *= exp(4*A*C);
    }

    return exp(-2*A*C - B) / (2*A*k) * (P1 + P2 + 1);
}

double AIFModel_spatialgaussdisp_alternate::kcblood(const double ti, const double deltblood,
    const double taub, const double T_1b, const bool casl, const ColumnVector dispparam) const
{
    // Gaussian dispersion arterial curve - in spatial rather than temporal domain
    // after Ozyurt ISMRM 2010 (p4065) for pASL
    // Using derivation from Thijs van Osch for cASL, but developed into closed
    // form solution
    double k = exp((dispparam.Row(1)).AsScalar());

    if (casl)
    {
        double A = sqrt(1 / (k*k) + (1/T_1b));
        double B = -2*deltblood / (k*k);
        double C = deltblood / k;
        return 2 * (integral(max(1e-6, ti), k, A, B, C) - integral(max(1e-6, ti-taub), k, A, B, C));
    }
    else
    {
        double kcblood = 0.5 * exp(-ti / T_1b);

        double erf1 = (ti - deltblood) / (k * sqrt(ti));
        double erf2 = (ti - deltblood - taub) / (k * sqrt(ti));
        return 2 * kcblood * (myerf(erf1) - myerf(erf2));
    }
}

double AIFModel_spatialgaussdisp::kcblood(const double ti, const double deltblood,
    const double taub, const double T_1b, const bool casl, const ColumnVector dispparam) const
{
    // Gaussian dispersion arterial curve - in spatial rather than temporal domain
    // after Ozyurt ISMRM 2010 (p4065) for pASL
    // using derrivation from Thijs van Osch for cASL

    double kcblood = 0.0;

    double k;
    k = exp((dispparam.Row(1)).AsScalar());
    double erf1;
    double erf2;

    if (casl)
    {
        double dt = 0.01;
        int ndels = floor((ti - max(0.0, ti - taub)) / dt) + 1;
        ColumnVector integrand(ndels);
        double lambda;
        for (int i = 1; i <= ndels; i++)
        {
            lambda = ti - (i - 1) * dt;
            if (lambda < 1e-12)
            {
                integrand(i) = 0;
            }
            else
            {
                integrand(i) = 1 / sqrt(lambda)
                    * exp(-1 / (k * k) * ((deltblood - lambda) * (deltblood - lambda) / lambda))
                    * exp(-lambda / T_1b);
            }
        }

        // final bit
        lambda = max(1e-12, ti - taub);
        double finalintegrand;
        finalintegrand = 1 / sqrt(lambda)
            * exp(-1 / (k * k) * ((deltblood - lambda) * (deltblood - lambda) / lambda))
            * exp(-lambda / T_1b);
        double finaldel = (ti - (ndels - 1) * dt) - lambda;

        double integral;
        integral = numerical_integration(integrand, dt, finalintegrand, finaldel, "trapezium");

        kcblood = 2 * integral / sqrt(M_PI) / k;
    }
    else
    {
        kcblood = 2 * exp(-ti / T_1b);

        erf1 = (ti - deltblood) / (k * sqrt(ti));
        erf2 = (ti - deltblood - taub) / (k * sqrt(ti));

        kcblood *= 0.5 * (myerf(erf1) - myerf(erf2));
    }

    return kcblood;
}

/*
 double kcblood_gallichan(const double ti,const double deltblood,const double
 taub,const double T_1b,const double xdivVm,const bool casl=false) {
 // Model of dispersion based on a geometrical argument from Gallichan MRM 2008
 // Taking equation [6] (so not including QUIPSSII style saturation)
 // including an 'extra' arrival time term as per the paper
 // bolus duration (taub) takes the place of X/V_m and we let it be a variable
 double kcblood = 0.0;

 assert(casl==false); // this model is pASL only

 // NOTE: the +xdivVm correction applied to the ti to shift the curve so that
 the
 // delay associated with the dispersion parameter has been removed, thus BAT is
 independent
 // of the dispersion.

 if(ti < deltblood)
 {
 kcblood = 0.0;
 }
 else if((ti >= deltblood) && ti <= (deltblood + taub))
 {
 if (casl) kcblood = 2 * exp(-deltblood/T_1b);
 else	    kcblood = 2 * exp(-ti/T_1b);

 kcblood *= 1 - xdivVm/(ti + xdivVm - deltblood);
 }
 else //(ti > deltblood + taub)
 {
 if (casl) kcblood = 2 * exp(-deltblood/T_1b);
 else	    kcblood = 2 * exp(-ti/T_1b);

 kcblood *= taub/(ti + xdivVm - deltblood) ;
 }
 return kcblood;
 }
 */
//-----------------------------------------
// Residue functions (these are primiarly specified for use with the numerical
// tissue model)
double ResidModel_wellmix::resid(const double ti, const double fcalib, const double T_1,
    const double T_1b, const double lambda, const ColumnVector residparam) const
{
    // Well mixed single compartment
    // Buxton (1998) model
    double T_1app = 1 / (1 / T_1 + fcalib / lambda);
    return exp(-ti / T_1app);
}

double ResidModel_simple::resid(const double ti, const double fcalib, const double T_1,
    const double T_1b, const double lambda, const ColumnVector residparam) const
{
    // Simple impermeable comparment
    // decays with T1b
    return exp(-ti / T_1b);
}

double ResidModel_imperm::resid(const double ti, const double fcalib, const double T_1,
    const double T_1b, const double lambda, const ColumnVector residparam) const
{
    // impermeable compartment with transit time
    // decays with T1b
    double transit = (residparam.Row(1)).AsScalar();
    double resid = exp(-ti / T_1b);
    if (ti > transit)
        resid = 0.0;
    return resid;
}

double ResidModel_twocpt::resid(const double ti, const double fcalib, const double T_1,
    const double T_1b, const double lambda, const ColumnVector residparam) const
{
    // Two compartment model - the simplest form of the two cpt model
    // No backflow from tissue to blood
    // no venous outflow
    // From Parkes & Tofts and also St. Lawrence 2000 - both models are the same
    // under these assumptions
    // extract residue function parameters
    double kw; // exchange rate = PS/vb
    kw = (residparam.Row(1)).AsScalar();

    // calculate the residue function
    double a = kw + 1 / T_1b;
    double b = (kw * T_1 * T_1b) / (kw * T_1 * T_1b + (T_1 - T_1b));
    return b * exp(-ti / T_1) + (1 - b) * exp(-a * ti);
}

double ResidModel_spa::resid(const double ti, const double fcalib, const double T_1,
    const double T_1b, const double lambda, const ColumnVector residparam) const
{
    // Two compartment model - Single Pass Approximation from St. Lawrence
    // (2000)
    // No backflow from tissue to blood
    // label starts to leave the cappilliary after a capilliary transit time
    // extract residue function parameters
    double PS;
    double vb;
    double tauc;
    PS = (residparam.Row(1)).AsScalar();
    vb = (residparam.Row(2)).AsScalar();
    tauc = (residparam.Row(3)).AsScalar();

    // calcualte residue function
    double a = PS / vb + 1 / T_1b;
    double b = (PS * T_1 * T_1b) / (PS * T_1 * T_1b + (T_1 - T_1b) * vb);
    double ER = 1 - exp(-PS / fcalib - (1 / T_1b - 1 / T_1) * tauc);

    if (ti < tauc)
    {
        return b * exp(-ti / T_1) + (1 - b) * exp(-a * ti);
    }
    else
    {
        return b * ER * exp(-ti / T_1);
    }
}

//----------------------------------
// Tissue Model
double TissueModel_nodisp_simple::kctissue(const double ti, const double fcalib,
    const double delttiss, const double tau, const double T_1b, const double T_1,
    const double lambda, const bool casl, const ColumnVector dispparam,
    const ColumnVector residparam)
{
    // Tissue kinetic curve - well mixed, but no outflow and decay with T1 blood only
    // (This is just the impermeable model with infinite residence time)
    double kctissue = 0.0;

    if (ti < delttiss)
    {
        kctissue = 0;
    }
    else if (ti >= delttiss && ti <= (delttiss + tau))
    {
        if (casl)
            kctissue = exp(-delttiss / T_1b) - exp(-ti / T_1b);
        else
            kctissue = ti - delttiss;
    }
    else //(ti > delttiss + tau)
    {
        if (casl)
            kctissue = exp(-ti / T_1b) * (exp(tau / T_1b) - 1);
        else
            kctissue = tau;
    }

    if (casl)
        kctissue *= T_1b;
    else
        kctissue *= exp(-ti / T_1b);
    kctissue *= 2;

    return kctissue;
}

double TissueModel_nodisp_wellmix::kctissue(const double ti, const double fcalib,
    const double delttiss, const double tau, const double T_1b, const double T_1,
    const double lambda, const bool casl, const ColumnVector dispparam,
    const ColumnVector residparam)
{
    // Tissue kinetic curve no dispersion
    // Buxton (1998) model
    double kctissue = 0.0;
    double T_1app = 1 / (1 / T_1 + fcalib / lambda);
    double R = 1 / T_1app - 1 / T_1b;
    double F = 2 * exp(-ti / T_1app);

    if (ti < delttiss)
    {
        kctissue = 0;
    }
    else if (ti >= delttiss && ti <= (delttiss + tau))
    {
        if (casl)
            kctissue = 2 * T_1app * exp(-delttiss / T_1b) * (1 - exp(-(ti - delttiss) / T_1app));
        else
            kctissue = F / R * ((exp(R * ti) - exp(R * delttiss)));
    }
    else //(ti > delttiss + tau)
    {
        if (casl)
            kctissue = 2 * T_1app * exp(-delttiss / T_1b) * exp(-(ti - tau - delttiss) / T_1app)
                * (1 - exp(-tau / T_1app));
        else
            kctissue = F / R * ((exp(R * (delttiss + tau)) - exp(R * delttiss)));
    }
    return kctissue;
}

double TissueModel_nodisp_imperm::kctissue(const double ti, const double fcalib,
    const double delttiss, const double tau, const double T_1b, const double T_1,
    const double lambda, const bool casl, const ColumnVector dispparam,
    const ColumnVector residparam)
{
    // Tissue kinetic curve no dispersion impermeable vessel
    double kctissue = 0.0;

    // extract the pre-cap residence time
    double taup = residparam.AsScalar();

    if (ti > delttiss && ti < delttiss + taup + tau)
    {
        if (ti < delttiss + tau && ti < delttiss + taup)
        {
            if (casl)
                kctissue = exp(-delttiss / T_1b) - exp(-ti / T_1b);
            else
                kctissue = ti - delttiss;
        }
        else if (ti < delttiss + tau && ti >= delttiss + taup)
        {
            if (casl)
                kctissue = exp(-delttiss / T_1b) * (1 - exp(-taup / T_1b));
            else
                kctissue = taup;
        }
        else if (ti >= delttiss + tau && ti >= delttiss + taup)
        {
            if (casl)
                kctissue = exp(-(ti - tau) / T_1b) - exp(-(delttiss + taup) / T_1b);
            else
                kctissue = delttiss + tau + taup - ti;
        }
        else if (ti >= delttiss + tau && ti < delttiss + taup)
        {
            if (casl)
                kctissue = exp(-ti / T_1b) * (exp(tau / T_1b) - 1);
            else
                kctissue = tau;
        }
    }

    if (casl)
        kctissue *= T_1b;
    else
        kctissue *= exp(-ti / T_1b);
    kctissue *= 2;

    return kctissue;
}

double TissueModel_nodisp_2cpt::kctissue(const double ti, const double fcalib,
    const double delttiss, const double tau, const double T_1b, const double T_1,
    const double lambda, const bool casl, const ColumnVector dispparam,
    const ColumnVector residparam)
{
    // Two compartment model - the simplest form of the two cpt model
    // No backflow from tissue to blood
    // no venous outflow
    // From Parkes & Tofts and also St. Lawrence 2000 - both models are the same
    // under these assumptions
    // extract residue function parameters
    double kw; // exchange rate = PS/vb
    kw = (residparam.Row(1)).AsScalar();

    // calculate the residue function
    double a, b, S, T;
    if (casl)
    {
        // For CASL, the following parameter substitutions makes the equations
        // for kctissue below equivalent to equations [22] and [23] in Parkes/Tofts 2002
        // taking J->T, A->kw, C->S and D->1/T1b and also assuming T1e = T1
        switch (m_solution)
        {
        case SLOW:
            T = kw + 1 / T_1b;
            break;
        case FAST:
            T = kw + 1 / T_1b + 1 / residparam(2);
            break;
        case DIST:
        default:
            T = kw + 1 / T_1b;
            if (ti >= residparam(2))
                T += 1 / residparam(2);
        }
        S = 1 / T_1;
        a = T;
        b = kw / (T - S);
    }
    else
    {
        S = 1 / T_1 - 1 / T_1b;
        T = kw;
        a = kw + 1 / T_1b;
        b = (kw * T_1 * T_1b) / (kw * T_1 * T_1b + (T_1 - T_1b));
    }

    double kctissue;
    if (ti <= delttiss)
    {
        kctissue = 0;
    }
    else if (ti <= (delttiss + tau))
    {
        kctissue = 2 * (b / S * exp(-ti / T_1) * (exp(S * ti) - exp(S * delttiss))
                           + (1 - b) / T * exp(-a * ti) * (exp(T * ti) - exp(T * delttiss)));
    }
    else
    {
        kctissue = 2 * (b / S * exp(-ti / T_1) * (exp(S * (delttiss + tau)) - exp(S * delttiss))
                           + (1 - b) / T * exp(-a * ti)
                               * (exp(T * (delttiss + tau)) - exp(T * delttiss)));
    }

    if (casl)
        kctissue *= exp(-delttiss / T_1b);

    return kctissue;
}

double TissueModel_nodisp_spa::kctissue(const double ti, const double fcalib, const double delttiss,
    const double tau, const double T_1b, const double T_1, const double lambda, const bool casl,
    const ColumnVector dispparam, const ColumnVector residparam)
{
    // Two compartment model
    // No backflow from tissue to blood
    // venous outflow (alhtough we dont model a venous component to the signal here)
    // St. Lawrence 2000
    assert(!casl);

    // extract residue function parameters
    double PS;
    double vb;
    double tauc;
    PS = (residparam.Row(1)).AsScalar();
    vb = (residparam.Row(2)).AsScalar(); // NB for SPA on the whole vb and PS
                                         // appear together (as kw), but PS is
                                         // on its own in ER
    tauc = (residparam.Row(3)).AsScalar();

    // calcualte residue function

    double kctissue = 0.0;
    if (ti > delttiss)
    {
        if (ti < delttiss + tauc && ti < delttiss + tau)
        {
            kctissue = Q(delttiss, ti, ti, PS, vb, tauc, fcalib, T_1, T_1b);
        }
        else if (ti > delttiss + tauc && ti < delttiss + tau)
        {
            kctissue = Q(delttiss, delttiss + tauc, ti, PS, vb, tauc, fcalib, T_1, T_1b)
                + R(delttiss + tauc, ti, ti, PS, vb, tauc, fcalib, T_1, T_1b);
        }
        else if (ti < delttiss + tauc && ti > delttiss + tau)
        {
            kctissue = Q(delttiss, delttiss + tau, ti, PS, vb, tauc, fcalib, T_1, T_1b);
        }
        else if (ti >= delttiss + tauc && ti >= delttiss + tau && ti < delttiss + tau + tauc)
        {
            if (tauc <= tau)
            {
                kctissue = Q(delttiss, delttiss + tauc, ti, PS, vb, tauc, fcalib, T_1, T_1b)
                    + R(delttiss + tauc, delttiss + tau, ti, PS, vb, tauc, fcalib, T_1, T_1b);
            }
            else if (tau < tauc)
            {
                kctissue = Q(delttiss, delttiss + tau, ti, PS, vb, tauc, fcalib, T_1, T_1b);
            }
        }
        else if (ti >= delttiss + tau + tauc)
        {
            kctissue = R(delttiss, delttiss + tau, ti, PS, vb, tauc, fcalib, T_1, T_1b);
        }
    }
    return kctissue;
}

double TissueModel_nodisp_spa::Q(const double t1, const double t2, const double t3, const double PS,
    const double vb, const double tauc, const double fcalib, const double T_1,
    const double T_1b) const
{
    double a = PS / vb + 1 / T_1b;
    double b = (PS * T_1 * T_1b) / (PS * T_1 * T_1b + (T_1 - T_1b) * vb);
    double S = 1 / T_1 - 1 / T_1b;
    double T = PS / vb;
    return b / S * exp(-t3 / T_1) * (exp(S * t2) - exp(S * t1))
        + (1 - b) / T * exp(-a * t3) * (exp(T * t2) - exp(T * t1));
}

double TissueModel_nodisp_spa::R(const double t1, const double t2, const double t3, const double PS,
    const double vb, const double tauc, const double fcalib, const double T_1,
    const double T_1b) const
{
    double b = (PS * T_1 * T_1b) / (PS * T_1 * T_1b + (T_1 - T_1b) * vb);
    double ER = 1 - exp(-PS / fcalib - (1 / T_1b - 1 / T_1) * tauc);
    double S = 1 / T_1 - 1 / T_1b;
    return b * ER / S * exp(-t3 / T_1) * (exp(S * t2) - exp(S * t1));
}

double TissueModel_gammadisp_wellmix::kctissue(const double ti, const double fcalib,
    const double delttiss, const double tau, const double T_1b, const double T_1,
    const double lambda, const bool casl, const ColumnVector dispparam,
    const ColumnVector residparam)
{
    double kctissue = 0.0;

    assert(!casl); // only pASL at the moment!

    // extract dispersion parameters
    double s;
    double p;
    s = (dispparam.Row(1)).AsScalar();
    s = exp(s);
    double sp = (dispparam.Row(2)).AsScalar();
    sp = exp(sp);
    if (sp > 10)
        sp = 10;
    p = sp / s;

    double k = 1 + p * s;
    double T_1app = 1 / (1 / T_1 + fcalib / lambda);
    double A = T_1app - T_1b;
    double B = A + s * T_1app * T_1b;
    if (B < 1e-12)
        B = 1e-12; // really shouldn't happen, but combination of parameters may
                   // arise in artefactual voxels?
    double C = pow(s - 1 / T_1app + 1 / T_1b, p * s);
    if (s - 1 / T_1app + 1 / T_1b <= 0)
        C = 1e-12; // really shouldn't happen, but combination of parameters may
                   // arise in artefactual voxels?

    if (ti < delttiss)
    {
        kctissue = 0;
    }
    else if (ti >= delttiss && ti <= (delttiss + tau))
    {
        kctissue = 2 * 1 / A * exp(-(T_1app * delttiss + (T_1app + T_1b) * ti) / (T_1app * T_1b))
            * T_1app * T_1b * pow(B, -k)
            * (exp(delttiss / T_1app + ti / T_1b) * pow(s * T_1app * T_1b, k)
                           * (1 - igamc(k, B / (T_1app * T_1b) * (ti - delttiss)))
                       + exp(delttiss / T_1b + ti / T_1app) * pow(B, k)
                           * (-1 + igamc(k, s * (ti - delttiss))));
    }
    else //(ti > delttiss + tau)
    {
        kctissue = 2 * 1 / (A * B)
            * (exp(-A / (T_1app * T_1b) * (delttiss + tau) - ti / T_1app) * T_1app * T_1b / C
                       * (pow(s, k) * T_1app * T_1b
                                 * (-1
                                       + exp((-1 / T_1app + 1 / T_1b) * tau)
                                           * (1 - igamc(k, B / (T_1app * T_1b) * (ti - delttiss)))
                                       + igamc(k, B / (T_1app * T_1b) * (ti - delttiss - tau)))
                             - exp(-A / (T_1app * T_1b) * (ti - delttiss - tau)) * C * B
                                 * (igamc(k, s * (ti - delttiss - tau))
                                       - igamc(k, s * (ti - delttiss)))));
    }
    return kctissue;
}

/*
 double kctissue_gvf(const double ti,const double delttiss,const double tau,
 const double T_1b,const double T_1app,const double s,const double p) {
 // Tissue KC with a GVF AIF
 // NOTE: tau onyl scales the magnitude and does not affacet the overall shape
 in this model (see also kcblood_gvf)
 double kctissue = 0.0;

 double k=1+p*s;
 double A = T_1app - T_1b;
 double B = A + s*T_1app*T_1b;
 double C = pow(s-1/T_1app+1/T_1b,p*s);
 double sps = pow(s,k);

 if(ti < delttiss)
 { kctissue = 0.0;}
 else //if(ti >= delttiss && ti <= (delttiss + tau))
 {
 kctissue = 2* 1/(B*C) * exp(-(ti-delttiss)/T_1app)*sps*T_1app*T_1b * (1 -
 igamc(k,(s-1/T_1app-1/T_1b)*(ti-delttiss)));
 }
 // bolus duraiton is specified by the CVF AIF shape and is not an explicit
 parameter
 //else //(ti > delttiss + tau)
 //	{
 //	  kctissue(it) = exp(-(ti-delttiss-tau)/T_1app) * 2* 1/(B*C) *
 exp(-(delttiss+tau)/T_1app)*sps*T_1app*T_1b * (1 -
 igamc(k,(s-1/T_1app-1/T_1b)*(delttiss+tau)));
 //	}

 kctissue /= tau; // the aif is scaled by the duration of the bolus
 return kctissue;
 }

 double kctissue_gaussdisp(const double ti,const double delttiss,const double
 tau,const double T_1b,const double T_1app,const double sig1,const double sig2)
 {
 // Tissue kinetic curve gaussian dispersion (pASL)
 // Hrabe & Lewis, MRM, 2004
 double kctissue = 0.0;

 double R = 1/T_1app - 1/T_1b;
 double sqrt2 = sqrt(2);
 double F = 2 * exp(-ti/T_1app);
 double u1 = (ti-delttiss)/(sqrt2*sig1);
 double u2 = (ti - delttiss - tau)/(sqrt2*sig2);

 kctissue = F/(2*R) * (  (erf(u1) - erf(u2))*exp(R*ti)
 - (1 + erf(u1 - (R*sig1)/sqrt2))*exp(R*(delttiss+(R*sig1*sig1)/2))
 + (1 + erf(u2 - (R*sig2)/sqrt2))*exp(R*(delttiss+tau+(R*sig2*sig2)/2))
 );

 return kctissue;
 }
 */

/**
 * Calculate Kc using numerical convolution integral of AIF and residue functions
 */
double TissueModel_aif_residue::kctissue(const double ti, const double fcalib,
    const double delttiss, const double tau, const double T_1b, const double T_1,
    const double lambda, const bool casl, const ColumnVector dispparam,
    const ColumnVector residparam)
{
    // Number of time point required, including one at zero
    int nt = floor(ti / m_delta) + 1;

    // Extra time interval required to reach the TI
    double dti = ti - (nt - 1) * m_delta;

    // Calculate the convolution integral of AIF(t)*RESID(TI-t) from 0 to TI
    // using the trapezium rule
    double kctissue = 0;
    double integrand = 0;
    if (nt > 1)
    {
        for (int i = 0; i < nt; i++)
        {
            double t = i * m_delta;
            double aif = aifmodel->kcblood(t, delttiss, tau, T_1b, casl, dispparam);
            double resid = residmodel->resid(ti - t, fcalib, T_1, T_1b, lambda, residparam);
            integrand = aif * resid;
            if ((i == 0) || (i == (nt - 1)))
                kctissue += integrand / 2;
            else
                kctissue += integrand;
        }
        kctissue *= m_delta;
    }
    // Last bit of integral to get from nt * m_delta to TI. Note that resid(0) = 1
    double final_integrand = aifmodel->kcblood(ti, delttiss, tau, T_1b, casl, dispparam);
    kctissue += (0.5 * final_integrand + 0.5 * integrand) * dti;

    return kctissue;
}

// --- useful general functions ---
double icgf(const double a, const double x)
{
    // incomplete gamma function with a=k, based on the incomplete gamma integral

    return MISCMATHS::gamma(a) * igamc(a, x);
}

double gvf(const double t, const double s, const double p)
{
    // The Gamma Variate Function (correctly normalised for area under curve)
    // Form of Rausch 2000
    // NB this is basically a gamma pdf
    if (t < 0)
        return 0.0;
    else
        return pow(s, 1 + s * p) / MISCMATHS::gamma(1 + s * p) * pow(t, s * p) * exp(-s * t);
}

double numerical_integration(
    ColumnVector integrand, double del, double finalval, double finaldel, string method)
{
    int ndel;
    ndel = integrand.Nrows();
    ColumnVector prod(ndel);
    prod = integrand;

    // do numerical intergration on a supplied equispaced vector - with possible
    // extra point at end with spacing: finaldel < del
    if (method == "rect")
    {
    }
    // rect intergration
    // nothing to do here - this is essentially done as the default option
    else if (method == "trapezoid")
    {
        ColumnVector trap(ndel);
        trap = 1;
        trap(1) = 0.5;
        trap(del) = 0.5;
        prod = SP(prod, trap);
    }

    double result;
    result = prod.Sum() * del;
    if (finaldel > 1e-12)
    {
        result += (0.5 * finalval + 0.5 * prod(ndel)) * finaldel;
    }

    return result;
}
} // end namespace
