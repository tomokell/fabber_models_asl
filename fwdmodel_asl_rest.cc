/*  fwdmodel_asl_rest.cc - Resting state ASL model

    Michael Chappell, QuBIc & FMRIB Image Analysis Group

    Copyright (C) 2011-12 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_asl_rest.h"

#include "miscmaths/miscprob.h"
#include "newimage/newimageall.h"
#include <iostream>
#include <newmatio.h>
#include <stdexcept>
using namespace NEWIMAGE;
#include "fabber_core/easylog.h"
#include "fabber_core/tools.h"

FactoryRegistration<FwdModelFactory, ASLFwdModel> ASLFwdModel::registration("aslrest");

string ASLFwdModel::ModelVersion() const
{
    string version = "fwdmodel_asl_rest.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

static OptionSpec OPTIONS[] = {
    { "disp", OPT_STR, "AIF dispersion type", OPT_NONREQ, "none" },
    { "exch", OPT_STR, "Type of exchange in tissue compartment", OPT_NONREQ, "mix" },
    { "forceconv", OPT_BOOL, "Force numerical convolution for evaluation of model", OPT_NONREQ,
        "" },
    { "conv-delta", OPT_FLOAT,
        "Time interval for numerical convolution (used if disp=none or forceconv specified",
        OPT_NONREQ, "0.1" },
    { "inctiss", OPT_BOOL, "Include tissue parameters", OPT_NONREQ, "" },
    { "infertiss", OPT_BOOL, "Infer tissue parameters", OPT_NONREQ, "" },
    { "incart", OPT_BOOL, "Include arterial parameters", OPT_NONREQ, "" },
    { "inferart", OPT_BOOL, "Infer arterial parameters", OPT_NONREQ, "" },
    { "incwm", OPT_BOOL, "Include white matter parameters", OPT_NONREQ, "" },
    { "incbat", OPT_BOOL, "Include bolus arrival time parameter", OPT_NONREQ, "" },
    { "inferbat", OPT_BOOL, "Infer bolus arrival time parameter", OPT_NONREQ, "" },
    { "auto-init-bat", OPT_BOOL,
        "Automatically initialize BAT posterior (only if inferring BAT and number of TIs > 1)",
        OPT_NONREQ, "" },
    { "incpc", OPT_BOOL, "Include pre-capillary parameters", OPT_NONREQ, "" },
    { "inferpc", OPT_BOOL, "Infer pre-capillary parameters", OPT_NONREQ, "" },
    { "inctau", OPT_BOOL, "Include bolus duration parameter", OPT_NONREQ, "" },
    { "infertau", OPT_BOOL, "Infer bolus duration parameter", OPT_NONREQ, "" },
    { "septau", OPT_BOOL, "Separate values of tau for each component", OPT_NONREQ, "" },
    { "inct1", OPT_BOOL, "Include T1 parameter", OPT_NONREQ, "" },
    { "infert1", OPT_BOOL, "Infer T1 parameter", OPT_NONREQ, "" },
    { "inferdisp", OPT_BOOL, "Infer dispersion parameters (if present in model)", OPT_NONREQ, "" },
    { "sepdisp", OPT_BOOL, "Separate tissue dispersion parameters for each component", OPT_NONREQ,
        "" },
    { "inferexch", OPT_BOOL, "Infer exchange parameters (if present in model)", OPT_NONREQ, "" },
    { "incpve", OPT_BOOL, "Include PVE parameters", OPT_NONREQ, "" },
    { "pvcorr", OPT_BOOL, "Partial volume correction", OPT_NONREQ, "" },
    { "incstattiss", OPT_BOOL, "Include static tissue parameters", OPT_NONREQ, "" },
    { "inferstattiss", OPT_BOOL, "Infer static tissue parameters", OPT_NONREQ, "" },
    { "ardoff", OPT_BOOL, "Turn off ARD", OPT_NONREQ, "" },
    { "repeats", OPT_INT, "Number of repeats in data", OPT_NONREQ, "1" },
    { "pretisat", OPT_FLOAT, "Deal with saturation of the bolus a fixed time pre TI measurement",
        OPT_NONREQ, "0.0" },
    { "slicedt", OPT_FLOAT, "Increase in TI per slice (s)", OPT_NONREQ, "0.0" },
    { "sliceband", OPT_INT, "Number of slices per band for multiband readouts", OPT_NONREQ, "1" },
    { "casl", OPT_BOOL, "Data is CASL (not PASL)", OPT_NONREQ, "PASL" },
    { "bat", OPT_FLOAT, "Bolus arrival time", OPT_NONREQ, "0.7" },
    { "batwm", OPT_FLOAT, "Bolus arrival time (white matter)", OPT_NONREQ, "bat+0.3" },
    { "batart", OPT_FLOAT, "Bolus arrival time (arterial)", OPT_NONREQ, "bat-0.3" },
    { "batsd", OPT_FLOAT, "Bolus arrival time standard deviation", OPT_NONREQ, "0.316" },
    { "bat-init", OPT_FLOAT, "Initial value of BAT to use in optimisation procedure. If not "
                             "specified use prior. Special values: step=use step function "
                             "auto-estimation, max=use max signal auto-estimation",
        OPT_NONREQ, "0.7" },
    { "iaf", OPT_STR, "Data information", OPT_NONREQ, "diff" },
    { "calib", OPT_BOOL, "Data has already been subjected to calibration", OPT_NONREQ, "" },
    { "t1", OPT_FLOAT, "T1 value", OPT_NONREQ, "1.3" },
    { "t1b", OPT_FLOAT, "T1b value", OPT_NONREQ, "1.65" },
    { "t1wm", OPT_FLOAT, "T1wm value", OPT_NONREQ, "1.1" },
    { "lambda", OPT_FLOAT, "lambda value", OPT_NONREQ, "0.9 (0.98 with WM component)" },
    { "lambdawm", OPT_FLOAT, "lambda value for WM", OPT_NONREQ, "0.82" },
    { "ti", OPT_FLOAT, "Single TI value (s)", OPT_NONREQ, "" },
    { "ti<n>", OPT_FLOAT, "List of TI values (s)", OPT_NONREQ, "" },
    { "pld", OPT_FLOAT, "Single PLD value (s)", OPT_NONREQ, "" },
    { "pld<n>", OPT_FLOAT, "List of PLD values (s)", OPT_NONREQ, "" },
    { "tiimg", OPT_TIMESERIES, "4D volume containing the effective TI (s) at each voxel and for each volume", OPT_NONREQ, ""},
    { "hadamard", OPT_INT, "Number of Hadamard encoding lines", OPT_NONREQ, "0" },
    { "fullhad", OPT_BOOL, "Activate full Hadamard matrix", OPT_NONREQ, "" },
    { "tau", OPT_FLOAT, "Single tau value", OPT_NONREQ, "" },
    { "tau<n>", OPT_FLOAT, "List of tau values", OPT_NONREQ, "" },
    { "crush<n>", OPT_STR, "List of vascular crushing specifications values", OPT_NONREQ, "" },
    { "FA", OPT_FLOAT, "Look-Locker correction", OPT_NONREQ, "" },
    { "facorr", OPT_BOOL, "Do FA correction", OPT_NONREQ, "" },
    { "2cpt-solution", OPT_BOOL, "For 2-compartment model, solution type (fast, slow or dist)",
        OPT_NONREQ, "slow" },
    { "mtt", OPT_BOOL,
        "For 2-compartment model and fast/dist solution method, mean transit time (s)", OPT_NONREQ,
        "1" },
    { "capria", OPT_BOOL, "CAPRIA acquisition", OPT_NONREQ, "" }, // T.O. Add CAPRIA options
    { "capriafa1", OPT_FLOAT, "CAPRIA starting flip angle", OPT_NONREQ, "" }, 
    { "capriafa2", OPT_FLOAT, "CAPRIA finishing flip angle", OPT_NONREQ, "" }, 
    { "capriatr", OPT_FLOAT, "CAPRIA TR", OPT_NONREQ, "" }, 
    { "dispsd", OPT_FLOAT, "Dispersion parameter standard deviation", OPT_NONREQ, "" }, // T.O. Add dispersion prior options
    { "" },
};

void ASLFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string ASLFwdModel::GetDescription() const { return "Resting state ASL model"; }
void ASLFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    assert(prior.means.Nrows() == NumParams());
    // Set priors
    SymmetricMatrix precisions = IdentityMatrix(NumParams())
        * 1e12;      // by default all paramerers are included as fully informative
    prior.means = 0; // set a default zero value - this should get overwritten
                     // if a non-zero mean is not applicable
    // Flow
    if (infertiss)
    {
        prior.means(flow_index()) = 0;
        precisions(flow_index(), flow_index()) = 1e-12;
    }
    if (inferwm)
    {
        prior.means(flow_index() + wmidx) = 0;
        precisions(flow_index() + wmidx, flow_index() + wmidx) = 1e-12;
    }
    if (inferart)
    {
        prior.means(flow_index() + artidx) = 0;
        precisions(flow_index() + artidx, flow_index() + artidx) = 1e-12;
    }
    // BAT
    if (incbat)
    {
        if (inctiss)
        {
            prior.means(bat_index()) = setdelt;
            if (infertiss & inferbat)
            {
                precisions(bat_index(), bat_index()) = deltprec;
            }
        }
        if (incwm)
        {
            prior.means(bat_index() + wmidx) = setdeltwm;
            if (inferwm & inferbat)
            {
                precisions(bat_index() + wmidx, bat_index() + wmidx) = deltprec;
            }
        }
        if (incart)
        {
            prior.means(bat_index() + artidx) = setdeltart;
            if (inferart & inferbat)
            {
                precisions(bat_index() + artidx, bat_index() + artidx) = deltartprec;
            }
        }
    }
    // tau
    if (inctau)
    {
        if (septau)
        {
            if (inctiss)
            {
                prior.means(tau_index()) = seqtau;
                if (infertau)
                {
                    precisions(tau_index(), tau_index()) = 10;
                }
            }
            if (incwm)
            {
                prior.means(tau_index() + wmidx) = seqtau;
                if (infertau)
                {
                    precisions(tau_index() + wmidx, tau_index() + wmidx) = 10;
                }
            }
            if (incart)
            {
                prior.means(tau_index() + artidx) = seqtau;
                if (infertau)
                {
                    precisions(tau_index() + artidx, tau_index() + artidx) = 10;
                }
            }
        }
        else
        {
            prior.means(tau_index()) = seqtau;
            if (infertau)
            {
                precisions(tau_index(), tau_index()) = 10;
            }
        }
    }
    // T1
    if (inct1)
    {
        if (inctiss)
        {
            prior.means(t1_index()) = t1;
            if (infert1)
            {
                precisions(t1_index(), t1_index()) = 100;
            }
        }
        if (incwm)
        {
            prior.means(t1_index() + wmidx) = t1wm;
            if (infert1)
            {
                precisions(t1_index() + wmidx, t1_index() + wmidx) = 100;
            }
        }
        prior.means(t1_index() + artidx)
            = t1b; // always include t1b, the artidx is still okay to use here
                   // even if we dont have an explicit arterial component
        if (infert1)
        {
            precisions(t1_index() + artidx, t1_index() + artidx) = 100;
        }
    }
    // PVE
    if (incpve)
    {
        // PVE if used are set as image priors (so this is overwridden)
        // PVE if only included should be pure GM
        prior.means(pv_index()) = 1.0;
        prior.means(pv_index() + 1) = 0.0;
    }
    // taupc
    if (incpc)
    {
        if (inctiss)
        {
            prior.means(taupc_index()) = 0.3;
            if (inferpc)
            {
                precisions(taupc_index(), taupc_index()) = 10;
            }
        }
        if (incwm)
        {
            prior.means(taupc_index() + wmidx) = 0.5; // WM longer than GM
            if (inferpc)
            {
                precisions(taupc_index() + wmidx, taupc_index() + wmidx) = 10;
            }
        }
    }
    // dispersion
    if (incdisp)
    {
        ColumnVector disppriors;
        ColumnVector dispmean;
        ColumnVector dispprec;
        if (sepdisp)
        {
            if (inctiss)
            {
                // load priors from tissue model
                int ndisp;
                ndisp = tiss_model->NumDisp();
                disppriors = tiss_model->DispPriors();
                dispmean.ReSize(ndisp);
                dispmean = disppriors.Rows(1, ndisp);
                dispprec.ReSize(ndisp);
                dispprec = disppriors.Rows(ndisp + 1, 2 * ndisp);
                for (int i = 1; i <= ndisp; i++)
                {
                    prior.means(disp_index() + i - 1) = dispmean(i);
                }
                if (inferdisp)
                {
                    for (int i = 1; i <= ndisp; i++)
                    {
                        if (dispprecset > 0)
                        { // Prior SD set by the user, so use this
                            LOG << "Updating tissue dispersion parameter " << i << " prior precision to " << dispprecset << endl;
                            precisions(disp_index() + i - 1, disp_index() + i - 1) = dispprecset;
                        } else { // No (sensible) prior SD set by the user, so use the defaults
                            precisions(disp_index() + i - 1, disp_index() + i - 1) = dispprec(i);
                        }
                    }
                }
            }
            if (incwm)
            {
                // load priors from tissue model
                int ndisp;
                ndisp = tiss_model->NumDisp();
                disppriors = tiss_model->DispPriors();
                dispmean.ReSize(ndisp);
                dispmean = disppriors.Rows(1, ndisp);
                dispprec.ReSize(ndisp);
                dispprec = disppriors.Rows(ndisp + 1, 2 * ndisp);
                for (int i = 1; i <= ndisp; i++)
                {
                    prior.means(disp_index() + ndisp * wmidx + i - 1) = dispmean(i);
                }
                if (inferdisp)
                {
                    for (int i = 1; i <= ndisp; i++)
                    {
                        if (dispprecset > 0)
                        { // Prior SD set by the user, so use this
                            LOG << "Updating WM dispersion parameter " << i << " prior precision to " << dispprecset << endl;
                            precisions(disp_index() + ndisp * wmidx + i - 1,
                            disp_index() + ndisp * wmidx + i - 1)
                            = dispprecset;
                        } else { // Use the defaults
                            precisions(disp_index() + ndisp * wmidx + i - 1,
                                disp_index() + ndisp * wmidx + i - 1)
                                = dispprec(i);
                        }
                    }
                }
            }
            if (incart)
            {
                // load priors from art model
                int ndisp;
                ndisp = art_model->NumDisp();
                dispmean.ReSize(ndisp);
                dispmean = disppriors.Rows(1, ndisp);
                dispprec.ReSize(ndisp);
                dispprec = disppriors.Rows(ndisp + 1, 2 * ndisp);

                for (int i = 1; i <= ndisp; i++)
                {
                    prior.means(disp_index() + ndisp * artidx + i - 1) = dispmean(i);
                }
                if (inferdisp)
                {
                    for (int i = 1; i <= ndisp; i++)
                    {

                        if (dispprecset > 0)
                        { // Prior SD set by the user, so use this
                            LOG << "Updating arterial dispersion parameter " << i << " prior precision to " << dispprecset << endl;
                            precisions(disp_index() + ndisp * artidx + i - 1,
                                disp_index() + ndisp * artidx + i - 1)
                                = dispprecset;
                        } else { // No (sensible) prior SD set by the user, so use the defaults
                            precisions(disp_index() + ndisp * artidx + i - 1,
                                disp_index() + ndisp * artidx + i - 1)
                                = dispprec(i);
                        }
                    }
                }
            }
        }
        else
        {
            // dispersion parameters are shared between arterial and tissue
            // models
            // load the defaults from the tissue model
            int ndisp;
            ndisp = tiss_model->NumDisp();
            disppriors = art_model->Priors();
            // cout << "DISPERSION PRIORS" << endl << disppriors << endl;
            disppriors = tiss_model->DispPriors();
            dispmean.ReSize(ndisp);
            dispmean = disppriors.Rows(1, ndisp);
            dispprec.ReSize(ndisp);
            dispprec = disppriors.Rows(ndisp + 1, 2 * ndisp);

            for (int i = 1; i <= ndisp; i++)
            {
                prior.means(disp_index() + i - 1) = dispmean(i);
            }
            if (inferdisp)
            {
                for (int i = 1; i <= ndisp; i++)
                {
                    if (dispprecset > 0)
                    { // Prior SD set by the user, so use this
                        LOG << "Updating combined arterial/tissue dispersion parameter " << i << " prior precision to " << dispprecset << endl;
                        precisions(disp_index() + i - 1, disp_index() + i - 1) = dispprecset;
                    } else { // No (sensible) prior SD set by the user, so use the defaults
                        precisions(disp_index() + i - 1, disp_index() + i - 1) = dispprec(i);
                    }
                }
            }
        }
    }
    // residue function (exchange) parameters
    if (incexch)
    {
        ColumnVector residpriors;
        ColumnVector residmean;
        ColumnVector residprec;
        int nresid;
        nresid = tiss_model->NumResid();
        residpriors = tiss_model->ResidPriors();
        residmean.ReSize(nresid);
        residmean = residpriors.Rows(1, nresid);
        residprec.ReSize(nresid);
        residprec = residpriors.Rows(nresid + 1, 2 * nresid);
        for (int i = 1; i <= nresid; i++)
        {
            prior.means(resid_index() + i - 1) = residmean(i);
        }
        if (inferexch)
        {
            for (int i = 1; i <= nresid; i++)
            {
                precisions(resid_index() + i - 1, resid_index() + i - 1) = residprec(i);
            }
        }
    }
    if (incfacorr)
    {
        prior.means(facorr_index()) = 1;                  // should be overwritten by an image
        precisions(facorr_index(), facorr_index()) = 100; // small uncertainty
    }
    // Static tissue contribution to the signal (e.g. non-subtracted data)
    if (inferstattiss)
    {
        prior.means(stattiss_index()) = 0;
        precisions(stattiss_index(), stattiss_index()) = 1e-12;
    }
    // Set precsions on priors
    prior.SetPrecisions(precisions);

    // Set initial posterior
    posterior = prior;
    // For parameters with uniformative prior chosoe more sensible inital
    // posterior
    // Tissue perfusion
    if (infertiss)
    {
        posterior.means(flow_index()) = 10;
        precisions(flow_index(), flow_index()) = 1;
    }
    // Arterial perfusion
    if (inferart)
    {
        posterior.means(flow_index() + artidx) = 10;
        precisions(flow_index() + artidx, flow_index() + artidx) = 1;
    }
    if (inferstattiss)
    {
        posterior.means(stattiss_index()) = 1000;
        precisions(stattiss_index(), stattiss_index()) = 1e2;
    }
    posterior.SetPrecisions(precisions);
}

/**
 * Fit a step function to data vector
 *
 * @return the data index (starting at 0) of the data point after which the data
 *         increases in a step-wise fashion giving the best fit
 */
static int fit_step(const ColumnVector &data)
{
    // MSC method - fit a step function
    double best_ssq = -1;
    int step_pos = 0;
    for (int pos = 1; pos < data.Nrows(); pos++)
    {
        double mean_left = data.Rows(1, pos).Sum() / pos;
        double mean_right = data.Rows(pos + 1, data.Nrows()).Sum() / (data.Nrows() - pos);
        if (mean_right > mean_left)
        {
            double ssq = 0;
            for (int t = 1; t <= pos; t++)
            {
                ssq += (data(t) - mean_left) * (data(t) - mean_left);
            }
            for (int t = pos + 1; t <= data.Nrows(); t++)
            {
                ssq += (data(t) - mean_right) * (data(t) - mean_right);
            }
            if ((ssq < best_ssq) || (best_ssq < 0))
            {
                best_ssq = ssq;
                step_pos = pos - 1;
            }
        }
    }

    return step_pos;
}

void ASLFwdModel::InitParams(MVNDist &posterior) const
{
    if (inferstattiss)
    {
        // if we have static tissue then this should be init from the raw data
        // intensities
        double dataval = data.Maximum();
        posterior.means(stattiss_index()) = dataval;

        if (infertiss)
        {
            // init the CBF  - 1% of static tissue
            posterior.means(flow_index()) = 0.01 * dataval;
        }

        if (inferart)
        {
            // init the aBV - 1% of static tissue
            posterior.means(flow_index() + artidx) = 0.01 * dataval;
        }
    }
    else
    {
        if (infertiss)
        {
            // init the CBF  - to max value in the data
            posterior.means(flow_index()) = data.Maximum();
        }

        if (inferart)
        {
            // init the aBV - use max value in data
            posterior.means(flow_index() + artidx) = data.Maximum();
        }
        if (inferbat)
        {
            if (bat_init == "step")
            {
                // MSC method: Average data across TIs and fit a step function to estimate BAT
                ColumnVector data_mean(tis.Nrows());
                data_mean = 0;
                int ti_start = 1;
                for (int ti = 0; ti < tis.Nrows(); ti++)
                {
                    for (int rpt = 0; rpt < repeats[ti]; rpt++)
                    {
                        data_mean(ti + 1) += data(ti_start + rpt);
                    }
                    data_mean(ti + 1) /= repeats[ti];
                    ti_start += repeats[ti];
                }
                int step_pos = fit_step(data_mean) + 1;

                double bat_post = (tis(step_pos) + tis(step_pos + 1)) / 2;
                if (multitau)
                {
                    bat_post -= (taus(step_pos) + taus(step_pos + 1)) / 2;
                }
                else if (seqtau < 10)
                {
                    bat_post -= seqtau;
                }

                posterior.means(bat_index()) = bat_post;
            }
            else if (bat_init == "max")
            {
                // Because this is a const method
                ColumnVector tis_local;
                if (have_tiimg)
                    // Initialize tis from voxelwise image
                    tis_local = tiimg.Column(voxel);
                else
                    tis_local = tis;

                // Federico's method - Time of maximum signal minus the bolus duration (only for
                // differenced data)
                int ind;
                data.Maximum1(ind); // find the point of maximum

                // Find the corresponding TI for this index
                double ti = 0;
                int idx = 1;
                for (int ti_num = 1; ti_num <= tis_local.Nrows(); ti_num++)
                {
                    for (int rpt = 0; rpt < repeats[ti]; rpt++)
                    {
                        if (idx == ind)
                            ti = tis_local(ti_num);
                        idx++;
                    }
                }

                // Calculate the actual TI of the data given slice timing
                double batinitval = ti + slicedt * coord_z;

                // Subtract bolus duration, unless we are assuming it is 'infinite'
                if (multitau)
                {
                    batinitval -= taus(ind);
                }
                else if (seqtau < 10)
                {
                    batinitval -= seqtau;
                }

                // Don't allow the BAT to get too long or too short
                if (batinitval > (timax - 0.5))
                    batinitval = timax - 0.5;
                if (batinitval < 0.01)
                    batinitval = 0.01;
                double tissbatinit = setdelt;
                if (inferart)
                {
                    // arterial component is earliest
                    posterior.means(bat_index() + artidx) = batinitval;
                    if (infertiss)
                    {
                        tissbatinit = batinitval + 0.3;
                        posterior.means(bat_index()) = tissbatinit;
                    }
                }
                else
                {
                    if (infertiss)
                    {
                        tissbatinit = batinitval;
                        posterior.means(bat_index()) = tissbatinit;
                    }
                }
                if (inferwm)
                {
                    posterior.means(bat_index() + wmidx) = tissbatinit + 0.3;
                }
            }
            else if (bat_init != "")
            {
                double bat_post = atof(bat_init.c_str());
                if (bat_post > 0)
                {
                    posterior.means(bat_index()) = bat_post;
                }
                else
                {
                    throw InvalidOptionValue(
                        "bat-init", bat_init, "Must be 'step', 'max' or a number > 0");
                }
            }
        }
    }
}

void ASLFwdModel::GetOutputs(std::vector<std::string> &outputs) const
{
    outputs.push_back("aif");
}

void ASLFwdModel::EvaluateModel(const NEWMAT::ColumnVector &params, 
                                NEWMAT::ColumnVector &result,
                                const std::string &key) const
{
    // ensure that values are reasonable
    // negative check
    ColumnVector paramcpy = params;
    for (int i = 1; i <= NumParams(); i++)
    {
        if (params(i) < 0)
        {
            paramcpy(i) = 0;
        }
    }
    // cout << params.t() << endl;

    // define parameters
    // set defaults to be used if the parameters are fixed (not 'included')
    double ftiss = 0.0;
    double fwm = 0.0;
    double fblood = 0.0;
    double delttiss = setdelt;
    double deltwm = setdeltwm;
    double deltblood = setdeltart;
    double tautiss = seqtau;
    double tauwm = seqtau;
    double taublood = seqtau;
    double T_1 = t1;
    double T_1wm = t1wm;
    double T_1b = t1b;
    double pvgm = 1.0;
    double pvwm = 0.0;
    double taupctiss = 0.0;
    ColumnVector pcvec(1);
    double taupcwm = 0.0;
    ColumnVector pcvecwm(1);
    ColumnVector disptiss;
    ColumnVector dispwm;
    ColumnVector dispart;
    ColumnVector residtiss;
    ColumnVector residwm;
    double g = 1.0;
    double stattiss = 0.0;
    // extract parameter values
    // Flow
    if (inctiss)
    {
        ftiss = paramcpy(flow_index());
    }
    if (incwm)
    {
        fwm = paramcpy(flow_index() + wmidx);
    }
    if (incart)
    {
        fblood = paramcpy(flow_index() + artidx);
    }
    // BAT
    if (incbat)
    {
        if (inctiss)
        {
            delttiss = paramcpy(bat_index());
            if (delttiss > timax - 0.2)
                delttiss = timax - 0.2;
        }
        if (incwm)
        {
            deltwm = paramcpy(bat_index() + wmidx);
            if (deltwm > timax - 0.2)
                deltwm = timax - 0.2;
        }
        if (incart)
        {
            deltblood = paramcpy(bat_index() + artidx);
            if (deltblood > timax - 0.2)
                deltblood = timax - 0.2;
        }
    }
    // Tau
    if (inctau)
    {
        if (septau)
        {
            if (inctiss)
            {
                tautiss = paramcpy(tau_index());
            }
            if (incwm)
            {
                tauwm = paramcpy(tau_index() + wmidx);
            }
            if (incart)
            {
                taublood = paramcpy(tau_index() + artidx);
            }
        }
        else
        {
            tautiss = paramcpy(tau_index());
            tauwm = paramcpy(tau_index());
            taublood = paramcpy(tau_index());
        }
    }
    // T1
    if (inct1)
    {
        if (inctiss)
        {
            T_1 = paramcpy(t1_index());
            if (T_1 < 0.01)
                T_1 = 0.01;
        }
        if (incwm)
        {
            T_1wm = paramcpy(t1_index() + wmidx);
            if (T_1wm < 0.01)
                T_1wm = 0.01;
        }
        T_1b = paramcpy(t1_index() + artidx); // always have T_1b even without
                                              // an explicit arterial component
        if (T_1b < 0.01)
            T_1b = 0.01;
    }
    // PVE
    if (incpve)
    {
        pvgm = paramcpy(pv_index());
        pvwm = paramcpy(pv_index() + 1);
    }
    // taupc
    if (incpc)
    {
        if (inctiss)
        {
            taupctiss = paramcpy(taupc_index());
            pcvec << taupctiss;
        }
        if (incwm)
        {
            taupcwm = paramcpy(taupc_index() + wmidx);
            pcvecwm << taupcwm;
        }
    }
    // dispersion
    if (incdisp)
    {
        if (sepdisp)
        {
            int end = disp_index();
            if (inctiss)
            {
                disptiss = params.Rows(disp_index(), disp_index() - 1 + tiss_model->NumDisp());
                end = disp_index() + tiss_model->NumDisp();
            }
            if (incwm)
            {
                dispwm = params.Rows(end + 1, end + tiss_model->NumDisp());
                end += tiss_model->NumDisp();
            }
            if (incart)
            {
                dispart = params.Rows(end + 1, end + art_model->NumDisp());
            }
        }
        else
        {
            disptiss = params.Rows(disp_index(), disp_index() - 1 + tiss_model->NumDisp());
            dispwm = disptiss;
            dispart = disptiss;
        }
    }
    // exchange (residue function) parameters
    if (incexch)
    {
        // tissue and white matter have same residue function parameters (at present)
        residtiss = params.Rows(resid_index(), resid_index() - 1 + tiss_model->NumResid());
        residwm = residtiss;
    }
    // flip angle correction
    if (incfacorr)
    {
        g = params(facorr_index());
    }
    if (incstattiss)
    {
        // we have a static tissue component
        stattiss = params(stattiss_index());
    }
    // look-locker readout
    // NOTE: this is only valid for even spaced readout
    double T_1ll = T_1b;
    double deltll = 0.0;
    if (looklocker)
    {
        // flip angle correction (if reqd)
        double FAtrue = (g + dg) * FA;
        // modify the T1 values according to the simplified relationship for LL
        // FA in radians
        T_1 = 1 / (1 / T_1 - log(cos(FAtrue)) / dti);
        T_1wm = 1 / (1 / T_1wm - log(cos(FAtrue)) / dti);
        // T1 for blood once it has entered imaging region
        T_1ll = 1 / (1 / T_1b - log(cos(FAtrue)) / dti);
        deltll = deltblood; // the arrival time of the blood within the readout
                            // region i.e. where it sees the LL pulses.
    }
    double f_calib;
    double f_calibwm;
    // if we are using calibrated data then we can use ftiss to calculate T_1app
    // otherwise assume sensible value (units of s^-1)
    if (calib)
    {
        f_calib = ftiss;
        f_calibwm = fwm;
    }
    else
    {
        f_calib = 0.01;
        f_calibwm = 0.003;
    }

    double kctissue;
    kctissue = 0.0;
    double kcblood;
    kcblood = 0.0;
    double kcwm;
    kcwm = 0.0;
    double kcpc = 0.0;
    double kcpcwm = 0.0;
    // Calcualte the Kinetic Model contirbutions at each TI
    // loop over tis
    // We need a local copy of TIs as we need to set this each time
    // when we have a voxelwise TI image and this method is declared const
    ColumnVector tis_local;
    if (have_tiimg)
        tis_local = tiimg.Column(voxel);
    else
        tis_local = tis;

    double ti;
    ColumnVector statcont(tis_local.Nrows()); // this stores the static tissue contribution at each TI
    statcont = 0.0;
    ColumnVector kctotal(tis_local.Nrows()); // this stores the total kinetic signal contribution
    kctotal = 0.0;

    for (int it = 1; it <= tis_local.Nrows(); it++)
    {
        // account here for an increase in the TI due to delays between slices
        int thisz = coord_z; // the slice number
        if (sliceband > 0)
        {
            // multiband setup in which we have specified the number of slices
            // per band
            div_t divresult;
            divresult = div(coord_z, sliceband);
            thisz = divresult.rem; // the number of sluices above the base of
                                   // this band (lowest slice in volume for
                                   // normal (non multi-band) data)
        }
        ti = tis_local(it) + slicedt * thisz; // calcualte the actual TI for this
                                        // slice
        if (multitau)
        {
            // overwrite default tau value with the TI specific one
            tautiss = taus(it);
            tauwm = tautiss;
            taublood = tautiss;
        }
        // look-locker correction for arterial blood
        if (looklocker)
        {
            if (ti > deltll)
                T_1b = T_1ll;
        }

        // crushers
        double artweight = 1.0;
        artweight = 1.0 - crush(it); // arterial weight is opposte of crsuh extent

        // Tissue
        if (inctiss)
            kctissue = pvgm * ftiss
                * tiss_model->kctissue(
                      ti, f_calib, delttiss, tautiss, T_1b, T_1, lambda, casl, disptiss, residtiss);
        // White matter
        if (incwm)
            kcwm = pvwm * fwm
                * tiss_model->kctissue(
                      ti, f_calibwm, deltwm, tauwm, T_1b, T_1wm, lamwm, casl, dispwm, residwm);
        // Arterial
        if (incart)
            kcblood = artweight * fblood
                * art_model->kcblood(ti, deltblood, taublood, T_1b, casl, dispart);

        if (incpc)
        {
            // pre-capilliary component
            // NB its arrival time is the tissue arrival time minus the pc
            // transit time
            if (inctiss)
                kcpc = pvgm * ftiss
                    * pc_model->kctissue(ti, f_calib, delttiss - taupctiss, tautiss, T_1b, T_1,
                          lambda, casl, disptiss, pcvec);
            if (incwm)
                kcpcwm = pvwm * fwm
                    * pc_model->kctissue(ti, f_calibwm, deltwm - taupcwm, tauwm, T_1b, T_1, lamwm,
                          casl, dispwm, pcvecwm);
        }

        // total kinetic contribution
        kctotal(it) = kctissue + kcblood + kcwm + kcpc + kcpcwm;
        // static tissue contribution
        if (incstattiss)
        {
            // a simple static contirbution at all TIs
            statcont(it) = stattiss;
        }

	// T.O. Include CAPRIA signal modulation if needed, which applies equally to all signal components
	if (capria) {
	  kctotal(it) *= capriaR(it);
	}
    }

    // Assemble result
    if (key == "aif") 
    {
        result.ReSize(tis_local.Nrows());
        int thisz = coord_z; // the slice number
        if (sliceband > 0)
        {
            // multiband setup in which we have specified the number of slices
            // per band
            div_t divresult;
            divresult = div(coord_z, sliceband);
            thisz = divresult.rem; // the number of sluices above the base of
                                   // this band (lowest slice in volume for
                                   // normal (non multi-band) data)
        }

        for (int it = 1; it <= tis_local.Nrows(); it++)
        {
            ti = tis_local(it) + slicedt * thisz; // calcualte the actual TI for this
            result(it) = art_model->kcblood(ti, deltblood, taublood, T_1b, casl, dispart);
        }
    }
    else if (key != "") 
    {
        throw InvalidOptionValue("output name", key, "No such output supported");
    }
    else if (hadamard)
    {
        // Hadamard encoded ASL data
        ColumnVector signal;
        signal = HadEncMatrix * kctotal; // collect the KC contirbutions from each block of label
        signal = statcont - signal;      // inlcude static tissue signal
        // loop over repeats - we assume we always get blocks of each set of encodings
        result = signal;
        for (int rpt = 2; rpt <= repeats[0]; rpt++)
        {
            result &= result;
        }
    }
    else if (raw)
    {
        // ASL data - unstracted (raw)
        result.ReSize(2 * tpoints);
        int ent = 0;
        for (int it = 1; it <= tis_local.Nrows(); it++)
        {
            // data is in blocks of repeated TIs
            // loop over the repeats
            for (int rpt = 1; rpt <= repeats[it - 1]; rpt++)
            {
                result(ent + 2 * rpt - 1) = statcont(it) + kctotal(it); // TAG
                result(ent + 2 * rpt) = statcont(it);                   // CONTROL
            }
            ent += 2 * repeats[it - 1];
        }
    }
    else
    {
        // normal (differenced) ASL data
        result.ReSize(tpoints);
        int ent = 0;
        for (int it = 1; it <= tis_local.Nrows(); it++)
        {
            // data is in blocks of repeated TIs
            // loop over the repeats
            for (int rpt = 1; rpt <= repeats[it - 1]; rpt++)
            {
                result(ent + rpt) = statcont(it) + kctotal(it);
            }
            ent += repeats[it - 1];
        }
    }

    if ((key != "aif") && (result.Nrows() != data.Nrows()))
    {
        throw InvalidOptionValue("num volumes", stringify(data.Nrows()),
            "Expected " + stringify(result.Nrows()) + " volumes - check the number of repeats/TIs");
    }
}

FwdModel *ASLFwdModel::NewInstance() { return new ASLFwdModel(); }
void ASLFwdModel::Initialize(ArgsType &args)
{
    //
    // Parameter inference and inclusion
    //

    // Tissue and arterial components - at least one is required
    inctiss = args.ReadBool("inctiss");
    infertiss = args.ReadBool("infertiss");
    inferart = args.ReadBool("inferart");
    incart = args.ReadBool("incart") || inferart;
    if (!inctiss & !incart)
    {
        throw FabberRunDataError("Neither tissue nor arterial components have been selected. "
                                 "At least one of --inctiss and --incart must be set");
    }

    // We only infer WM if we are doing PV correction (below)
    incwm = args.ReadBool("incwm");
    inferwm = false;

    // Bolus arrival time (ATT)
    incbat = args.ReadBool("incbat");
    inferbat = args.ReadBool("inferbat");

    // Method to use when automatically initializing the BAT posterior
    bat_init = args.GetStringDefault("bat-init", "");

    // Pre-capillary component
    inferpc = args.ReadBool("inferpc");
    incpc = args.ReadBool("incpc") || inferpc;

    // Bolus duration
    infertau = args.ReadBool("infertau");
    inctau = args.ReadBool("inctau") || infertau;

    // Separate values of tau for tissue, WM and blood components
    septau = args.ReadBool("septau");

    // Tissue T1 value
    infert1 = args.ReadBool("infert1");
    inct1 = args.ReadBool("inct1") || infert1;

    // incdisp is set based on whether there are dispersion parameters in the model
    inferdisp = args.ReadBool("inferdisp");

    // Separate dispersion parameters for tissue, WM and blood components
    sepdisp = args.ReadBool("sepdisp");

    // incexch is set based on whether there are residue function parameters in the model
    inferexch = args.ReadBool("inferexch");

    // Partial volume estimates - make sure we have a WM component in the model
    incpve = args.ReadBool("incpve");
    incwm = incwm || incpve;

    // Partial volume correction - seems to be just a shorthand for including partial
    // volume estimates and inferring the WM partial volume
    pvcorr = args.ReadBool("pvcorr");
    if (pvcorr)
    {
        incpve = true;
        incwm = true;
        inferwm = true;
    }

    // Static tissue contribution - constant at all times
    incstattiss = args.ReadBool("incstattiss");
    inferstattiss = args.ReadBool("inferstattiss");

    // Parameter indices, relative to flow_index(), for WM and arterial components
    ncomps = (inctiss ? 1 : 0) + (incart ? 1 : 0) + (incwm ? 1 : 0);
    wmidx = 0;
    if (inctiss)
        wmidx++;
    artidx = 0;
    if (inctiss)
        artidx++;
    if (incwm)
        artidx++;

    // ARD selection for aBV. ARD seems to be disabled currently?
    bool ardoff = false;
    ardoff = args.ReadBool("ardoff");
    doard = false;
    if (inferart == true && ardoff == false)
    {
        ardindices.push_back(flow_index() + artidx);
    }

    //
    // Tissue parameters
    //

    // Bolus arrival time (ATT)
    setdelt = args.GetDoubleDefault("bat", 0.7);

    // By default BAT in WM is longer then GM
    setdeltwm = args.GetDoubleDefault("batwm", setdelt + 0.3);

    // By default BAT in blood is shorter then GM
    setdeltart = args.GetDoubleDefault("batart", setdelt - 0.3);

    // BAT std dev (same for all tissue)
    double deltsd = args.GetDoubleDefault("batsd", 0.316);
    deltprec = 1 / (deltsd * deltsd);

    // Arterial BAT std dev - by default the arterial BAT SD is same as tissue
    deltsd = args.GetDoubleDefault("batartsd", deltsd);
    deltartprec = 1 / (deltsd * deltsd);

    // Possible different default for T1 if we have a WM component (since then the 'tissue' is
    // GM only rather than mixed WM/GM)
    t1 = args.GetDoubleDefault("t1", incwm ? 1.3 : 1.3);
    t1b = args.GetDoubleDefault("t1b", 1.65);
    t1wm = args.GetDoubleDefault("t1wm", 1.1);

    // Different default for lambda if we have a WM component (since then the 'tissue' is
    // GM only rather than mixed WM/GM)
    lambda = args.GetDoubleDefault("lambda", incwm ? 0.98 : 0.9);
    lamwm = args.GetDoubleDefault("lambdawm", 0.82);

    // T.O. Add dispersion parameter prior SD
    double dispsd = args.GetDoubleDefault("dispsd", -1.0);
    if (dispsd == -1.0) // No (sensible) value set by user, so use defaults later
    {
        dispprecset = -1.0;
    } else {
        dispprecset = 1 / (dispsd * dispsd);
    }

    //
    // Data acquisition parameters
    //

    // Whether data has been tag-control subtracted or not
    string iaf = args.GetStringDefault("iaf", "diff");
    raw = false;
    if ((iaf == "tc") || (iaf == "ct"))
    {
        // data is in raw (non-subtracted) form
        raw = true;
    }
    tagfirst = (iaf != "ct");

    // Data is CASL or PASL (default)
    casl = args.GetBool("casl");

    // data has already been subjected to calibration
    calib = args.GetBool("calib");

    // Increase in TI per slice
    slicedt = args.GetDoubleDefault("slicedt", 0);

    // Number of slices in a band in a multi-band setup (zero implies single band)
    sliceband = args.GetIntDefault("sliceband", 0);

    // Saturation of the bolus a fixed time pre TI measurement
    pretisat = args.GetDoubleDefault("pretisat", 0);

    // Timing parameters (TIs / PLDs)
    have_tiimg = false;
    vector<double> ti_list = args.GetDoubleList("ti");
    bool ti_set = (ti_list.size() >= 1);
    unsigned int num_times = ti_list.size();

    vector<double> pld_list = args.GetDoubleList("pld");
    bool pld_set = (pld_list.size() >= 1);
    if (pld_set)
        num_times = pld_list.size();

    if (ti_set && pld_set)
        throw FabberRunDataError("Cannot specify TIs and PLDs at the same time");
    else if (num_times == 0)
    {
        if (!args.HaveKey("tiimg"))
        {
            throw FabberRunDataError("No TIs/PLDs specified");
        }
        else 
        {
            // A TI image is a voxelwise list of TIs. Note that this may be one
            // per data volume (requiring repeats=1) or it may interact with
            // the repeats in the normal way.
            tiimg  = args.GetVoxelData("tiimg");
            have_tiimg = true;
            num_times = tiimg.Nrows();
            if (slicedt != 0)
            {
                throw InvalidOptionValue("slicedt", stringify(slicedt),
                                         "Can't specify slicedt when providing voxelwise TI image");
            }
        }
    }

    // Repeats - may be single repeat or multiple (one per TI/PLD)
    if (args.HaveKey("rpt1"))
    {
        // Multiple repeats have been specified
        repeats = args.GetIntList("rpt");
        if (repeats.size() != num_times)
        {
            throw InvalidOptionValue("Number of repeats", stringify(repeats.size()),
                "Mismatch between number of inflow times (TIs/PLDs) and repeats - should be equal");
        }
    }
    else
    {
        // Multiple repeats not specified - use single value defaulting to 1
        for (unsigned int it = 0; it < num_times; it++)
        {
            repeats.push_back(args.GetIntDefault("repeats", 1));
        }
    }

    // Bolus durations - may be single value or multiple (one per TI/PLD)
    // Default of 1000 is effectively infinite
    seqtau = args.GetDoubleDefault("tau", 1000);
    vector<double> taus_list = args.GetDoubleList("tau");
    multitau = taus_list.size() > 1;
    if (multitau)
    {
        if (inctau)
        {
            throw InvalidOptionValue("Number of taus", stringify(taus_list.size()),
                "Inference/Inclusion of (variable) bolus duration "
                "is not compatible with multiple bolus durations "
                "use a single value (--tau=)");
        }
        if (septau)
        {
            throw InvalidOptionValue("Number of taus", stringify(taus_list.size()),
                "Separate bolus duration for different components "
                "is not compatible with multiple bolus durations "
                "use a single value (--tau=)");
        }
        if (taus_list.size() != num_times)
        {
            throw InvalidOptionValue("Number of taus", stringify(taus_list.size()),
                "Mismatch between number of inflow times (TIs/PLDs) and "
                "bolus durations - these should be equal");
        }
        taus.ReSize(taus_list.size());
        for (unsigned int i = 0; i < taus_list.size(); i++)
            taus(i + 1) = taus_list[i];
    }

    // Time encoded data
    // if nothing is stated, no time encoding is assumed.

    // If hadamard option is set to an integer N, N Hadamard encoding lines are assumed.
    hadamard = args.HaveKey("hadamard");

    // In some cases the full Hadamard matrix is needed,
    // i.e. all N rows (and not only N-1). This is activated by this command.
    // FIXME this looks wrong, why the OR?
    bool FullHad = hadamard || args.ReadBool("fullhad");

    // Indicates if a Walsh-sorted Hadamard matrix was used for encoding.
    bool Walsh = args.ReadBool("walsh");

    // Indicates if a non-hadamard matrix was used for encoding.
    bool NonHadamardMatrix = args.ReadBool("NonHadamardMatrix");

    // Initialize the TIs list
    if (hadamard)
    {
        if (repeats.size() > 1)
        {
            throw InvalidOptionValue("Number of repeats", stringify(repeats.size()),
                "Cannot specify more than one set of repeats with Hadmard encoding");
        }
        if (have_tiimg)
        {
            throw InvalidOptionValue("TI image", stringify(num_times),
                "Hadamard time encoding cannot be used with a voxelwise TI image");
        }
        if (num_times != 1)
        {
            throw InvalidOptionValue("Number of TIs/PLDs", stringify(num_times),
                "Hadamard time encoding requires exactly one TI/PLD");
        }

        // List of TIs to skip in Hadamard encoding
        // FIXME Could this be implemented using masked timepoints?
        vector<int> skip_list = args.GetIntList("skip");

        // Size must be power of 2 or 12 (special case). Note clever
        // bitwise method for identifying powers of 2!
        HadamardSize = args.GetInt("hadamard");
        if (((HadamardSize & (HadamardSize - 1)) != 0) && (HadamardSize != 12))
        {
            throw InvalidOptionValue("hadamard", stringify(HadamardSize),
                "Hadamard encoding is only possible "
                "with a number of encodings that are "
                "modulo 2 (2,4,8,16...) or number of "
                "encodings equal to 12");
        }

        NumberOfSubBoli = HadamardSize;
        if (!FullHad)
            NumberOfSubBoli--;

        HadEncMatrix = HadamardMatrix(HadamardSize, FullHad, Walsh, NonHadamardMatrix, skip_list);

        // we will usually ingore the first column (all control), unless we are doing FullHad
        int start_col = HadamardSize - (NumberOfSubBoli - 1);
        HadEncMatrix = HadEncMatrix.SubMatrix(1, HadamardSize, start_col, HadamardSize);

        // FIXME skip list?
        HadEncMatrix
            = HadEncMatrix.SubMatrix(1, HadamardSize - skip_list.size(), 1, NumberOfSubBoli);

        // calculate the TIs corresponding to the indiviual subboli
        tis.ReSize(NumberOfSubBoli);
        for (int i = 0; i < NumberOfSubBoli; i++)
        {
            ColumnVector tmp(1);
            if (multitau)
                tis(i + 1) = pld_list[0] + (NumberOfSubBoli - i) * taus_list[i];
            else
                tis(i + 1) = pld_list[0] + (NumberOfSubBoli - i) * seqtau;
        }
    }
    else
    {
        // normal ASL data. Note that if we have a TI image this
        // code will leave the tis vector uninitialized - this is
        // fine we will do it voxelwise in the Evaluate method
        tis.ReSize(num_times);
        if (ti_set)
        {
            for (unsigned int i = 0; i < num_times; i++)
                tis(i + 1) = ti_list[i];
        }
        if (pld_set)
        {
            if (casl)
            {
                for (unsigned int i = 0; i < num_times; i++)
                {
                    if (multitau)
                        tis(i + 1) = pld_list[i] + taus_list[i];
                    else
                        tis(i + 1) = pld_list[i] + seqtau;
                }
            }
            else
            {
                // unlikely to happen, but permits the user to supply PLDs for a pASL acquisition
                for (unsigned int i = 0; i < num_times; i++)
                    tis(i + 1) = pld_list[i];
            }
        }
    }

    // Initialise BAT only if number of TIs > 1
    if (num_times == 1)
    {
        bat_init = "";
    }

    // Total number of time points in data
    tpoints = 0;
    for (int it = 0; it < num_times; it++)
    {
        if (repeats.size() > 1)
        {
            tpoints += repeats[it];
        }
        else
        {
            tpoints += repeats[0];
        }
    }

    // The final TI
    if (have_tiimg)
        timax = tiimg.Maximum();
    else
        timax = tis.Maximum();

    // Vascular crushing
    // FIXME use GetStringList

    // To specify a custom combination of vascular crushing
    string crush_temp = args.ReadWithDefault("crush1", "notsupplied");

    // If crush_temp = none then we assume all data has same crushing
    // parameters we will represent this as no crushers
    crush.ReSize(num_times);
    crush = 0.0; // default is no crusher
    crushdir.ReSize(num_times, 3);
    crushdir = 0.0;
    crushdir.Column(3) = 1.0; // default (which should remain ignored normally) is z-only
    if (crush_temp != "notsupplied")
    {
        // we need to assemble crusher information
        int N = 1;
        while (true)
        {
            if (N > 1)
            {
                crush_temp = args.ReadWithDefault("crush" + stringify(N), "stop!");
                if (crush_temp == "stop!")
                    break; // we have run out of crusher specifications
            }
            // determine what crusher type we have
            if (crush_temp == "off" || crush_temp == "none")
            {
                crush(N) = 0.0;
            }
            else if (crush_temp == "on")
            {
                crush(N) = 1.0;
            }
            else if (crush_temp == "xyz")
            {
                crush(N) = 1.0;
                crushdir.Row(N) << 1 << 1 << 1;
                crushdir.Row(N) /= sqrt(3);
            }
            else if (crush_temp == "-xyz")
            {
                crush(N) = 1.0;
                crushdir.Row(N) << -1 << 1 << 1;
                crushdir.Row(N) /= sqrt(3);
            }
            else if (crush_temp == "x-yz")
            {
                crush(N) = 1.0;
                crushdir.Row(N) << 1 << -1 << 1;
                crushdir.Row(N) /= sqrt(3);
            }
            else if (crush_temp == "-x-yz")
            {
                crush(N) = 1.0;
                crushdir.Row(N) << -1 << -1 << 1;
                crushdir.Row(N) /= sqrt(3);
            }
            N++;
        }
    }

    // Look-Locker correction
    string FAin = args.ReadWithDefault("FA", "none");
    if (FAin == "none")
        looklocker = false;
    else
    {
        looklocker = true;
        FA = convertTo<double>(FAin);
        // FA required in radians
        FA *= M_PI / 180;
        // NOTE LL correction is only valid with evenly spaced TIs
        if (have_tiimg)
        {
            throw InvalidOptionValue("FA", FAin, "Look-locker correction cannot be used with voxelwise TI image");
        }
        if (num_times < 2)
        {
            throw InvalidOptionValue("FA", FAin, "Look-locker correction can only be used with >1 evenly spaced TIs");
        }
        dti = tis(2) - tis(1);
    }
    // Indicate that we want to do FA correction - the g image will
    // need to be separately loaded in as an image prior
    incfacorr = args.ReadBool("facorr");
    dg = convertTo<double>(args.ReadWithDefault("dg", "0.0"));

    // ---------------------------
    // T.O. CAPRIA Initialisation
    // ---------------------------
    capria = args.ReadBool("capria"); // Use CAPRIA options?
    if (capria) {
      // Read in arguments
      capriafa1 = convertTo<double>(args.ReadWithDefault("capriafa1", "-1.0")); // CAPRIA starting flip angle
      capriafa2 = convertTo<double>(args.ReadWithDefault("capriafa2", "-1.0")); // CAPRIA finishing flip angle
      capriatr = convertTo<double>(args.ReadWithDefault("capriatr", "-1.0"));  // CAPRIA repetition time

      // Make sure parameters are defined
      if (capriafa1 < 0.0)
      {
	throw InvalidOptionValue("capriafa1",stringify(capriafa1),
				 "CAPRIA FA1 must be set and greater than zero");
      }
      if (capriafa2 < 0.0)
      {
	throw InvalidOptionValue("capriafa2",stringify(capriafa2),
				 "CAPRIA FA2 must be set and greater than zero");
      }
      if (capriatr < 0.0)
      {
	throw InvalidOptionValue("capriatr",stringify(capriatr),
				 "CAPRIA TR must be set and greater than zero");
      }

      // Calculate some useful parameters
      dti = tis(2) - tis(1); // Assume constant spacing
      Nsegs = std::round( dti / capriatr ); // Number of segments (excitations per frame per ASL prep)
      int Ntrs = Nsegs * tis.Nrows(); // Number of TR periods after each ASL prep

      // Calculate the flip angle schedule
      ColumnVector capriaFAs;
      capriaFAs.ReSize(Ntrs);
      //      LOG << "CAPRIA flip angles:" << endl;
      for (int ii=1; ii<=Ntrs; ii++) {
	capriaFAs(ii) = capriafa1 + (capriafa2-capriafa1)*pow((double)(ii-1)/(double)(Ntrs-1),2);
	//	LOG << capriaFAs(ii) << endl;
      }

      // Calculate the CAPRIA attenuation factors at each TR
      ColumnVector Ract;
      Ract.ReSize(Ntrs);
      Ract(1) = 1.0; // No attenuation prior to the first pulse
      //      LOG << "Ract:" << endl;
      //      LOG << Ract(1) << endl;
      for (int ii=2; ii<=Ntrs; ii++) {
	// Each pulse leaves cos(alpha) of the previous Mz available
	Ract(ii) = Ract(ii-1) * cos(capriaFAs(ii-1)*M_PI/180.0);
	//	LOG << Ract(ii) << endl;
	// NB. Agrees with Matlab implementation
      }

      // Incoporate the sin(alpha) term which scales the signal strength
      //      LOG << "CAPRIA signal scaling:" << endl;
      for (int ii=1; ii<=Ntrs; ii++) {
	// The signal is scaled by sin(alpha)
	Ract(ii) = Ract(ii) * sin(capriaFAs(ii)*M_PI/180.0);
	// LOG << Ract(ii) << endl;
	// NB. Agrees with Matlab implementation
      }

      // Downsample to the reconstructed temporal resolution
      capriaR.ReSize(tis.Nrows());
      //      LOG << "CAPRIA signal scaling at downsampled resolution:" << endl;
      for (int jj=1; jj<=tis.Nrows(); jj++) {
	int StartIdx = (jj-1)*Nsegs + 1;
	int EndIdx = jj*Nsegs;

	// Sum up all the actual signal strength scaling factors in this frame
	capriaR(jj) = 0;
	for (int ii=StartIdx; ii<=EndIdx; ii++){
	  capriaR(jj) += Ract(ii);
	}

	// Divide by the number of elements to get the mean signal scaling
	capriaR(jj) = capriaR(jj) / ((double) Nsegs );
	// NB. Multiplying this factor by the Buxton model evaluated at the midpoint 
	// of each frame is not exactly equivalent to calculating the Buxton model at
	// each TR, scaling it by the CAPRIA factors and then averaging, but it's a 
	// pretty close approximation

	//	LOG << capriaR(jj) << endl;
	// NB. Agrees with Matlab implementation
      }
      

    } // end capria

    // Kinetic model definition
    //
    // Up to three kinetic models are used - one for the arterial component,
    // one for the pre-capillary component and one for the tissue component.
    // Each kinetic model is a convolution of an input function (AIF) and a
    // residue function.
    //
    // The 'disp' option controls the residue function for the arterial
    // component - no dispersion is the default. The input function for
    // the arterial component is the ASL bolus - this is accounted for in
    // the AIF model.
    //
    // The pre-capillary model uses the arterial component as an input and
    // impermeable residue model. This is used to model a separate pre-capillary
    // transit time. Note that the ATT (tau) remains the tissue transit time.
    //
    // The input function for the tissue component is the arterial component.
    // The 'exch' option controls the residue function. The default is 'mix'
    // for a well-mixed model.
    //
    // Without dispersion, analytic solutions are available for the tissue
    // and pre-capillary models. An analytic solution for the tissue model
    // is also available for disp=gamma and exch=mix. If an analytic solution
    // is not available a numerical model is used which is provided directly
    // with the input and residue functions.
    //
    // The same tissue model is uesd for GM and WM
    //
    // NB it is feasible to have different dispersion for arterial and
    // tissue models, but this is not implemented

    // Set the AIF dispersion type
    string disptype = args.ReadWithDefault("disp", "none");

    // Set the type of exchange in tissue compartment
    string exchtype = args.ReadWithDefault("exch", "mix");

    // Force numerical convolution for the evaluation of the model
    bool forceconv = args.ReadBool("forceconv");
    double conv_delta = args.GetDoubleDefault("conv-delta", 0.1);

    // Arterial Model (only depends upon dispersion type)
    art_model = NULL;
    if (disptype == "none")
    {
        art_model = new AIFModel_nodisp(args.GetDoubleDefault("aif-leadscale", 0.01));
    }
    else if (disptype == "gamma")
    {
        art_model = new AIFModel_gammadisp();
    }
    else if (disptype == "gauss")
    {
        art_model = new AIFModel_gaussdisp();
    }
    else if (disptype == "sgauss")
    {
        art_model = new AIFModel_spatialgaussdisp();
    }
    else if (disptype == "sgauss_alt")
    {
        art_model = new AIFModel_spatialgaussdisp_alternate();
    }
    else
    {
        throw InvalidOptionValue("disptype", disptype, "Dispersion type not recognized");
    }

    // PC model
    pc_model = NULL;
    if (disptype == "none")
    {
        pc_model = new TissueModel_nodisp_imperm();
    }
    // Default is to use convolution model with the impermeable residue function
    if ((pc_model == NULL) || forceconv)
    {
        ResidModel *imperm_resid = new ResidModel_imperm();
        pc_model = new TissueModel_aif_residue(art_model, imperm_resid, conv_delta);
    }

    // Tissue model
    //
    // This should ALWAYS set a residual model and maybe a MATCHING tissue model
    tiss_model = NULL;
    resid_model = NULL;

    if (exchtype == "mix")
    {
        // Well mixed single compartment
        resid_model = new ResidModel_wellmix();
        if (disptype == "none")
        {
            tiss_model = new TissueModel_nodisp_wellmix();
        }
        else if ((disptype == "gamma")) // T.O. Allow gamma dispersion with the (P)CASL models now
        {
            if(casl) // T.O. Debugging
            {
                LOG << "Using gamma dispersion model with (P)CASL." << endl; 
            }
            tiss_model = new TissueModel_gammadisp_wellmix();
        }
    }
    else if (exchtype == "simple")
    {
        // Simple single impermeable compartment that decays with T1b
        resid_model = new ResidModel_simple();
        if (disptype == "none")
        {
            tiss_model = new TissueModel_nodisp_simple();
        }
    }
    else if (exchtype == "2cpt")
    {
        // 2 compartment exchange (simplest model - no backflow, no venous outflow)
        resid_model = new ResidModel_twocpt();
        if (disptype == "none")
        {
            string solution = args.GetStringDefault("2cpt-solution", "slow");
            double mtt_prior = args.GetDoubleDefault("mtt", 1.0, 0);
            tiss_model = new TissueModel_nodisp_2cpt(solution, mtt_prior);
        }
    }
    else if (exchtype == "spa")
    {
        // SPA 2 compartment model
        resid_model = new ResidModel_spa();
        if ((disptype == "none") && !casl)
        {
            tiss_model = new TissueModel_nodisp_spa();
        }
    }

    // Default is to use the convolution model with the residue function
    // if no analytical tissue model can be found
    if ((tiss_model == NULL) || args.ReadBool("forceconv"))
    {
        // note the 'forceconv' option, this forces the model to use the
        // convolution formulation voer any analytic form it has
        if (resid_model == NULL)
        {
            // we cannot do a convolution in this case as a residue function
            // has not been found either!
            throw InvalidOptionValue("exch", exchtype,
                "A residue function model for this exchange type cannot be found");
        }
        else
        {
            tiss_model = new TissueModel_aif_residue(art_model, resid_model, conv_delta);
        }
    }

    // Include dispersion parameters if the model has them
    if ((art_model->NumDisp() > 0) | (tiss_model->NumDisp() > 0))
    {
        incdisp = true;
    }

    // Include resdue function (i.e. exchange) parameters is the model has them
    if (tiss_model->NumResid() > 0)
    {
        incexch = true;
    }

    // Dispersion model - load priors from command line
    // NB this could now be replaced with PSP_byname_mean etc
    for (int i = 1; i <= art_model->NumDisp(); i++)
    {
        string priormean = args.ReadWithDefault("disp_prior_mean_" + stringify(i), "null");
        if (priormean != "null")
        {
            // ASSUME that dispersion
            // model same for arterial
            // and tissue model (and they
            // share the same priors)
            art_model->SetPriorMean(i, convertTo<double>(priormean));
            tiss_model->SetDispPriorMean(i, convertTo<double>(priormean));
        }
    }

    // Set default voxel coordinates to a default of 0 (for the times we
    // call the model before we start handling data)
    coord_x = 0;
    coord_y = 0;
    coord_z = 0;

    // Write information about the parameters to the log
    LOG << "Inference using resting state ASL model" << endl;
    LOG << "INCLUSIONS:" << endl;
    if (inctiss)
        LOG << "Tissue" << endl;
    if (incart)
        LOG << "Arterial/MV" << endl;
    if (incwm)
        LOG << "White matter" << endl;
    LOG << "Number of components: " << ncomps << endl;
    if (incbat)
        LOG << "Bolus arrival time" << endl;
    if (incpc)
        LOG << "Pre capilliary" << endl;
    if (inctau)
        LOG << "Bolus duration (tau)" << endl;
    if (inct1)
        LOG << "T1 values" << endl;
    if (incdisp)
        LOG << "Dispersion" << endl;
    if (incexch)
        LOG << "Restricted exchange" << endl;
    if (incpve)
        LOG << "Partial volume estimates" << endl;
    if (incstattiss)
        LOG << "Statis tissue" << endl;
    LOG << "INFERENCE:" << endl;
    if (infertiss)
        LOG << "Tissue" << endl;
    if (inferart)
        LOG << "Arterial/MV" << endl;
    if (inferwm)
        LOG << "White matter" << endl;
    if (inferbat)
        LOG << "Bolus arrival time" << endl;
    if (inferpc)
        LOG << "Pre capilliary" << endl;
    if (infertau)
        LOG << "Bolus duration (tau)" << endl;
    if (septau)
        LOG << "  Separate values of tau for each component" << endl;
    if (infert1)
        LOG << "T1 values" << endl;
    if (inferdisp)
        LOG << "Dispersion" << endl;
    if (sepdisp)
        LOG << "  Separate dispersion parameters for each component" << endl;
    if (inferexch)
        LOG << "Restricted exchange" << endl;
    if (pvcorr)
        LOG << "Partial volume correction" << endl;
    if (inferstattiss)
        LOG << "Static tissue" << endl;
    LOG << "------" << endl;
    // kinetic model
    LOG << "Kinetic model:" << endl;
    if (!casl)
        LOG << "Data being analysed using PASL inversion profile" << endl;
    if (casl)
        LOG << "Data being analysed using CASL inversion profile" << endl;
    LOG << "Tissue model:" << tiss_model->Name() << endl;
    if (tiss_model->NumDisp() > 0)
        LOG << "Dispersion parameter priors (means then precisions): " << tiss_model->DispPriors()
            << endl;
    if (tiss_model->NumResid() > 0)
        LOG << "Residue function parameter priors (means then precisions): "
            << tiss_model->ResidPriors() << endl;
    if (incart)
        LOG << "Arterial model:" << art_model->Name() << endl;
    if (art_model->NumDisp() > 0)
        LOG << "Dispersion parameter priors (means then precisions): " << art_model->Priors()
            << endl;
    if (incpc)
        LOG << "Pre-capillary model:" << pc_model->Name() << endl;
    LOG << "------" << endl;
    // scan parameters
    if (pretisat > 0)
        LOG << "Saturation of " << pretisat << " s before TI has been specified" << endl;
    if (calib)
        LOG << "Input data is in physioligcal units, using estimated CBF "
               "in T_1app calculation"
            << endl;
    LOG << "Data parameters: #repeats = " << repeats << endl;
    LOG << " t1 = " << t1 << ", t1b = " << t1b;
    if (incwm)
        LOG << ", t1wm= " << t1wm << endl;
    LOG << " bolus duration (tau) = " << seqtau << endl;
    // Hadamard
    if (hadamard)
    {
        LOG << "Hadamard time encoding:" << endl;
        if (FullHad)
            LOG << "  Full nxn-Hadamard matrix is used." << endl;
        if (Walsh)
            LOG << "  Walsh-sorted Hadamard matrix is used." << endl;
        if (NonHadamardMatrix)
            LOG << "  Non-Hadamard matrix is used." << endl;
        LOG << "  Number of Hadamard encoded images: " << HadamardSize << endl;
        LOG << "  Number of Hadamard subboli: " << NumberOfSubBoli << endl;
        LOG << "  Subbolus length: " << seqtau << endl;
        LOG << "  Post labeling delay (PLD): " << pld_list[1] << endl;
    }

    // T.O. CAPRIA
    if (capria)
    {
        LOG << "CAPRIA is on:" << endl;
        LOG << "  Starting flip  angle/degs: " << capriafa1 << endl;
        LOG << "  Finishing flip angle/degs: " << capriafa2 << endl;
        LOG << "  ASSUME QUADRATIC FLIP ANGLE VARIATION! "  << endl;
        LOG << "  Excitation pulse TR/s: " << capriatr << endl;
        LOG << "  Calculated Nsegs: " << Nsegs << endl;
    }

    if (doard)
    {
        LOG << "ARD has been set on arterial compartment " << endl;
    }
    if (have_tiimg)
    {
        LOG << "Using voxelwise TI image" << endl;
    }
    else
    {
        LOG << "TIs: ";
        for (int i = 1; i <= tis.Nrows(); i++)
            LOG << tis(i) << " ";
        LOG << endl;
    }
}

void ASLFwdModel::NameParams(vector<string> &names) const
{
    names.clear();

    if (inctiss)
        names.push_back("ftiss");
    if (incwm)
        names.push_back("fwm");
    if (incart)
        names.push_back("fblood");
    if (incbat)
    {
        if (inctiss)
            names.push_back("delttiss");
        if (incwm)
            names.push_back("deltwm");
        if (incart)
            names.push_back("deltblood");
    }
    if (inctau)
    {
        if (septau)
        {
            if (inctiss)
                names.push_back("tautiss");
            if (incwm)
                names.push_back("tauwm");
            if (incart)
                names.push_back("taublood");
        }
        else
        {
            names.push_back("tau");
        }
    }
    if (inct1)
    {
        if (inctiss)
            names.push_back("T_1");
        if (incwm)
            names.push_back("T_1wm");
        names.push_back("T_1b");
    }
    if (incpve)
    {
        names.push_back("pvgm");
        names.push_back("pvwm");
    }
    if (incpc)
    {
        if (inctiss)
            names.push_back("taupctiss");
        if (incwm)
            names.push_back("taupcwm");
    }
    if (incdisp)
    {
        if (sepdisp)
        {
            if (inctiss)
            {
                for (int i = 1; i <= tiss_model->NumDisp(); i++)
                {
                    names.push_back("disptiss" + stringify(i));
                }
            }
            if (incwm)
            {
                for (int i = 1; i <= tiss_model->NumDisp(); i++)
                {
                    names.push_back("dispwm" + stringify(i));
                }
            }
            if (incart)
            {
                for (int i = 1; i <= art_model->NumDisp(); i++)
                {
                    names.push_back("dispart" + stringify(i));
                }
            }
        }
        else
        {
            for (int i = 1; i <= tiss_model->NumDisp(); i++)
            {
                names.push_back("disp" + stringify(i));
            }
        }
    }
    if (incexch)
    {
        for (int i = 1; i <= tiss_model->NumResid(); i++)
        {
            names.push_back("exch" + stringify(i));
        }
    }
    if (incstattiss)
    {
        names.push_back("stattiss");
    }
}

// Short functions used for Walsh ordering of the Hadamard matrix

static long BitReverse(long input, int maxInt)
{
    long result = 0;
    for (int i = 0; i < maxInt; ++i)
    {
        if (1 << (maxInt - 1 - i) & input)
            result |= 1 << i;
    }
    return result;
}

static long GrayCode(long input)
{
    long result = 0;
    result = input ^ (input >> 1);
    return result;
}

static long WalshSorting(long input, int maxInt)
{
    long result = 0;
    result = BitReverse(GrayCode(input), maxInt);
    return result;
}

static long ComplementaryMatrix(long input)
{
    long result = 0;
    result = long(input / 2);
    return result;
}

// Useful other functions
Matrix ASLFwdModel::HadamardMatrix(const int size, const bool FullHad, const bool Walsh,
    const bool NonHadamardMatrix, const vector<int> &skip_list) const
{
    // generate a Hadamard matrix
    // This has zeros (for control) and ones (for label) in place of the classic +1 and -1

    Matrix matrix, matrix_temp, norm_matrix, comp_matrix;

    if (size == 12)
    {
        // the 12x12 Hadamard matrix cannot be generated by  Sylvester's costruction and is
        // therefore
        // hardcoded here as a special case
        Real b[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1,
            1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1,
            0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0,
            0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1,
            1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0,
            1, 1 };
        matrix.ReSize(12, 12);
        matrix << b;
    }
    else if (NonHadamardMatrix)
    {
        Real b[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
            1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
            1, 0, 1, 0, 1, 0, 1, 0, 1 };
        matrix.ReSize(size, size);
        matrix << b;
    }
    else
    {
        // Build the Hadamard matrix - using Sylvester's construction
        Matrix H2k(2, 2);
        H2k << 1.0 << 1.0 << 1.0 << -1.0;
        matrix_temp = H2k;
        // TODO check that size is a power of 2 first
        for (int i = 4; i <= size; i *= 2)
        {
            matrix_temp = KP(matrix_temp, H2k);
        }

        matrix.ReSize(size, size);

        if (Walsh)
        {
            // Walsh-sorting of the Hadamard matrix
            int Nbits = log2(matrix_temp.Nrows());
            for (long l = 0; l < matrix_temp.Nrows(); l++)
            {
                matrix.Row(l + 1) << matrix_temp.Row(WalshSorting(l, Nbits) + 1);
            }
        }
        else
        {
            matrix << matrix_temp;
        }

        norm_matrix.ReSize(size, size);
        norm_matrix = abs(0.5 * (matrix - 1.0));

        if (FullHad)
        {
            // For the case, when a mirrored, interleaved, 2NxN matrix was used for encoding
            // (FullHad==true), the normal matrix norm_matrix and the mirrored matrix mirr_matrix
            // have to
            // be taken into account
            matrix.ReSize(2 * size, size);
            comp_matrix.ReSize(size, size);
            comp_matrix = abs(0.5 * (-matrix - 1.0));
            for (int k = 1; k <= matrix.Nrows(); k++)
            {
                if (k % 2 == 1)
                {
                    matrix.Row(k) << norm_matrix.Row(ComplementaryMatrix(k + 1));
                }
                else
                {
                    matrix.Row(k) << comp_matrix.Row(ComplementaryMatrix(k + 1));
                }
            }
        }
        else
        {
            // Otherwise the mirrored matrix is not used
            if (skip_list.size() > 0)
            {
                std::cout << "Skipped images: " << std::endl;
                std::cout << skip_list << std::endl;
                matrix.ReSize(size - skip_list.size(), size);

                unsigned int i = 0;
                unsigned int k = 0;
                while (i < size - skip_list.size())
                {
                    if (k < skip_list.size())
                    {
                        if (int(i + k + 1) == skip_list[k])
                        {
                            k++;
                        }
                        else
                        {
                            matrix.Row(i + 1) << norm_matrix.Row(i + k + 1);
                            std::cout << i + 1 << "  " << i + k + 1 << std::endl;
                            i++;
                        }
                    }
                    else
                    {
                        matrix.Row(i + 1) << norm_matrix.Row(i + k + 1);
                        std::cout << i + 1 << "  " << i + k - 1 << std::endl;
                        i++;
                    }
                }
            }
            else
            {
                matrix << norm_matrix;
            }
        }
    }

    return matrix;
}

void ASLFwdModel::SetupARD(const MVNDist &theta, MVNDist &thetaPrior, double &Fard) const
{
    int ardindex = ard_index();
    if (doard)
    {
        SymmetricMatrix PriorPrec;
        PriorPrec = thetaPrior.GetPrecisions();

        PriorPrec(ardindex, ardindex) = 1e-12; // set prior to be initally non-informative

        thetaPrior.SetPrecisions(PriorPrec);
        thetaPrior.means(ardindex) = 0;
        // set the Free energy contribution from ARD term
        SymmetricMatrix PostCov = theta.GetCovariance();
        double b
            = 2 / (theta.means(ardindex) * theta.means(ardindex) + PostCov(ardindex, ardindex));
        Fard = -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5)
            - 0.5 * log(b); // taking c as 0.5 - which it will be!
    }
}

void ASLFwdModel::UpdateARD(const MVNDist &theta, MVNDist &thetaPrior, double &Fard) const
{
    int ardindex = ard_index();
    if (doard)
    {
        SymmetricMatrix PriorCov;
        SymmetricMatrix PostCov;
        PriorCov = thetaPrior.GetCovariance();
        PostCov = theta.GetCovariance();
        PriorCov(ardindex, ardindex)
            = theta.means(ardindex) * theta.means(ardindex) + PostCov(ardindex, ardindex);

        thetaPrior.SetCovariance(PriorCov);
        // Calculate the extra terms for the free energy
        double b
            = 2 / (theta.means(ardindex) * theta.means(ardindex) + PostCov(ardindex, ardindex));
        Fard = -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5)
            - 0.5 * log(b); // taking c as 0.5 - which it will be!
    }
}
