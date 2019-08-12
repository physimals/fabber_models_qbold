/**
 * fwdmodel_qbold.cc - ASE qBOLD curve fitting model
 *
 * Matthew Cherukara, IBME
 *
 * Copyright (C) 2019 University of Oxford  
 */

#include "fwdmodel_qbold.h"

#include <fabber_core/fwdmodel.h>
#include <fabber_core/priors.h>

#include <newmat.h>

#include <vector>
#include <string>
#include <complex>

using namespace std;
using namespace NEWMAT;

// Complex number i
const complex<double> I(0, 1);

// Nuclear gyromagnetic ratio of H (10^6 rad s^-1 T^-1)
const double GAMMA = 267.522;

// Susceptibility difference between deoxy and oxy blood
const double D_CHI0 = 0.264;

FactoryRegistration<FwdModelFactory, R2primeFwdModel>
    R2primeFwdModel::registration("qboldR2p");

FwdModel *R2primeFwdModel::NewInstance() // unchanged
{
    return new R2primeFwdModel();
}

string R2primeFwdModel::GetDescription() const 
{
    return "ASE qBOLD model R2-prime version";
}

static OptionSpec OPTIONS[] = {
    { "tau<n>", OPT_FLOAT, "Tau values (s)", OPT_REQ, "" },
    { "te", OPT_FLOAT, "Single TE value", OPT_NONREQ, "" },
    { "te<n>", OPT_FLOAT, "Sequence of TE values, alternative to --te", OPT_NONREQ, "" },
    { "tr", OPT_FLOAT, "TR value", OPT_REQ, "3.0" },
    { "ti", OPT_FLOAT, "TI value", OPT_REQ, "0.0" },
    { "inferoef", OPT_BOOL, "Infer oxygen extraction fraction", OPT_NONREQ, "" },
    { "inferr2p", OPT_BOOL, "Infer OEF-related R2 rate rather than OEF", OPT_NONREQ, "" },
    { "inferdbv", OPT_BOOL, "Infer deoxygenated blood volume", OPT_NONREQ, "" },
    { "inferr2t", OPT_BOOL, "Infer T2 relaxation rate of tissue", OPT_NONREQ, "" },
    { "infersig0", OPT_BOOL, "Infer baseline signal", OPT_NONREQ, "" },
    { "inferhct", OPT_BOOL, "Infer heamatocrit value", OPT_NONREQ, "" },
    { "inferr2e", OPT_BOOL, "Infer T2 relaxation rate of CSF", OPT_NONREQ, "" },
    { "inferdf", OPT_BOOL, "Infer CSF frequency shift df", OPT_NONREQ, "" },
    { "inferlam", OPT_BOOL, "Infer CSF fractional volume", OPT_NONREQ, "" },
    { "incintra", OPT_BOOL, "Include intravascular signal", OPT_NONREQ, "" },
    { "inccsf", OPT_BOOL, "Include CSF signal", OPT_NONREQ, "" },
    { "ignore-t1", OPT_BOOL, "Ignore T1", OPT_NONREQ, "" },
    { "motion-narrowing", OPT_BOOL, "Use motional narrowing model for intravascular signal", OPT_NONREQ, "" },
    { "dbv", OPT_FLOAT, "Default deoxygenated blood volume fraction", OPT_NONREQ, "0.03" },
    { "r2t", OPT_FLOAT, "Default T2 relaxation rate of tissue (s^-1)", OPT_NONREQ, "11.5" },
    { "sig0", OPT_FLOAT, "Default signal offset", OPT_NONREQ, "500" },
    { "r2e", OPT_FLOAT, "Default T2 relaxation rate of CSF (s^-1)", OPT_NONREQ, "4.0" },
    { "hct", OPT_FLOAT, "Default haematocrit", OPT_NONREQ, "0.4" },
    { "df", OPT_FLOAT, "Default CSF frequency shift df", OPT_NONREQ, "5.0" },
    { "lam", OPT_FLOAT, "Default CSF fractional volume (if including CSF component)", OPT_NONREQ, "0.1" },
    { "b0", OPT_FLOAT, "Field strength", OPT_NONREQ, "3.0" },
    { "" },
};

void R2primeFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string R2primeFwdModel::ModelVersion() const
{
    string version = "fwdmodel_qbold.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void R2primeFwdModel::Initialize(FabberRunData &rundata)
{
    // Default parameter values. These are used as the prior mean
    // when inferring the parameter, or as the fixed value if not
    m_dbv = rundata.GetDoubleDefault("dbv", 0.03);
    m_r2t = rundata.GetDoubleDefault("r2t", 11.5);
    m_sig0 = rundata.GetDoubleDefault("sig0", 500);
    m_r2e = rundata.GetDoubleDefault("r2e", 4.0);
    m_hct = rundata.GetDoubleDefault("hct", 0.4);
    m_df = rundata.GetDoubleDefault("df", 5.0);
    m_lam = rundata.GetDoubleDefault("lam", 0.1);
    m_b0 = rundata.GetDoubleDefault("b0", 3.0);
    m_tc_factor = rundata.GetDoubleDefault("tc-factor", 1.5);
    
    // Inference flags
    m_infer_oef = rundata.GetBool("inferoef");
    m_infer_r2p = rundata.GetBool("inferr2p");
    m_infer_dbv = rundata.GetBool("inferdbv");
    m_infer_r2t = rundata.GetBool("inferr2t");
    m_infer_sig0  = rundata.GetBool("infersig0");
    m_infer_hct = rundata.GetBool("inferhct");
    m_infer_r2e = rundata.GetBool("inferr2e");
    m_infer_df  = rundata.GetBool("inferdf");
    m_infer_lam = rundata.GetBool("inferlam");

    if (m_infer_oef && m_infer_r2p) 
    {
        throw FabberRunDataError("Can't infer both OEF and R2p");
    }
    else if (!m_infer_oef && !m_infer_r2p)
    {
        throw FabberRunDataError("Must infer either OEF or R2p");
    }

    // Model options and inclusions
    m_inc_intra = rundata.GetBool("incintra");
    m_inc_csf = rundata.GetBool("inccsf");
    m_ignore_t1 = rundata.GetBool("ignore-t1");
    m_motion_narr = rundata.GetBool("motion-narrowing");

    // Motion narrowing implies an intravascular component
    m_inc_intra = m_inc_intra || m_motion_narr;

    // Inferring any CSF parameters implies a CSF component
    m_inc_csf = m_inc_csf || m_infer_lam || m_infer_df;

    // If not including a CSF component, set fractional CSF volume to zero
    if (!m_inc_csf) m_lam = 0.0;

    // First, read tau values, these will always be specified
    vector<double> tausv = rundata.GetDoubleList("tau");
    m_taus.ReSize(tausv.size());
    for (unsigned int i=0; i<tausv.size(); i++) m_taus(i+1) = tausv[i];

    // Then read TE values. There might be a sequence or just one
    vector<double> tes = rundata.GetDoubleList("te");
    m_tes.ReSize(tausv.size());
    if (tes.size() == 1) 
    {
        m_tes = tes[0];
    }
    else
    {
        for (unsigned int i=0; i<tes.size(); i++)
        {
            m_tes(i+1) = tes[i];
        }
    }

    // Read TR and TI
    m_tr = rundata.GetDoubleDefault("tr", 3.0);
    m_ti = rundata.GetDoubleDefault("ti", 0.0);

    // add information to the log
    LOG << "Inference using development model" << endl;
    LOG << "Using TR = " << m_tr << "s, and TI = " << m_ti << "s" << endl;
    for (int ii = 1; ii <= m_taus.Nrows(); ii++)
    {
        LOG << "    TE(" << ii << ") = " << m_tes(ii) << "    tau(" << ii << ") = " << m_taus(ii) << endl;
    }
    if (m_inc_intra)
    {
        if (m_motion_narr)
        {
            LOG << "Using two-compartment model with motion narrowing intravascular signal " << endl;
        }
        else
        {
            LOG << "Using two-compartment model with static (powder) intravascular signal " << endl;
        }
    }
    else
    {
        LOG << "Using single-compartment model (ignoring intravascular signal) " << endl;
    }
    if (m_infer_oef)
    {
        LOG << "Inferring on OEF " << endl;
    }
    if (m_infer_r2p)
    {
        LOG << "Inferring on R2p " << endl;
    }
    if (m_infer_dbv)
    {
        LOG << "Inferring on DBV " << endl;
    }
    if (m_infer_r2t)
    {
        LOG << "Inferring on R2/T2 of tissue " << endl;
    }
    if (m_infer_sig0)
    {
        LOG << "Inferring on scaling parameter sig0 " << endl;
    }
    if (m_infer_hct)
    {
        LOG << "Inferring on fractional hematocrit " << endl;
    }
    if (m_infer_r2e)
    {
        LOG << "Inferring on R2 of CSF " << endl;
    }
    if (m_infer_df)
    {
        LOG << "Inferring on CSF frequency shift dF " << endl;
    }
    if (m_infer_lam)
    {
        LOG << "Inferring on CSF volume fraction lambda " << endl;
    }
}

void R2primeFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    if (m_infer_sig0) params.push_back(Parameter(p++, "sig0", DistParams(m_sig0, 1e6), DistParams(m_sig0, 100), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_oef) params.push_back(Parameter(p++, "oef", DistParams(0.4, 10), DistParams(0.4, 10), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_r2p) params.push_back(Parameter(p++, "r2p", DistParams(4.0, 1e3), DistParams(4.0, 1e2), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_dbv) params.push_back(Parameter(p++, "dbv", DistParams(m_dbv, 10), DistParams(m_dbv, 10), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_r2t) params.push_back(Parameter(p++, "r2t", DistParams(m_r2t, 1e2), DistParams(m_r2t, 1e2), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_hct) params.push_back(Parameter(p++, "hct", DistParams(m_hct, 1e-3), DistParams(m_hct, 1e-3), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_r2e) params.push_back(Parameter(p++, "r2e", DistParams(m_r2e, 1e2), DistParams(m_r2e, 1e2), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_df) params.push_back(Parameter(p++, "df", DistParams(m_df, 1e4), DistParams(m_df, 100), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_lam) params.push_back(Parameter(p++, "lam", DistParams(m_lam, 10), DistParams(m_lam, 10), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
}

/**
 * Evaluate the quantitative BOLD model
 * 
 * Seems to largely follow He and Yablonskiy 2007: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3971521/
 * 
 */
void R2primeFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    result.ReSize(data.Nrows());

    // Default values for parameters - values are overridden if the parameter is inferred
    double sig0 = m_sig0;
    double oef = 0;         // If not inferred, derived from R2P
    double r2p = 0;         // If not inferred, derived from OEF
    double dbv = m_dbv;
    double r2t = m_r2t;
    double hct = m_hct;
    double r2e = m_r2e;
    double df = m_df;
    double lam = m_lam;

    // Assign values to parameters which are being inferred
    int p = 1;
    if (m_infer_sig0) sig0 = params(p++);
    if (m_infer_oef) oef = params(p++);
    if (m_infer_r2p) r2p = params(p++);
    if (m_infer_dbv) dbv = params(p++);
    if (m_infer_r2t) r2t = params(p++);
    if (m_infer_hct) hct = params(p++);
    if (m_infer_r2e) r2e = params(p++);
    if (m_infer_df) df = params(p++);
    if (m_infer_lam) lam = params(p++);

    if (m_infer_dbv)
    {
        // This bit makes sure the value for dbv isn't ridiculous
        if (dbv < 0.0001)
        {
            dbv = 0.0001;
        }
        else if (dbv > 0.5)
        {
            dbv = 0.5;
        }
    }

    // OEF and R2p are equivalent ways to describe the effect of
    // oxygen extraction on the T2 relaxation rate.
    //
    // Relationship below is based on Eq 2 in He and Yablonskiy 2007
    //
    // dw = gamma * 4/3 * pi * delta_chi0 * hct * oef * B0
    // r2p = dbv * dw

    double dw;          // 1/tc in He and Yablonskiy 2007
    
    // Conversion factor from OEF to dw
    double oef_dw_factor = GAMMA * 4/3 * M_PI * D_CHI0 * hct * m_b0;

    if (m_infer_oef)
    {
        dw = oef_dw_factor*oef;
        r2p = dw*dbv;
    }
    else
    {
        if (r2p < 0.01)
        {
            r2p = 0.01;
        }

        dw = r2p / dbv;
        oef = dw / oef_dw_factor;
    }

    // Characteristic time. Note that He and Yablonskiy 2007 define tc = 1/dw and
    // use regimes separated by 1.5tc - hence we default to a TC factor of 1.5.
    // Matt Cherukara has work demonstrating that 1.76tc is better, when this is 
    // published we will probably move the default factor to 1.76
    double tc = m_tc_factor/dw;

    // Cerebral blood volume fraction, calculated below
    double CBV;

    // Apparent lambda (CSF volume fraction - opposite way round from literature!)
    double lam0;

    // Simulated data doesn't have any T1 contrast
    if (m_ignore_t1)
    {
        lam0 = lam;
        CBV = dbv;
    }
    else
    {
        // This mostly follows He and Yablonskiy 2007 eqns 14-17

        // Here are some more constants we will need. Note they are
        // not the same as He and Yablonskiy!
        //
        // T1 times
        double T1t = 1.20;
        double T1e = 3.87;
        double T1b = 1.58;

        // Spin densities
        double nt = 0.723;
        double ne = 1.000;
        double nb = 0.775;

        // Calculate magnetizations in tissue, blood and extracellular space
        // Related to Eq 16 but not identical
        double mt = 1.0 - ( ( 2.0 - exp(-(m_tr-m_ti)/T1t) ) * exp(-m_ti/T1t) );
        double mb = 1.0 - ( ( 2.0 - exp(-(m_tr-m_ti)/T1b) ) * exp(-m_ti/T1b) );
        double me = 1.0 - ( ( 2.0 - exp(-(m_tr-m_ti)/T1e) ) * exp(-m_ti/T1e) );

        // Calculate tissue compartment weightings

        // Extracellular (CSF/ISF) (Eq 15)
        lam0 = (ne*me*lam) / ( (nt*mt*(1-lam)) + (ne*me*lam) );

        // Intravascular (Eq 17)
        CBV = nb*mb*(1-lam0)*dbv;
    }
    
    // loop through taus
    result.ReSize(m_taus.Nrows());
    for (int ii = 1; ii <= m_taus.Nrows(); ii++)
    {
        double tau = m_taus(ii);
        double te = m_tes(ii);

        // Static dephasing signal (will be used to make St and Se)
        double Ss;

        // Calculate tissue signal. This is always included
        if (tau < -tc)
        {
            Ss = exp(dbv + (r2p*tau));          // SDR model
        }
        else if (tau > tc)
        {
            Ss = exp(dbv - (r2p*tau));          // SDR model
        }
        else
        {
            Ss = exp(-0.3*pow(r2p*tau,2.0)/dbv);          // SDR model
        }

        // Add T2 effect to tissue compartment to get tissue signal
        double St = Ss*exp(-r2t*te);

        // Calculate intravascular signal if this is included
        double Sb = 0.0;
        if (m_motion_narr)
        {
            double td   = 0.0045067;       // (based on rc=2.6 um and D=1.5 um^2 / ms)
            double gm   = 2.67513e8;
            double dChi = ((0.27*oef) + 0.14)*1e-6;
            double G0   = (4/45)*hct*(1-hct)*pow((dChi*m_b0),2.0);
            double kk   = 0.5*pow(gm,2.0)*G0*pow(td,2.0);
            double R2b = 5.291;    // fixed value (Berman, 2017) 

            // motion narrowing model
            Sb = exp(-kk* ( (te/td) + pow((0.25 + (te/td)),0.5) + 1.5 - 
                            (2*pow((0.25 + (pow((te+tau),2.0)/td) ),0.5)) - 
                            (2*pow((0.25 + (pow((te-tau),2.0)/td) ),0.5)) ) );

            // T2 effect 
            Sb *= exp(-R2b*te); 
        }
        else if (m_inc_intra)
        {
            // Blood relaxation rate
            double R2b  = ( 4.5 + (16.4*hct)) + ( ((165.2*hct) + 55.7)*pow(oef,2.0) );
    
            // Linear model
            //double R2bp = (10.2 - ( 1.5*hct)) + ( ((136.9*hct) - 13.9)*pow(oef,2.0) );
            //Sb = exp(-R2b*te)*exp(-R2bp*abs(tau));

            // Powder model

            // Threshold constant (similar to tc but different?)
            double pp = 1.5*dw*tau;

            // Signum of tau
            double st = (tau > 0) - (tau < 0);

            // Complex version of the blood signal for powder model
            complex<double> Sbc; 
            if (abs(pp) > 1)
            {
                // large pp
                Sbc = 0.5*sqrt(M_PI/abs(pp))*exp(I*((pp/3.0)-(st*M_PI/4.0)));
            }
            else
            {
                // small pp
                Sbc = 1.0 - ((2.0/45.0)*pow(pp,2.0)) + ((8.0*I*pow(pp,3.0))/2835.0);
            }

            Sb = real(Sbc);

            // T2 effect
            Sb *= exp(-R2b*te);
        }

        // Calculate CSF signal if this is included
        double Se = 0.0;
        if (m_inc_csf)
        {
            // Complex version of the extracellular signal (may be unnecessary)
            complex<double> Sec = Ss*exp(-r2e*te)*exp(-2.0*I*M_PI*df*abs(tau));
            Se = abs(Sec);
        }

        // Add up the compartments
        result(ii) = sig0*(((1-CBV-lam0)*St) + (CBV*Sb) + (lam0*Se)); 
    }
}
