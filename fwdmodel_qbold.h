/**
 * fwdmodel_qbold_R2p.h - ASE qBOLD curve fitting model
 *
 * Matthew Cherukara, IBME
 *
 * Copyright (C) 2019 University of Oxford  
 */
#pragma once

#include <fabber_core/fwdmodel.h>

#include <newmat.h>

#include <string>
#include <vector>

class R2primeFwdModel : public FwdModel 
{
public:
    static FwdModel* NewInstance();

    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;
    void Initialize(ArgsType &args);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

protected:
    // Scan Parameters
    double m_tr;
    double m_ti;
    NEWMAT::ColumnVector m_taus;
    NEWMAT::ColumnVector m_tes;

    // Default parameter values
    double m_dbv;
    double m_r2t;
    double m_sig0;
    double m_hct;
    double m_r2e;
    double m_df;
    double m_lam;
    double m_b0;
    double m_tc_factor;

    // Which parameters will we infer on
    bool m_infer_oef;
    bool m_infer_r2p;
    bool m_infer_dbv;
    bool m_infer_r2t;
    bool m_infer_sig0;
    bool m_infer_hct;
    bool m_infer_r2e;
    bool m_infer_df;
    bool m_infer_lam;

    // Model options defining what signal contributions are included in the moodel
    bool m_motion_narr;
    bool m_inc_intra;
    bool m_inc_csf;
    bool m_ignore_t1;

private:
    static FactoryRegistration<FwdModelFactory, R2primeFwdModel> registration;
};
