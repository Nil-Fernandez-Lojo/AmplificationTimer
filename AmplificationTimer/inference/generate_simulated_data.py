import numpy as np
import copy

from AmplificationTimerObjects import SNV, Position

def coverage(segment, rho, nrpcc):
    cn = segment.major_cn + segment.minor_cn
    n_reads = round(nrpcc * (cn + ((1 - rho) / rho) * segment.get_ploidy_healthy()))
    return n_reads

def generate_simulated_data_under_prior(sample_,
                                     amplification_,
                                     model_idx,
                                     t=None,
                                     mu=None,
                                     min_reads_detect_SNV=0,
                                     subclonality_modelled_simulation=True,
                                     nrpcc='same'):
    # TODO: usage of constant nrpcc while actually should sample from distribution with that mean
    # TODO: put rest of attributes SNV
    # TODO: keep APOBEC like mutations?
    # TODO: actually use another prior for mutation rate

    # Some preprocessing
    sample = copy.deepcopy(sample_)
    amplification = copy.deepcopy(amplification_)

    rho = sample.clinical_data['purity']
    for segment in sample.segments:
        segment.snvs = []
    for segment in amplification.segments:
        segment.snvs = []

    if nrpcc == "same":
        nrpcc = sample.clinical_data["nrpcc"]
        print("nrpcc", nrpcc)

    fSNV = sample.subclonal_structure["n_snvs"].values / sample.subclonal_structure["n_snvs"].sum()
    subclone_ccf = sample.subclonal_structure["fraction_cancer_cells"].values

    # Sample latent variables
    if mu is None:
        mu = np.random.beta(1, 10 ** 6)
    if t is None:
        t = np.random.uniform()

    if model_idx == 3:
        u = np.random.uniform()
    elif model_idx == 4:
        u = np.random.uniform(size=len(amplification.segments))
    else:
        u = None

    # generate SNVs on 1+0 and 1+1 segments with NRPCC
    for segment in sample.segments:
        if segment.major_cn == 1 and segment.minor_cn <= 1:
            if segment.minor_cn == 0:
                p = mu
                q = rho / (rho + segment.get_ploidy_healthy() * (1 - rho))
            else:
                p = 2 * mu
                q = rho / (2 * rho + segment.get_ploidy_healthy() * (1 - rho))
            D = coverage(segment, rho, nrpcc)
            n_SNVs = np.random.binomial(n=segment.get_length(), p=p)

            for i in range(n_SNVs):
                d = np.random.binomial(n=D, p=q)
                # SNV detected
                if d >= min_reads_detect_SNV:
                    # TODO: sample position at random in amplified segments
                    segment.add_snv(SNV(segment.chromosome,
                                        Position(segment.chromosome, segment.start.position + 1, sample.config),
                                        D,
                                        D - d,
                                        "N",
                                        "N",
                                        sample.config,
                                        previous_ref_base='N',
                                        next_ref_base='N'))

    # generate SNVs on Amplified segments with NRPCC
    for i, segment in enumerate(amplification.segments):
        include_u = segment.minor_cn >= 2
        p_all = mu * t
        if model_idx in [1, 2]:
            p_1 = segment.major_cn * mu * (1 - t)
        elif model_idx == 3:
            p_1 = segment.major_cn * mu * (1 - t) + segment.minor_cn * mu * (u ** include_u)
        else:
            p_1 = segment.major_cn * mu * (1 - t) + segment.minor_cn * mu * (u[i] ** include_u)

        q_clonal_1 = rho / ((segment.major_cn + segment.minor_cn) * rho + segment.get_ploidy_healthy() * (
                    1 - rho))  # 1 copy clonal
        q_clonal_all = q_clonal_1 * segment.major_cn

        if subclonality_modelled_simulation and len(fSNV) > 1:
            q_subclonal_1 = q_clonal_1 * subclone_ccf[1:]

            p_all = p_all * fSNV[0]
            p_1 = p_1 * fSNV[0]
            p_1_subclonal = mu * fSNV[1:]
            p_0 = 1 - (p_all + p_1 + sum(p_1_subclonal))
            p = np.concatenate(([p_all, p_1], p_1_subclonal, [p_0]))
        else:
            q_subclonal_1 = []

            p_0 = 1 - p_all - p_1
            p = np.array([p_all, p_1, p_0])

        n = np.random.multinomial(segment.get_length(), pvals=p)
        D = coverage(segment, rho, nrpcc)
        q = np.concatenate(([q_clonal_all, q_clonal_1], q_subclonal_1))

        for snv_type_idx in range(len(n) - 1):
            for snv_idx in range(n[snv_type_idx]):
                d = np.random.binomial(n=D, p=q[snv_type_idx])
                if d >= min_reads_detect_SNV:
                    segment.add_snv(SNV(segment.chromosome,
                                        Position(segment.chromosome, segment.start.position + 1, sample.config),
                                        D - d,
                                        d,
                                        "N",
                                        "N",
                                        sample.config,
                                        previous_ref_base='N',
                                        next_ref_base='N'))

    amplification.mutation_rate = sample.get_mutation_rate(sample.segments)

    return t,u,amplification
