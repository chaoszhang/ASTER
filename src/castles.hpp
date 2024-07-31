#include<cmath>

struct CastlesNode{
    static length_t safe_div(length_t n, length_t d){
        if (d == 0) return 0.0;
        else return n / d;
    }

    static length_t lambertw(length_t x, length_t err = 1e-9){
        length_t L, U;
        if (x > exp(1.0)){
            L = 1.0;
            U = log(x);
        }
        else {
            L = -1;
            U = 1;
        }
        while (L + err < U && U - L > L * err){
            length_t M = (L + U) / 2.0;
            if (M * exp(M) < x) L = M;
            else U = M;
        }
        return L;
    }

    length_t edgeLengthChildOfRootSiblingLeaf, edgeLengthOtherwise;
    length_t leftEdgeLength, rightEdgeLength, siblingEdgeLength;

    CastlesNode(const CustomizedAnnotation &annot, length_t s = 1000, length_t ngts = 1000){
        length_t LR_SO_quartetCnt = annot.ab_cd.quartetCnt, LR_SO_sumInternal = annot.ab_cd.sumInternalLength;
        length_t LR_SO_sumL = annot.ab_cd.sumLengthD, LR_SO_sumR = annot.ab_cd.sumLengthC, LR_SO_sumS = annot.ab_cd.sumLengthB, LR_SO_sumO = annot.ab_cd.sumLengthA;
        length_t LS_RO_quartetCnt = annot.ac_bd.quartetCnt, LS_RO_sumInternal = annot.ac_bd.sumInternalLength;
        length_t LS_RO_sumL = annot.ac_bd.sumLengthD, LS_RO_sumR = annot.ac_bd.sumLengthB, LS_RO_sumS = annot.ac_bd.sumLengthC, LS_RO_sumO = annot.ac_bd.sumLengthA;
        length_t LO_RS_quartetCnt = annot.ad_bc.quartetCnt, LO_RS_sumInternal = annot.ad_bc.sumInternalLength;
        length_t LO_RS_sumL = annot.ad_bc.sumLengthB, LO_RS_sumR = annot.ad_bc.sumLengthD, LO_RS_sumS = annot.ad_bc.sumLengthC, LO_RS_sumO = annot.ad_bc.sumLengthA;

        length_t num_m_gts = LR_SO_quartetCnt;
        length_t num_n_gts = LS_RO_quartetCnt + LO_RS_quartetCnt;
        length_t p_est = (num_m_gts - 0.5 * (1.0 + num_n_gts)) / (num_n_gts + num_m_gts + 1);
        length_t d_est = -log(1.0 - p_est);

        length_t lm_i = safe_div(LR_SO_sumInternal, LR_SO_quartetCnt);
        length_t ln_i = safe_div(LS_RO_sumInternal + LO_RS_sumInternal, LS_RO_quartetCnt + LO_RS_quartetCnt);
        length_t lm_a = safe_div(LR_SO_sumL, LR_SO_quartetCnt);
        length_t ln_a = safe_div(LS_RO_sumL + LO_RS_sumL, LS_RO_quartetCnt + LO_RS_quartetCnt);
        length_t lm_b = safe_div(LR_SO_sumR, LR_SO_quartetCnt);
        length_t ln_b = safe_div(LS_RO_sumR + LO_RS_sumR, LS_RO_quartetCnt + LO_RS_quartetCnt);
        length_t lm_c = safe_div(LR_SO_sumS, LR_SO_quartetCnt);
        length_t ln_c = safe_div(LS_RO_sumS + LO_RS_sumS, LS_RO_quartetCnt + LO_RS_quartetCnt);
        length_t lm_d = safe_div(LR_SO_sumO, LR_SO_quartetCnt);
        length_t ln_d = safe_div(LS_RO_sumO + LO_RS_sumO, LS_RO_quartetCnt + LO_RS_quartetCnt);

        length_t delta = safe_div(lm_i - ln_i, ln_i + 1.0 / s), l_est;

        if (delta > 0 && ln_i > 0) l_est = (delta + lambertw(-1.0 / 3.0 * exp(-delta - 1.0) * (2.0 * delta + 3.0)) + 1.0) * ln_i;
        else{
            delta = 1e-03;
            length_t l_lambert = (delta + lambertw(-1.0 / 3.0 * exp(-delta - 1.0) * (2.0 * delta + 3.0)) + 1.0) * ln_i;
            length_t threshold = log10(ngts);
            length_t w_mean = threshold * d_est, w_lambert = 1.0 / (threshold * d_est);
            l_est = (w_mean * lm_i + w_lambert * l_lambert) / (w_mean + w_lambert);
        }
        length_t mu1_est = l_est / d_est;
        
        length_t l_a_est = ln_a + (mu1_est * (d_est - p_est) + (lm_a - ln_a) * (1.0 - 2.0 / 3.0 * (1.0 - p_est))) / (1.0 - 4.0 / 5.0 * (1.0 - p_est)) - l_est;
        length_t l_b_est = ln_b + (mu1_est * (d_est - p_est) + (lm_b - ln_b) * (1.0 - 2.0 / 3.0 * (1.0 - p_est))) / (1.0 - 4.0 / 5.0 * (1.0 - p_est)) - l_est;
        length_t l_c_est = ln_c - 1.0 / 3.0 * (2.0 - 1.0 / (p_est + 1.0)) * (lm_c - ln_c);
        length_t l_d_est = ln_d - 2.0 / 3.0 * (2.0 + 1.0 / p_est) * (lm_d - ln_d);

        edgeLengthChildOfRootSiblingLeaf = ((l_d_est > 0) ? l_d_est : 0); // TODO: replace zero with average terminal bl
        edgeLengthOtherwise = abs(l_est);
        leftEdgeLength = abs(l_a_est);
        rightEdgeLength = abs(l_b_est);
        siblingEdgeLength = abs(l_c_est);
    }
};