#include<cmath>

struct CastlesNode{
    static length_t safe_div(length_t n, length_t d){
        if (d == 0) return 0;
        else return n / d;
    }

    length_t edgeLengthChildOfRootSiblingLeaf, edgeLengthOtherwise;
    length_t leftEdgeLength, rightEdgeLength, siblingEdgeLength;

    CastlesNode(const CustomizedAnnotation &annot){
        length_t LR_SO_quartetCnt = annot.ab_cd.quartetCnt, LR_SO_sumInternal = annot.ab_cd.sumInternalLength;
        length_t LR_SO_sumL = annot.ab_cd.sumLengthD, LR_SO_sumR = annot.ab_cd.sumLengthC, LR_SO_sumS = annot.ab_cd.sumLengthB, LR_SO_sumO = annot.ab_cd.sumLengthA;
        length_t LS_RO_quartetCnt = annot.ac_bd.quartetCnt, LS_RO_sumInternal = annot.ac_bd.sumInternalLength;
        length_t LS_RO_sumL = annot.ac_bd.sumLengthD, LS_RO_sumR = annot.ac_bd.sumLengthB, LS_RO_sumS = annot.ac_bd.sumLengthC, LS_RO_sumO = annot.ac_bd.sumLengthA;
        length_t LO_RS_quartetCnt = annot.ad_bc.quartetCnt, LO_RS_sumInternal = annot.ad_bc.sumInternalLength;
        length_t LO_RS_sumL = annot.ad_bc.sumLengthB, LO_RS_sumR = annot.ad_bc.sumLengthD, LO_RS_sumS = annot.ad_bc.sumLengthC, LO_RS_sumO = annot.ad_bc.sumLengthA;

        length_t num_m_gts = LR_SO_quartetCnt;
        length_t num_n_gts = LS_RO_quartetCnt + LO_RS_quartetCnt;
        length_t p_est = (num_m_gts - 0.5 * (1 + num_n_gts)) / (num_n_gts + num_m_gts + 1);
        length_t d_est = -log(1 - p_est);

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
        length_t delta = ((lm_i > ln_i) ? safe_div(lm_i - ln_i, ln_i) : 1e-03);
        length_t l_est = ((ln_i > 0) ? 1 / 6.0 * (3 * delta + sqrt(3 * delta * (4 + 3 * delta))) * ln_i : 1e-06);
        length_t mu1_est = l_est / d_est;
        length_t l_a_est = ln_a + (mu1_est * (d_est - p_est) + (lm_a - ln_a) * (1 - 2 / 3.0 * (1 - p_est))) / (1 - 4 / 5.0 * (1 - p_est)) - l_est;
        length_t l_b_est = ln_b + (mu1_est * (d_est - p_est) + (lm_b - ln_b) * (1 - 2 / 3.0 * (1 - p_est))) / (1 - 4 / 5.0 * (1 - p_est)) - l_est;
        length_t l_c_est = ln_c - 1 / 3.0 * (2 - 1 / (p_est + 1)) * (lm_c - ln_c);
        length_t l_d_est = ln_d - 2 / 3.0 * (2 + 1 / p_est) * (lm_d - ln_d);

        edgeLengthChildOfRootSiblingLeaf = ((l_d_est > 0) ? l_d_est : 0); // TODO: replace zero with average terminal bl
        edgeLengthOtherwise = abs(l_est);
        leftEdgeLength = abs(l_a_est);
        rightEdgeLength = abs(l_b_est);
        siblingEdgeLength = abs(l_c_est);
    }
};