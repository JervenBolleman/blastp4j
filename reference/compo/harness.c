#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/composition_adjustment/composition_adjustment.h>
#include <algo/blast/composition_adjustment/compo_mode_condition.h>
#include <algo/blast/composition_adjustment/nlm_linear_algebra.h>
#include "freqratios.inc"
static Int4 BLAST_Gcd(Int4 a, Int4 b) { Int4 c; b = ABS(b); if (b > a) c=a, a=b, b=c; while (b != 0) { c = a%b; a = b; b = c; } return a; }
static double
#include "nr.inc"
/* Blast_KarlinLambdaNR on the probabilities, as s_CalcLambda does */
static double calc_lambda(double probs[], int min_score, int max_score, double lambda0) {
    int i, d, itn; double avg = 0.0; double *sprob = &probs[-min_score];
    for (i = 0; i < max_score - min_score + 1; i++) avg += (min_score + i) * probs[i];
    if (avg >= 0.) return -1.0;
    if (min_score >= 0 || max_score <= 0) return -1.;
    for (i = 1, d = -min_score; i <= max_score - min_score && d > 1; ++i) if (sprob[i + min_score] != 0.0) d = BLAST_Gcd(d, i);
    return NlmKarlinLambdaNR(sprob, d, min_score, max_score, lambda0, 1.e-5, 20, 20 + 17, &itn);
}
static double rd(void) { double v; if (scanf("%la", &v) != 1) exit(2); return v; }
static int ri(void) { int v; if (scanf("%d", &v) != 1) exit(2); return v; }
int main(void) {
    int cases = ri();
    Blast_CompositionWorkspace *ws = Blast_CompositionWorkspaceNew();
    Blast_CompositionWorkspaceInit(ws, "BLOSUM62");
    double lambda = rd();
    Blast_MatrixInfo *info = Blast_MatrixInfoNew(28, 28, 0);
    info->matrixName = "BLOSUM62";
    info->ungappedLambda = lambda;
    for (int i = 0; i < 28; i++) for (int j = 0; j < 28; j++) info->startFreqRatios[i][j] = FREQRATIOS[i][j];
    Blast_Int4MatrixFromFreq(info->startMatrix, 28, info->startFreqRatios, lambda);
    int **matrix = Nlm_Int4MatrixNew(28, 28);
    for (int c = 0; c < cases; c++) {
        Blast_AminoAcidComposition q, s; int qlen, slen; EMatrixAdjustRule rule; double pvalue = -1, ratio = 0;
        for (int i = 0; i < 28; i++) q.prob[i] = rd(); q.numTrueAminoAcids = ri(); qlen = ri();
        for (int i = 0; i < 28; i++) s.prob[i] = rd(); s.numTrueAminoAcids = ri(); slen = ri();
        for (int i = 0; i < 28; i++) for (int j = 0; j < 28; j++) matrix[i][j] = 0;
        int status = Blast_AdjustScores(matrix, &q, qlen, &s, slen, info, eCompositionMatrixAdjust, 20, ws, &rule, calc_lambda, &pvalue, 0, &ratio);
        printf("%d %d %a\n", status, (int) rule, ratio);
        for (int i = 0; i < 28; i++) { for (int j = 0; j < 28; j++) printf("%d ", matrix[i][j]); printf("\n"); }
    }
    return 0;
}
