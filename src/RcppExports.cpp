// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ATCtoNumeric
void ATCtoNumeric(DataFrame& patients, const DataFrame& tree);
RcppExport SEXP _emcAdr_ATCtoNumeric(SEXP patientsSEXP, SEXP treeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame& >::type patients(patientsSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type tree(treeSEXP);
    ATCtoNumeric(patients, tree);
    return R_NilValue;
END_RCPP
}
// histogramToDitribution
Rcpp::NumericVector histogramToDitribution(const std::vector<int>& vec);
RcppExport SEXP _emcAdr_histogramToDitribution(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(histogramToDitribution(vec));
    return rcpp_result_gen;
END_RCPP
}
// incorporateOustandingRRToDistribution
Rcpp::NumericVector incorporateOustandingRRToDistribution(const std::vector<double>& outstandingRR, int RRmax);
RcppExport SEXP _emcAdr_incorporateOustandingRRToDistribution(SEXP outstandingRRSEXP, SEXP RRmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type outstandingRR(outstandingRRSEXP);
    Rcpp::traits::input_parameter< int >::type RRmax(RRmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(incorporateOustandingRRToDistribution(outstandingRR, RRmax));
    return rcpp_result_gen;
END_RCPP
}
// EMC
Rcpp::List EMC(int n, const DataFrame& ATCtree, const DataFrame& observations, double P_type1, double P_type2, double P_crossover, int nbIndividuals, int nbResults, double alpha, Rcpp::Nullable<Rcpp::List> startingIndividuals, Rcpp::Nullable<Rcpp::NumericVector> startingTemperatures);
RcppExport SEXP _emcAdr_EMC(SEXP nSEXP, SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP P_type1SEXP, SEXP P_type2SEXP, SEXP P_crossoverSEXP, SEXP nbIndividualsSEXP, SEXP nbResultsSEXP, SEXP alphaSEXP, SEXP startingIndividualsSEXP, SEXP startingTemperaturesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< double >::type P_type1(P_type1SEXP);
    Rcpp::traits::input_parameter< double >::type P_type2(P_type2SEXP);
    Rcpp::traits::input_parameter< double >::type P_crossover(P_crossoverSEXP);
    Rcpp::traits::input_parameter< int >::type nbIndividuals(nbIndividualsSEXP);
    Rcpp::traits::input_parameter< int >::type nbResults(nbResultsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type startingIndividuals(startingIndividualsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type startingTemperatures(startingTemperaturesSEXP);
    rcpp_result_gen = Rcpp::wrap(EMC(n, ATCtree, observations, P_type1, P_type2, P_crossover, nbIndividuals, nbResults, alpha, startingIndividuals, startingTemperatures));
    return rcpp_result_gen;
END_RCPP
}
// DistributionApproximation
Rcpp::List DistributionApproximation(int epochs, const DataFrame& ATCtree, const DataFrame& observations, int temperature, int nbResults, int Smax, double p_type1, int beta, int max_Metric, int num_thread);
RcppExport SEXP _emcAdr_DistributionApproximation(SEXP epochsSEXP, SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP temperatureSEXP, SEXP nbResultsSEXP, SEXP SmaxSEXP, SEXP p_type1SEXP, SEXP betaSEXP, SEXP max_MetricSEXP, SEXP num_threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type temperature(temperatureSEXP);
    Rcpp::traits::input_parameter< int >::type nbResults(nbResultsSEXP);
    Rcpp::traits::input_parameter< int >::type Smax(SmaxSEXP);
    Rcpp::traits::input_parameter< double >::type p_type1(p_type1SEXP);
    Rcpp::traits::input_parameter< int >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type max_Metric(max_MetricSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    rcpp_result_gen = Rcpp::wrap(DistributionApproximation(epochs, ATCtree, observations, temperature, nbResults, Smax, p_type1, beta, max_Metric, num_thread));
    return rcpp_result_gen;
END_RCPP
}
// GeneticAlgorithm
Rcpp::List GeneticAlgorithm(int epochs, int nbIndividuals, const DataFrame& ATCtree, const DataFrame& observations, double p_crossover, double p_mutation, int nbElite, int tournamentSize);
RcppExport SEXP _emcAdr_GeneticAlgorithm(SEXP epochsSEXP, SEXP nbIndividualsSEXP, SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP p_crossoverSEXP, SEXP p_mutationSEXP, SEXP nbEliteSEXP, SEXP tournamentSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< int >::type nbIndividuals(nbIndividualsSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< double >::type p_crossover(p_crossoverSEXP);
    Rcpp::traits::input_parameter< double >::type p_mutation(p_mutationSEXP);
    Rcpp::traits::input_parameter< int >::type nbElite(nbEliteSEXP);
    Rcpp::traits::input_parameter< int >::type tournamentSize(tournamentSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(GeneticAlgorithm(epochs, nbIndividuals, ATCtree, observations, p_crossover, p_mutation, nbElite, tournamentSize));
    return rcpp_result_gen;
END_RCPP
}
// trueDistributionSizeTwoCocktail
Rcpp::List trueDistributionSizeTwoCocktail(const DataFrame& ATCtree, const DataFrame& observations, int beta, int num_thread);
RcppExport SEXP _emcAdr_trueDistributionSizeTwoCocktail(SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP betaSEXP, SEXP num_threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    rcpp_result_gen = Rcpp::wrap(trueDistributionSizeTwoCocktail(ATCtree, observations, beta, num_thread));
    return rcpp_result_gen;
END_RCPP
}
// MetricCalc
std::vector<double> MetricCalc(const std::vector<int>& cocktail, const std::vector<int>& ATClength, const std::vector<int>& upperBounds, const std::vector<std::vector<int>>& observationsMedication, const Rcpp::LogicalVector& observationsADR, int ADRCount, int num_thread);
RcppExport SEXP _emcAdr_MetricCalc(SEXP cocktailSEXP, SEXP ATClengthSEXP, SEXP upperBoundsSEXP, SEXP observationsMedicationSEXP, SEXP observationsADRSEXP, SEXP ADRCountSEXP, SEXP num_threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type cocktail(cocktailSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type ATClength(ATClengthSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type upperBounds(upperBoundsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type observationsMedication(observationsMedicationSEXP);
    Rcpp::traits::input_parameter< const Rcpp::LogicalVector& >::type observationsADR(observationsADRSEXP);
    Rcpp::traits::input_parameter< int >::type ADRCount(ADRCountSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    rcpp_result_gen = Rcpp::wrap(MetricCalc(cocktail, ATClength, upperBounds, observationsMedication, observationsADR, ADRCount, num_thread));
    return rcpp_result_gen;
END_RCPP
}
// computeMetrics
Rcpp::DataFrame computeMetrics(const Rcpp::DataFrame& df, const DataFrame& ATCtree, const DataFrame& observations, int num_thread);
RcppExport SEXP _emcAdr_computeMetrics(SEXP dfSEXP, SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP num_threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::DataFrame& >::type df(dfSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMetrics(df, ATCtree, observations, num_thread));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_emcAdr_ATCtoNumeric", (DL_FUNC) &_emcAdr_ATCtoNumeric, 2},
    {"_emcAdr_histogramToDitribution", (DL_FUNC) &_emcAdr_histogramToDitribution, 1},
    {"_emcAdr_incorporateOustandingRRToDistribution", (DL_FUNC) &_emcAdr_incorporateOustandingRRToDistribution, 2},
    {"_emcAdr_EMC", (DL_FUNC) &_emcAdr_EMC, 11},
    {"_emcAdr_DistributionApproximation", (DL_FUNC) &_emcAdr_DistributionApproximation, 10},
    {"_emcAdr_GeneticAlgorithm", (DL_FUNC) &_emcAdr_GeneticAlgorithm, 8},
    {"_emcAdr_trueDistributionSizeTwoCocktail", (DL_FUNC) &_emcAdr_trueDistributionSizeTwoCocktail, 4},
    {"_emcAdr_MetricCalc", (DL_FUNC) &_emcAdr_MetricCalc, 7},
    {"_emcAdr_computeMetrics", (DL_FUNC) &_emcAdr_computeMetrics, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_emcAdr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
