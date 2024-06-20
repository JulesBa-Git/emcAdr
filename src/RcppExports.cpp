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
Rcpp::List GeneticAlgorithm(int epochs, int nbIndividuals, const DataFrame& ATCtree, const DataFrame& observations, int num_thread, bool diversity, double p_crossover, double p_mutation, int nbElite, int tournamentSize, double alpha, bool summary);
RcppExport SEXP _emcAdr_GeneticAlgorithm(SEXP epochsSEXP, SEXP nbIndividualsSEXP, SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP num_threadSEXP, SEXP diversitySEXP, SEXP p_crossoverSEXP, SEXP p_mutationSEXP, SEXP nbEliteSEXP, SEXP tournamentSizeSEXP, SEXP alphaSEXP, SEXP summarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< int >::type nbIndividuals(nbIndividualsSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    Rcpp::traits::input_parameter< bool >::type diversity(diversitySEXP);
    Rcpp::traits::input_parameter< double >::type p_crossover(p_crossoverSEXP);
    Rcpp::traits::input_parameter< double >::type p_mutation(p_mutationSEXP);
    Rcpp::traits::input_parameter< int >::type nbElite(nbEliteSEXP);
    Rcpp::traits::input_parameter< int >::type tournamentSize(tournamentSizeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type summary(summarySEXP);
    rcpp_result_gen = Rcpp::wrap(GeneticAlgorithm(epochs, nbIndividuals, ATCtree, observations, num_thread, diversity, p_crossover, p_mutation, nbElite, tournamentSize, alpha, summary));
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
// trueDistributionSizeThreeCocktail
Rcpp::List trueDistributionSizeThreeCocktail(const DataFrame& ATCtree, const DataFrame& observations, int beta, int num_thread);
RcppExport SEXP _emcAdr_trueDistributionSizeThreeCocktail(SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP betaSEXP, SEXP num_threadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    rcpp_result_gen = Rcpp::wrap(trueDistributionSizeThreeCocktail(ATCtree, observations, beta, num_thread));
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
// hyperparam_test_genetic_algorithm
void hyperparam_test_genetic_algorithm(int epochs, int nb_individuals, const DataFrame& ATCtree, const DataFrame& observations, int nb_test_desired, const std::vector<double>& mutation_rate, const std::vector<int>& nb_elite, const std::vector<double>& alphas, const std::string& path, int num_thread);
RcppExport SEXP _emcAdr_hyperparam_test_genetic_algorithm(SEXP epochsSEXP, SEXP nb_individualsSEXP, SEXP ATCtreeSEXP, SEXP observationsSEXP, SEXP nb_test_desiredSEXP, SEXP mutation_rateSEXP, SEXP nb_eliteSEXP, SEXP alphasSEXP, SEXP pathSEXP, SEXP num_threadSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type epochs(epochsSEXP);
    Rcpp::traits::input_parameter< int >::type nb_individuals(nb_individualsSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< int >::type nb_test_desired(nb_test_desiredSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type mutation_rate(mutation_rateSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type nb_elite(nb_eliteSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type path(pathSEXP);
    Rcpp::traits::input_parameter< int >::type num_thread(num_threadSEXP);
    hyperparam_test_genetic_algorithm(epochs, nb_individuals, ATCtree, observations, nb_test_desired, mutation_rate, nb_elite, alphas, path, num_thread);
    return R_NilValue;
END_RCPP
}
// analyse_resultats
void analyse_resultats(const std::vector<std::vector<int>>& reponses, const std::string& input_filename, int repetition, const DataFrame& ATCtree);
RcppExport SEXP _emcAdr_analyse_resultats(SEXP reponsesSEXP, SEXP input_filenameSEXP, SEXP repetitionSEXP, SEXP ATCtreeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type reponses(reponsesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type input_filename(input_filenameSEXP);
    Rcpp::traits::input_parameter< int >::type repetition(repetitionSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    analyse_resultats(reponses, input_filename, repetition, ATCtree);
    return R_NilValue;
END_RCPP
}
// analyse_resultats_2
void analyse_resultats_2(const std::vector<std::vector<int>>& reponses, const std::string& input_filename, int repetition, const DataFrame& ATCtree, bool have_solution);
RcppExport SEXP _emcAdr_analyse_resultats_2(SEXP reponsesSEXP, SEXP input_filenameSEXP, SEXP repetitionSEXP, SEXP ATCtreeSEXP, SEXP have_solutionSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type reponses(reponsesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type input_filename(input_filenameSEXP);
    Rcpp::traits::input_parameter< int >::type repetition(repetitionSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< bool >::type have_solution(have_solutionSEXP);
    analyse_resultats_2(reponses, input_filename, repetition, ATCtree, have_solution);
    return R_NilValue;
END_RCPP
}
// get_dissimilarity_from_list
std::vector<std::vector<double>> get_dissimilarity_from_list(const Rcpp::List& genetic_results, const DataFrame& ATCtree);
RcppExport SEXP _emcAdr_get_dissimilarity_from_list(SEXP genetic_resultsSEXP, SEXP ATCtreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type genetic_results(genetic_resultsSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dissimilarity_from_list(genetic_results, ATCtree));
    return rcpp_result_gen;
END_RCPP
}
// get_dissimilarity
std::vector<std::vector<double>> get_dissimilarity(const std::string& filename, const DataFrame& ATCtree, bool normalization);
RcppExport SEXP _emcAdr_get_dissimilarity(SEXP filenameSEXP, SEXP ATCtreeSEXP, SEXP normalizationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type ATCtree(ATCtreeSEXP);
    Rcpp::traits::input_parameter< bool >::type normalization(normalizationSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dissimilarity(filename, ATCtree, normalization));
    return rcpp_result_gen;
END_RCPP
}
// get_answer_class
Rcpp::DataFrame get_answer_class(const std::string& filename, const std::vector<std::string>& answer);
RcppExport SEXP _emcAdr_get_answer_class(SEXP filenameSEXP, SEXP answerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type answer(answerSEXP);
    rcpp_result_gen = Rcpp::wrap(get_answer_class(filename, answer));
    return rcpp_result_gen;
END_RCPP
}
// ATC_idx_to_string
std::vector<std::vector<std::string>> ATC_idx_to_string(const std::vector<std::vector<int>>& patients, const std::vector<std::string>& ATCName);
RcppExport SEXP _emcAdr_ATC_idx_to_string(SEXP patientsSEXP, SEXP ATCNameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<int>>& >::type patients(patientsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type ATCName(ATCNameSEXP);
    rcpp_result_gen = Rcpp::wrap(ATC_idx_to_string(patients, ATCName));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma
arma::vec dmvnrm_arma(const arma::mat& X, const arma::rowvec& mean, const arma::mat& sigma_k, const arma::vec& w);
RcppExport SEXP _emcAdr_dmvnrm_arma(SEXP XSEXP, SEXP meanSEXP, SEXP sigma_kSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma_k(sigma_kSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma(X, mean, sigma_k, w));
    return rcpp_result_gen;
END_RCPP
}
// FWD_EM
Rcpp::List FWD_EM(const arma::mat& X, int K, double eps, const arma::vec& w_i, int max_steps);
RcppExport SEXP _emcAdr_FWD_EM(SEXP XSEXP, SEXP KSEXP, SEXP epsSEXP, SEXP w_iSEXP, SEXP max_stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type w_i(w_iSEXP);
    Rcpp::traits::input_parameter< int >::type max_steps(max_stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(FWD_EM(X, K, eps, w_i, max_steps));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_emcAdr_ATCtoNumeric", (DL_FUNC) &_emcAdr_ATCtoNumeric, 2},
    {"_emcAdr_histogramToDitribution", (DL_FUNC) &_emcAdr_histogramToDitribution, 1},
    {"_emcAdr_incorporateOustandingRRToDistribution", (DL_FUNC) &_emcAdr_incorporateOustandingRRToDistribution, 2},
    {"_emcAdr_EMC", (DL_FUNC) &_emcAdr_EMC, 11},
    {"_emcAdr_DistributionApproximation", (DL_FUNC) &_emcAdr_DistributionApproximation, 10},
    {"_emcAdr_GeneticAlgorithm", (DL_FUNC) &_emcAdr_GeneticAlgorithm, 12},
    {"_emcAdr_trueDistributionSizeTwoCocktail", (DL_FUNC) &_emcAdr_trueDistributionSizeTwoCocktail, 4},
    {"_emcAdr_trueDistributionSizeThreeCocktail", (DL_FUNC) &_emcAdr_trueDistributionSizeThreeCocktail, 4},
    {"_emcAdr_MetricCalc", (DL_FUNC) &_emcAdr_MetricCalc, 7},
    {"_emcAdr_computeMetrics", (DL_FUNC) &_emcAdr_computeMetrics, 4},
    {"_emcAdr_hyperparam_test_genetic_algorithm", (DL_FUNC) &_emcAdr_hyperparam_test_genetic_algorithm, 10},
    {"_emcAdr_analyse_resultats", (DL_FUNC) &_emcAdr_analyse_resultats, 4},
    {"_emcAdr_analyse_resultats_2", (DL_FUNC) &_emcAdr_analyse_resultats_2, 5},
    {"_emcAdr_get_dissimilarity_from_list", (DL_FUNC) &_emcAdr_get_dissimilarity_from_list, 2},
    {"_emcAdr_get_dissimilarity", (DL_FUNC) &_emcAdr_get_dissimilarity, 3},
    {"_emcAdr_get_answer_class", (DL_FUNC) &_emcAdr_get_answer_class, 2},
    {"_emcAdr_ATC_idx_to_string", (DL_FUNC) &_emcAdr_ATC_idx_to_string, 2},
    {"_emcAdr_dmvnrm_arma", (DL_FUNC) &_emcAdr_dmvnrm_arma, 4},
    {"_emcAdr_FWD_EM", (DL_FUNC) &_emcAdr_FWD_EM, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_emcAdr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
