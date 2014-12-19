// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// diffCalcHarm
NumericVector diffCalcHarm(NumericVector idt, NumericMatrix pw);
RcppExport SEXP diveRsity_diffCalcHarm(SEXP idtSEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type idt(idtSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type pw(pwSEXP );
        NumericVector __result = diffCalcHarm(idt, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// genos2mat
NumericMatrix genos2mat(NumericMatrix mat, IntegerVector ip, NumericVector na);
RcppExport SEXP diveRsity_genos2mat(SEXP matSEXP, SEXP ipSEXP, SEXP naSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type mat(matSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type ip(ipSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type na(naSEXP );
        NumericMatrix __result = genos2mat(mat, ip, na);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// glbWCcpp
List glbWCcpp(IntegerVector hsum, NumericMatrix af, NumericVector indtyp);
RcppExport SEXP diveRsity_glbWCcpp(SEXP hsumSEXP, SEXP afSEXP, SEXP indtypSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< IntegerVector >::type hsum(hsumSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type indtyp(indtypSEXP );
        List __result = glbWCcpp(hsum, af, indtyp);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// myTab
NumericVector myTab(CharacterVector x);
RcppExport SEXP diveRsity_myTab(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP );
        NumericVector __result = myTab(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pwHCalc
List pwHCalc(NumericMatrix af, NumericVector sHarm, IntegerMatrix pw);
RcppExport SEXP diveRsity_pwHCalc(SEXP afSEXP, SEXP sHarmSEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sHarm(sHarmSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type pw(pwSEXP );
        List __result = pwHCalc(af, sHarm, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pwHt
List pwHt(NumericMatrix af, IntegerMatrix pw);
RcppExport SEXP diveRsity_pwHt(SEXP afSEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type pw(pwSEXP );
        List __result = pwHt(af, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pwTabMerge
List pwTabMerge(List hsum, NumericMatrix pw);
RcppExport SEXP diveRsity_pwTabMerge(SEXP hsumSEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type hsum(hsumSEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type pw(pwSEXP );
        List __result = pwTabMerge(hsum, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// pwWCcpp
List pwWCcpp(List hsum1, NumericMatrix af1, NumericVector indtyp1, IntegerMatrix pw);
RcppExport SEXP diveRsity_pwWCcpp(SEXP hsum1SEXP, SEXP af1SEXP, SEXP indtyp1SEXP, SEXP pwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type hsum1(hsum1SEXP );
        Rcpp::traits::input_parameter< NumericMatrix >::type af1(af1SEXP );
        Rcpp::traits::input_parameter< NumericVector >::type indtyp1(indtyp1SEXP );
        Rcpp::traits::input_parameter< IntegerMatrix >::type pw(pwSEXP );
        List __result = pwWCcpp(hsum1, af1, indtyp1, pw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// Tab
IntegerVector Tab(CharacterVector x);
RcppExport SEXP diveRsity_Tab(SEXP xSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP );
        IntegerVector __result = Tab(x);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// tabMerge
NumericVector tabMerge(List hsum);
RcppExport SEXP diveRsity_tabMerge(SEXP hsumSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< List >::type hsum(hsumSEXP );
        NumericVector __result = tabMerge(hsum);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// varFunc
List varFunc(NumericMatrix af, double sHarm);
RcppExport SEXP diveRsity_varFunc(SEXP afSEXP, SEXP sHarmSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type af(afSEXP );
        Rcpp::traits::input_parameter< double >::type sHarm(sHarmSEXP );
        List __result = varFunc(af, sHarm);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}