/*
 Copyright (C) 2022 Matthias Lungwitz

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_marketmodels_i
#define quantlib_marketmodels_i

%include common.i
%include linearalgebra.i
%include observer.i
%include volatilities.i

%{
using QuantLib::AbcdMathFunction;
%}

%shared_ptr(AbcdMathFunction);
class AbcdMathFunction {
  private:
    AbcdMathFunction();
  public:
    typedef Time argument_type;
    typedef Real result_type;

    AbcdMathFunction(Real a = 0.002,
                     Real b = 0.001, 
                     Real c = 0.16,
                     Real d = 0.0005);
    AbcdMathFunction(std::vector<Real> abcd);

    //! function value at time t: \f[ f(t) \f]
    Real operator()(Time t) const;

    //! time at which the function reaches maximum (if any)
    Time maximumLocation() const;

    //! maximum value of the function
    Real maximumValue() const;

    //! function value at time +inf: \f[ f(\inf) \f]
    Real longTermValue() const { return d_; }

    /*! first derivative of the function at time t
        \f[ f'(t) = [ (b-c*a) + (-c*b)*t) ] e^{-c*t} \f] */
    Real derivative(Time t) const;
    
    /*! indefinite integral of the function at time t
        \f[ \int f(t)dt = [ (-a/c-b/c^2) + (-b/c)*t ] e^{-c*t} + d*t \f] */
    Real primitive(Time t) const;
    
    /*! definite integral of the function between t1 and t2
        \f[ \int_{t1}^{t2} f(t)dt \f] */
    Real definiteIntegral(Time t1, Time t2) const;

    /*! Inspectors */
    Real a() const;
    Real b() const;
    Real c() const;
    Real d() const;
    const std::vector<Real>& coefficients();
    const std::vector<Real>& derivativeCoefficients();
    // the primitive is not abcd

    /*! coefficients of a AbcdMathFunction defined as definite
        integral on a rolling window of length tau, with tau = t2-t */
    std::vector<Real> definiteIntegralCoefficients(Time t,
                                                   Time t2) const;

    /*! coefficients of a AbcdMathFunction defined as definite
        derivative on a rolling window of length tau, with tau = t2-t */
    std::vector<Real> definiteDerivativeCoefficients(Time t,
                                                     Time t2) const;

    static void validate(Real a,
                         Real b,
                         Real c,
                         Real d);
};

%{
using QuantLib::AbcdFunction;
%}

%shared_ptr(AbcdFunction);
class AbcdFunction : public AbcdMathFunction {
  public:
    AbcdFunction(Real a = -0.06,
                 Real b =  0.17,
                 Real c =  0.54,
                 Real d =  0.17);

    //! maximum value of the volatility function
    Real maximumVolatility() const;

    //! volatility function value at time 0: \f[ f(0) \f]
    Real shortTermVolatility() const;

    //! volatility function value at time +inf: \f[ f(\inf) \f]
    Real longTermVolatility() const;

    /*! instantaneous covariance function at time t between T-fixing and
        S-fixing rates \f[ f(T-t)f(S-t) \f] */
    Real covariance(Time t, Time T, Time S) const;

    /*! integral of the instantaneous covariance function between
        time t1 and t2 for T-fixing and S-fixing rates
        \f[ \int_{t1}^{t2} f(T-t)f(S-t)dt \f] */
    Real covariance(Time t1, Time t2, Time T, Time S) const;

     /*! average volatility in [tMin,tMax] of T-fixing rate:
        \f[ \sqrt{ \frac{\int_{tMin}^{tMax} f^2(T-u)du}{tMax-tMin} } \f] */
    Real volatility(Time tMin, Time tMax, Time T) const;

    /*! variance between tMin and tMax of T-fixing rate:
        \f[ \frac{\int_{tMin}^{tMax} f^2(T-u)du}{tMax-tMin} \f] */
    Real variance(Time tMin, Time tMax, Time T) const;
    

    
    // INSTANTANEOUS
    /*! instantaneous volatility at time t of the T-fixing rate:
        \f[ f(T-t) \f] */
    Real instantaneousVolatility(Time t, Time T) const;

    /*! instantaneous variance at time t of T-fixing rate:
        \f[ f(T-t)f(T-t) \f] */
    Real instantaneousVariance(Time t, Time T) const;

    /*! instantaneous covariance at time t between T and S fixing rates:
        \f[ f(T-u)f(S-u) \f] */
    Real instantaneousCovariance(Time u, Time T, Time S) const;

    // PRIMITIVE
    /*! indefinite integral of the instantaneous covariance function at
        time t between T-fixing and S-fixing rates
        \f[ \int f(T-t)f(S-t)dt \f] */
    Real primitive(Time t, Time T, Time S) const;

};

%{
using QuantLib::AbcdSquared;
%}

%shared_ptr(AbcdSquared);
class AbcdSquared {
  public:
    typedef Real argument_type;
    typedef Real result_type;

    AbcdSquared(Real a, Real b, Real c, Real d, Time T, Time S);
    Real operator()(Time t) const;
};

%{
using QuantLib::EvolutionDescription;
%}

%shared_ptr(EvolutionDescription);

class EvolutionDescription {
  public:
    EvolutionDescription();
    explicit EvolutionDescription(
        const std::vector<Time>& rateTimes,
        const std::vector<Time>& evolutionTimes = std::vector<Time>(),
        const std::vector<std::pair<Size,Size> >& relevanceRates =
                                                    std::vector<range>());
    const std::vector<Time>& rateTimes() const;
    const std::vector<Time>& rateTaus() const;
    const std::vector<Time>& evolutionTimes() const;
    //const Matrix& effectiveStopTimes() const;
    const std::vector<Size>& firstAliveRate() const;
    const std::vector<std::pair<Size,Size> >& relevanceRates() const;
    Size numberOfRates() const;
    Size numberOfSteps() const;
};

%{
using QuantLib::MarketModel;
%}

%shared_ptr(MarketModel);
class MarketModel {
  private:
    MarketModel();
  public:
    virtual const std::vector<Rate>& initialRates() const;
    virtual const std::vector<Spread>& displacements() const;
    virtual const EvolutionDescription& evolution() const;
    virtual Size numberOfRates() const;
    virtual Size numberOfFactors() const;
    virtual Size numberOfSteps() const;
    virtual const Matrix& pseudoRoot(Size i) const;
    virtual const Matrix& covariance(Size i) const;
    virtual const Matrix& totalCovariance(Size endIndex) const;
    std::vector<Volatility> timeDependentVolatility(Size i) const;
};

%{
using QuantLib::MarketModelFactory;
%}

%shared_ptr(MarketModelFactory);
//! base class for market-model factories
class MarketModelFactory : public Observable {
  private:
    MarketModelFactory();
  public:
    virtual ext::shared_ptr<MarketModel> create(
                                          const EvolutionDescription&,
                                          Size numberOfFactors) const;
};

%{
using QuantLib::PiecewiseConstantCorrelation;
%}

%shared_ptr(PiecewiseConstantCorrelation);
class PiecewiseConstantCorrelation {
  private:
    PiecewiseConstantCorrelation();
  public:
    virtual const std::vector<Time>& times() const;
    virtual const std::vector<Time>& rateTimes() const;
    virtual const std::vector<Matrix>& correlations() const;
    virtual const Matrix& correlation(Size i) const;
    virtual Size numberOfRates() const;
};

%{
using QuantLib::ExponentialForwardCorrelation;
%}

%shared_ptr(ExponentialForwardCorrelation);
class ExponentialForwardCorrelation :
    public PiecewiseConstantCorrelation {
  public:
    ExponentialForwardCorrelation(const std::vector<Time>& rateTimes,
                                  Real longTermCorr = 0.5,
                                  Real beta = 0.2,
                                  Real gamma = 1.0,
                                  std::vector<Time> times = std::vector<Time>());
    const std::vector<Time>& times() const override;
    const std::vector<Time>& rateTimes() const override;
    const std::vector<Matrix>& correlations() const override;
    Size numberOfRates() const override;
};

%{
using QuantLib::AbcdVol;
%}

%shared_ptr(AbcdVol);
//! %Abcd-interpolated volatility structure
class AbcdVol : public MarketModel {
  public:
    AbcdVol(Real a,
            Real b,
            Real c,
            Real d,
            const std::vector<Real>& ks,
            const ext::shared_ptr<PiecewiseConstantCorrelation>& corr,
            const EvolutionDescription& evolution,
            Size numberOfFactors,
            const std::vector<Rate>& initialRates,
            const std::vector<Spread>& displacements);
    //! \name MarketModel interface
    //@{
    const std::vector<Rate>& initialRates() const override;
    const std::vector<Spread>& displacements() const override;
    const EvolutionDescription& evolution() const override;
    Size numberOfRates() const override;
    Size numberOfFactors() const override;
    Size numberOfSteps() const override;
    const Matrix& pseudoRoot(Size i) const override;
    //@}
};

#endif