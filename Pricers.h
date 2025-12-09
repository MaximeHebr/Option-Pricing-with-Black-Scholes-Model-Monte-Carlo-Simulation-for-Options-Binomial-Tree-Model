#ifndef PRICERS_H
#define PRICERS_H

#include "Options.h"
#include "Market.h"
#include <memory>
#include <vector>
#include <random>
#include <cmath>
#include <numbers>

// Classe de base pour les pricers
class Pricer {
protected:
    std::shared_ptr<Option> option_; // option à pricer
    std::shared_ptr<MarketData> marketData_; // données de marché associées à cette option 
    
public:
    Pricer(std::shared_ptr<Option> option, std::shared_ptr<MarketData> marketData)
        : option_(option), marketData_(marketData) {}
    
    virtual ~Pricer() = default;
    virtual double calculatePrice() = 0;
    virtual std::string getMethodName() const = 0;
};

// Black-Scholes pricer
class BlackScholesPricer : public Pricer {
private:
    mutable double d1_, d2_;
    mutable bool calculated_; // pour savoir si d1 et d2 ont été calculés
    
    static double normalCDF(double x) {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
    }

    static double normalPDF(double x) {
        return (1.0 / std::sqrt(2.0 * std::numbers::pi)) * std::exp(-0.5 * x * x); 
    }
    
    void calculateD1D2() const {
        if (calculated_) return;
        
        double S = marketData_->getSpotPrice();
        double K = option_->getStrike();
        double T = option_->getMaturity();
        double r = marketData_->getRiskFreeRate();
        double q = marketData_->getDividendYield();
        double sigma = marketData_->getVolatility();
        
        double sigmaSqrtT = sigma * std::sqrt(T);
        d1_ = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / sigmaSqrtT;
        d2_ = d1_ - sigmaSqrtT;
        calculated_ = true;
    }
    
public:
    BlackScholesPricer(std::shared_ptr<Option> option, std::shared_ptr<MarketData> marketData)
        : Pricer(option, marketData), calculated_(false) {
        if (!option->isEuropean()) {
            throw std::invalid_argument("Black-Scholes only works for European options");// on vérifie que l'option est européenne car le modèle de Black-Scholes ne s'applique qu'aux options européennes.
        }
    }
    
    double calculatePrice() override {
        calculateD1D2();
        
        double S = marketData_->getSpotPrice();
        double K = option_->getStrike();
        double T = option_->getMaturity();
        double r = marketData_->getRiskFreeRate();
        double q = marketData_->getDividendYield();
        
        double Nd1 = normalCDF(d1_);
        double Nd2 = normalCDF(d2_);
        double discountFactor = std::exp(-r * T);
        double dividendFactor = std::exp(-q * T);
        
        if (option_->isCall()) {
            return S * dividendFactor * Nd1 - K * discountFactor * Nd2; // formule de Black-Scholes pour un call européen
        } else {
            return K * discountFactor * normalCDF(-d2_) - S * dividendFactor * normalCDF(-d1_); // formule de Black-Scholes pour un put européen
        }
    }
    
    PricingResult calculateFullResult() {
        calculateD1D2();
        PricingResult result;
        
        double S = marketData_->getSpotPrice();
        double K = option_->getStrike();
        double T = option_->getMaturity();
        double r = marketData_->getRiskFreeRate();
        double q = marketData_->getDividendYield();
        double sigma = marketData_->getVolatility();
        double dividendFactor = std::exp(-q * T);
        double discountFactor = std::exp(-r * T);
        
        result.price = calculatePrice();
        
        // Delta
        if (option_->isCall()) {
            result.delta = dividendFactor * normalCDF(d1_);
        } else {
            result.delta = -dividendFactor * normalCDF(-d1_);
        }
        
        // Gamma
        result.gamma = (dividendFactor * normalPDF(d1_)) / (S * sigma * std::sqrt(T));
        
        // Vega
        result.vega = S * dividendFactor * normalPDF(d1_) * std::sqrt(T) / 100.0;
        
        // Theta
        double term1 = -(S * dividendFactor * normalPDF(d1_) * sigma) / (2.0 * std::sqrt(T));
        if (option_->isCall()) {
            double term2 = q * S * dividendFactor * normalCDF(d1_);
            double term3 = -r * K * discountFactor * normalCDF(d2_);
            result.theta = (term1 + term2 + term3) / 365.0;
        } else {
            double term2 = -q * S * dividendFactor * normalCDF(-d1_);
            double term3 = r * K * discountFactor * normalCDF(-d2_);
            result.theta = (term1 + term2 + term3) / 365.0;
        }
        
        // Rho
        if (option_->isCall()) {
            result.rho = K * T * discountFactor * normalCDF(d2_) / 100.0;
        } else {
            result.rho = -K * T * discountFactor * normalCDF(-d2_) / 100.0;
        }
        
        return result;
    }
    
    std::string getMethodName() const override { return "Black-Scholes"; }
};

// Binomial tree pricer
class BinomialTreePricer : public Pricer {
private:
    int numSteps_; // nombre d'étapes dans l'arbre binomial (plus il y en a, plus le prix est précis mais le calcul est long)
    std::vector<double> optionValues_; // valeurs de l'option aux nœuds de l'arbre
    double dt_, u_, d_, p_, discountFactor_; 
    
    void calculateTreeParameters() {
        double T = option_->getMaturity();
        double r = marketData_->getRiskFreeRate();
        double q = marketData_->getDividendYield();
        double sigma = marketData_->getVolatility();
        
        dt_ = T / numSteps_;
        u_ = std::exp(sigma * std::sqrt(dt_));
        d_ = 1.0 / u_;
        
        double growth = std::exp((r - q) * dt_);
        p_ = (growth - d_) / (u_ - d_);
        discountFactor_ = std::exp(-r * dt_);
    }
    
    double checkEarlyExercise(double spotPrice, double continuationValue) const {
        if (option_->isEuropean()) return continuationValue;
        double exerciseValue = option_->payoff(spotPrice);
        return std::max(exerciseValue, continuationValue);
    }
    
public:
    BinomialTreePricer(std::shared_ptr<Option> option, std::shared_ptr<MarketData> marketData, int numSteps = 100)
        : Pricer(option, marketData), numSteps_(numSteps) {
        if (numSteps <= 0) throw std::invalid_argument("Number of steps must be positive");
        optionValues_.reserve(numSteps + 1);
        calculateTreeParameters();
    }
    
    double calculatePrice() override {
        double S = marketData_->getSpotPrice();
        optionValues_.resize(numSteps_ + 1);
        
        // Terminal values
        for (int i = 0; i <= numSteps_; ++i) {
            double spotAtNode = S * std::pow(u_, i) * std::pow(d_, numSteps_ - i);// prix du sous-jacent au nœud final
            optionValues_[i] = option_->payoff(spotAtNode);// valeur de l'option au nœud final
        }
        
        // Backward induction
        for (int step = numSteps_ - 1; step >= 0; --step) {
            for (int i = 0; i <= step; ++i) {
                double continuationValue = discountFactor_ * (p_ * optionValues_[i + 1] + (1.0 - p_) * optionValues_[i]);// valeur actualisée attendue
                double spotAtNode = S * std::pow(u_, i) * std::pow(d_, step - i);
                optionValues_[i] = checkEarlyExercise(spotAtNode, continuationValue);// vérifier l'exercice anticipé pour les options américaines
            }
        }
        
        return optionValues_[0];
    }
    
    std::string getMethodName() const override { return "Binomial Tree (N=" + std::to_string(numSteps_) + ")"; }
};

// Monte Carlo simulation pricing
class MonteCarloPricer : public Pricer {
private:
    int numSimulations_;
    int numTimeSteps_;
    unsigned int seed_;// graine pour le générateur de nombres aléatoires
    mutable std::mt19937 generator_; // générateur de nombres aléatoires Mersenne Twister ( on en prend pas rand() car il est moins performant et moins fiable statistiquement)
    mutable std::normal_distribution<double> normalDist_; 
    
protected:
    double generateNormal() const { return normalDist_(generator_); }
    
    double simulatePath() const {
        double S = marketData_->getSpotPrice();
        double T = option_->getMaturity();
        double r = marketData_->getRiskFreeRate();
        double q = marketData_->getDividendYield();
        double sigma = marketData_->getVolatility();
        
        double dt = T / numTimeSteps_; // pas de temps
        double drift = (r - q - 0.5 * sigma * sigma) * dt; 
        double diffusion = sigma * std::sqrt(dt);
        
        double currentPrice = S;
        for (int i = 0; i < numTimeSteps_; ++i) {
            double Z = generateNormal();
            currentPrice *= std::exp(drift + diffusion * Z);
        }
        
        return currentPrice;
    }
    
    std::vector<double> simulateFullPath() const {
        std::vector<double> path;
        path.reserve(numTimeSteps_ + 1);
        
        double S = marketData_->getSpotPrice();
        double T = option_->getMaturity();
        double r = marketData_->getRiskFreeRate();
        double q = marketData_->getDividendYield();
        double sigma = marketData_->getVolatility();
        
        double dt = T / numTimeSteps_;
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double diffusion = sigma * std::sqrt(dt);
        
        path.push_back(S);
        double currentPrice = S;
        
        for (int i = 0; i < numTimeSteps_; ++i) {
            double Z = generateNormal();
            currentPrice *= std::exp(drift + diffusion * Z);
            path.push_back(currentPrice);
        }
        
        return path;
    }
    
public:
    MonteCarloPricer(std::shared_ptr<Option> option, std::shared_ptr<MarketData> marketData,
                     int numSimulations = 100000, int numTimeSteps = 252, unsigned int seed = std::random_device{}())
        : Pricer(option, marketData), numSimulations_(numSimulations), 
          numTimeSteps_(numTimeSteps), seed_(seed), generator_(seed), normalDist_(0.0, 1.0) {
        if (numSimulations <= 0 || numTimeSteps <= 0) {
            throw std::invalid_argument("Simulations and time steps must be positive"); // on vérifie que le nombre de simulations et le nombre de pas de temps sont positifs pour que le calcul soit cohérent.
        }
    }
    
    double calculatePrice() override {
        double r = marketData_->getRiskFreeRate();
        double T = option_->getMaturity();
        double discountFactor = std::exp(-r * T);
        
        bool isPathDependent = (dynamic_cast<AsianOption*>(option_.get()) != nullptr); // vérifier si l'option est path-dependent (asiatique)
        
        double sumPayoffs = 0.0;
        for (int i = 0; i < numSimulations_; ++i) {
            double payoff = 0.0;
            
            if (isPathDependent) { // option asiatique on doit simuler le chemin complet 
                auto path = simulateFullPath();
                auto asianOption = std::dynamic_pointer_cast<AsianOption>(option_);
                if (asianOption) {
                    asianOption->setPricePath(path);
                    payoff = asianOption->payoff(path.back());
                }
            } else {// option vanille on simule juste le prix final
                double finalPrice = simulatePath();
                payoff = option_->payoff(finalPrice);
            }
            
            sumPayoffs += payoff;
        }
        
        return discountFactor * sumPayoffs / numSimulations_; // prix moyen actualisé
    }
    
    std::pair<double, double> calculatePriceWithError() {
        double r = marketData_->getRiskFreeRate();
        double T = option_->getMaturity();
        double discountFactor = std::exp(-r * T);
        
        bool isPathDependent = (dynamic_cast<AsianOption*>(option_.get()) != nullptr);
        
        std::vector<double> payoffs;
        payoffs.reserve(numSimulations_);
        
        for (int i = 0; i < numSimulations_; ++i) {
            double payoff = 0.0;
            
            if (isPathDependent) {
                auto path = simulateFullPath();
                auto asianOption = std::dynamic_pointer_cast<AsianOption>(option_);
                if (asianOption) {
                    asianOption->setPricePath(path);
                    payoff = asianOption->payoff(path.back());
                }
            } else {
                double finalPrice = simulatePath();
                payoff = option_->payoff(finalPrice);
            }
            
            payoffs.push_back(payoff);
        }
        
        double mean = std::accumulate(payoffs.begin(), payoffs.end(), 0.0) / numSimulations_;
        
        double variance = 0.0;
        for (double payoff : payoffs) {
            variance += (payoff - mean) * (payoff - mean);
        }
        variance /= (numSimulations_ - 1);
        double stdDev = std::sqrt(variance);
        double standardError = stdDev / std::sqrt(numSimulations_); // erreur standard du prix estimé
        
        return {discountFactor * mean, discountFactor * standardError};
    }
    
    std::string getMethodName() const override { return "Monte Carlo (N=" + std::to_string(numSimulations_) + ")"; }
};

#endif
