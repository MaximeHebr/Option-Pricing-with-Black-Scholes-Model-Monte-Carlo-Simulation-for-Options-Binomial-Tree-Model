#ifndef MARKET_H
#define MARKET_H

#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>

// Market parameters
class MarketData {
private:
    double spotPrice_;
    double riskFreeRate_;
    double volatility_;
    double dividendYield_;
    
public:
    MarketData(double spotPrice, double riskFreeRate, double volatility, double dividendYield = 0.0)
        : spotPrice_(spotPrice), riskFreeRate_(riskFreeRate), 
          volatility_(volatility), dividendYield_(dividendYield) {
        if (spotPrice <= 0) throw std::invalid_argument("Spot price must be positive");
        if (volatility < 0) throw std::invalid_argument("Volatility cannot be negative");
    }
    
    double getSpotPrice() const { return spotPrice_; }
    double getRiskFreeRate() const { return riskFreeRate_; }
    double getVolatility() const { return volatility_; }
    double getDividendYield() const { return dividendYield_; }
    double getDrift() const { return riskFreeRate_ - dividendYield_ - 0.5 * volatility_ * volatility_; }
};

// Pricing results with Greeks
struct PricingResult {
    double price;
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
    
    PricingResult() : price(0), delta(0), gamma(0), vega(0), theta(0), rho(0) {}
    
    std::string toString() const {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(4);
        ss << "Price: " << price << "\n"
           << "Delta: " << delta << "\n"
           << "Gamma: " << gamma << "\n"
           << "Vega:  " << vega << "\n"
           << "Theta: " << theta << "\n"
           << "Rho:   " << rho;
        return ss.str();
    }
};

#endif
